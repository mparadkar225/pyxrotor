"""RotorPair(interaction="weights") against CROTOR's counter-rotation iteration.

The oracle is CROTOR's own code (crotor.f v0.85, CROTOR 755-ES1.3), transcribed into the
reference shim (REFCRLOAD, REFCRVCLR, REFCRRPM): both rotors share XROTOR's COMMON, and
LOADCR swaps the geometry and slipstream in before each solve. The port keeps two States.
The thrust-change history, the solve count and both rotors' results must agree bit for
bit, and the whole State of the rear rotor (the last solved) must match the COMMON."""
import shutil
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor import Rotor
from pyxrotor.pair import PairConvergenceError, RotorPair
from reference.ref755 import Ref755
from statesync import assert_states_match
from test_restart import fresh, short_dir  # noqa: F401  (fixture)

pytestmark = pytest.mark.differential
from conftest import EXPENSIVE  # noqa: E402

CROTOR = Path(__file__).resolve().parents[1] / "data" / "crotor"
NHIST = 200


def files(d):
    f1, f2 = fresh(d), fresh(d)
    shutil.copy(CROTOR / "cr1.xr", f1)
    shutil.copy(CROTOR / "cr2.xr", f2)
    return f1, f2


def crotor_ref(f1, f2, rpm, rrat, vel, crit, ccon, w, vclr, powe=None):
    """LFWD f1, LAFT f2, (VCLR), INPU, then RPM, or POWE when powe = (front power, rear
    power, FIX, SYNC): CROTOR's command sequence."""
    ref = Ref755()
    for icr, f in ((1, f1), (2, f2)):
        ref["FNAME"] = str(f)
        ref.call("REFCRLOAD", icr)
    if vclr:
        ref.call("REFCRVCLR")
    tail = (0, np.zeros(NHIST), False, False, np.zeros(150), 0.0, 0.0)
    if powe is None:
        out = ref.call("REFCRRPM", rpm, rrat, vel, crit, ccon, *w, *tail)[9:16]
    else:
        pf, pa, lpitch, lsrpm = powe
        out = ref.call("REFCRPOWE", pf, pa, rpm, rrat, vel, lpitch, lsrpm,
                       crit, ccon, *w, *tail)[13:20]
    nsolv, hist, lcon, lok, gamf, ttf, ptf = out
    return ref, int(nsolv), list(hist[:int(nsolv)]), bool(lcon), bool(lok), gamf, ttf, ptf


def pair(f1, f2, rpm, rrat, vel, crit, ccon, w, vclr):
    front, rear = Rotor(), Rotor()
    front.load(f1)
    rear.load(f2)
    for r in (front, rear):
        r.s.vel = 0.0                       # LOADF/LOADA: VEL = VELCR, zero before INPU
        if vclr:                            # VCLR zeroes the counts only
            r.s.nadd = 0
            r.s.conv = False
    p = RotorPair(front, rear, 0.0, gear_ratio=rrat, interaction="weights", weights=w,
                  thrust_tol=ccon, max_solves=crit)
    p.velocity(vel)
    return p


def check(d, rpm=3000.0, rrat=1.0, vel=35.0, crit=20, ccon=1e-3,
          w=(0.5, 0.0, 0.5, -1.0), vclr=False, power=None, split=0.5, fix="rpm",
          locked=True):
    """RPM (power None) or POWE (total power, front share split, FIX = fix == "pitch",
    SYNC = locked), CROTOR against the port."""
    f1, f2 = files(d)
    powe = None
    if power is not None:
        powe = (split * power, (1.0 - split) * power, fix == "pitch", locked)
    ref, nsolv, hist, lcon, lok, gamf, ttf, ptf = crotor_ref(
        f1, f2, rpm, rrat, vel, crit, ccon, w, vclr, powe)
    p = pair(f1, f2, rpm, rrat, vel, crit, ccon, w, vclr)

    def run():
        if power is None:
            return p.solve_rpm(rpm)
        return p.solve_power_split(power, rpm, split, fix=fix, locked=locked)
    if not lok:
        with pytest.raises(PairConvergenceError, match=f"solve {nsolv}\\)"):
            run()
        failed = p.front if nsolv % 2 else p.rear       # blade angles restored (label 590)
        np.testing.assert_array_equal(failed.s.beta[:failed.s.ii], ref["BETA"][:failed.s.ii])
        return None
    res = run()
    assert res.iterations == nsolv
    assert res.converged == lcon
    assert res.history == hist
    s = p.front.s
    assert s.ttot == ttf and s.ptot == ptf
    np.testing.assert_array_equal(s.gam[:s.ii], gamf[:s.ii])
    assert_states_match(p.rear.s, ref)
    return res


def test_crotor_sample(short_dir):  # noqa: F811
    """CROTOR's README case: 80 kW design point, 3000 rpm, 35 m/s."""
    res = check(short_dir)
    assert res.converged and res.iterations == 6
    assert abs(res.power / 80e3 - 1.0) < 1e-3              # its design point
    assert abs(res.front["thrust"] - 844.5) < 0.1 and abs(res.rear["thrust"] - 864.9) < 0.1


def test_crotor_from_cleared_slipstreams(short_dir):  # noqa: F811
    """VCLR first: the front's first solve is the isolated rotor."""
    assert check(short_dir, vclr=True).converged


def test_crotor_iteration_limit(short_dir):  # noqa: F811
    """CRIT reached: CROTOR still finishes on a front-then-rear pair, unconverged."""
    res = check(short_dir, crit=3, ccon=1e-9)
    assert not res.converged and res.iterations == 4


@EXPENSIVE
@given(rpm=st.floats(2200.0, 4000.0), rrat=st.floats(0.85, 1.15), vel=st.floats(20.0, 50.0),
       crit=st.integers(1, 30), ccon=st.sampled_from([1e-2, 1e-3, 1e-5, 1e-9]),
       uf=st.floats(0.0, 1.0), vf=st.sampled_from([0.0, 1.0, -1.0]),
       ua=st.floats(0.0, 1.0), va=st.sampled_from([-1.0, 0.0, 1.0]), vclr=st.booleans())
def test_crotor_rpm(short_dir, rpm, rrat, vel, crit, ccon, uf, vf, ua, va, vclr):  # noqa: F811
    check(short_dir, rpm, rrat, vel, crit, ccon, (uf, vf, ua, va), vclr)


def test_crotor_power_split_fixed_rpm(short_dir):  # noqa: F811
    """POWE with the default FIX (fixed rpm, variable pitch), README case: 40 kW each."""
    res = check(short_dir, power=80e3)
    assert res.converged
    assert abs(res.front["power"] - 40e3) < 1e-6 * 40e3
    assert abs(res.rear["power"] - 40e3) < 1e-6 * 40e3


def test_crotor_power_free_rotors(short_dir):  # noqa: F811
    """FIX on, SYNC off: fixed pitch, each rotor finds its own rpm for 40 kW."""
    res = check(short_dir, power=80e3, fix="pitch", locked=False)
    assert res.converged and abs(res.power - 80e3) < 1e-6 * 80e3


def test_crotor_power_locked_rpm(short_dir):  # noqa: F811
    """FIX on, SYNC on: CROTOR's two-phase procedure ends at one rpm, near the power."""
    res = check(short_dir, power=80e3, fix="pitch", locked=True, rrat=1.0)
    assert res.converged
    assert res.front["rpm"] == pytest.approx(res.rear["rpm"], rel=1e-12)
    assert abs(res.power / 80e3 - 1.0) < 0.02


@EXPENSIVE
@given(rpm=st.floats(2400.0, 3600.0), rrat=st.floats(0.9, 1.1), vel=st.floats(25.0, 45.0),
       crit=st.integers(1, 30), ccon=st.sampled_from([1e-2, 1e-3, 1e-5, 1e-9]),
       power=st.floats(30e3, 120e3), split=st.floats(0.3, 0.7),
       fix=st.sampled_from(["rpm", "pitch"]), locked=st.booleans(),
       ua=st.floats(0.0, 1.0), va=st.sampled_from([-1.0, 0.0, 1.0]), vclr=st.booleans())
def test_crotor_powe(short_dir, rpm, rrat, vel, crit, ccon, power, split, fix, locked,  # noqa: F811
                     ua, va, vclr):
    check(short_dir, rpm, rrat, vel, crit, ccon, (0.5, 0.0, ua, va), vclr,
          power=power, split=split, fix=fix, locked=locked)
