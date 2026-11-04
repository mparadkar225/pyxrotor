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


def crotor_ref(f1, f2, rpm, rrat, vel, crit, ccon, w, vclr):
    """LFWD f1, LAFT f2, (VCLR), INPU, RPM: CROTOR's command sequence."""
    ref = Ref755()
    for icr, f in ((1, f1), (2, f2)):
        ref["FNAME"] = str(f)
        ref.call("REFCRLOAD", icr)
    if vclr:
        ref.call("REFCRVCLR")
    out = ref.call("REFCRRPM", rpm, rrat, vel, crit, ccon, *w,
                   0, np.zeros(NHIST), False, False, np.zeros(150), 0.0, 0.0)
    nsolv, hist, lcon, lok, gamf, ttf, ptf = out[9:16]
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
          w=(0.5, 0.0, 0.5, -1.0), vclr=False):
    f1, f2 = files(d)
    ref, nsolv, hist, lcon, lok, gamf, ttf, ptf = crotor_ref(
        f1, f2, rpm, rrat, vel, crit, ccon, w, vclr)
    p = pair(f1, f2, rpm, rrat, vel, crit, ccon, w, vclr)
    if not lok:
        with pytest.raises(PairConvergenceError, match=f"solve {nsolv}\\)"):
            p.solve_rpm(rpm)
        return None
    res = p.solve_rpm(rpm)
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
