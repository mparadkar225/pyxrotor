"""VPUT (SAVVEL) and VGET (GETVEL) against the reference; LOAD's slipstream spline fix."""
import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from pyxrotor.fortio import write_list
from pyxrotor.restart import load, save
from pyxrotor.rotor import Rotor
from pyxrotor.slipstream import getvel, profile, savvel
from reference.ref755 import Ref755
from rotors import perturbed_states, solved_state
from statesync import assert_states_match, state_to_ref
from test_rotor import rotol_both
from strategies import finite
from test_restart import fresh, short_dir  # noqa: F401  (fixture)

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
from conftest import EXPENSIVE  # noqa: E402


@SETTINGS
@given(s=perturbed_states())
def test_vput(ref_module, short_dir, s):  # noqa: F811
    """The same bytes as SAVVEL (ducts, nacelle and imposed inflow included), same State."""
    s.conv = True
    fr, fp = fresh(short_dir, ".vel"), fresh(short_dir, ".vel")
    state_to_ref(s, ref_module)
    ref_module["FNAME"] = str(fr)
    ref_module.call("REFSAVVEL")
    savvel(s, fp)
    assert fp.read_bytes() == fr.read_bytes()
    assert_states_match(s, ref_module)


def test_vput_needs_solution(short_dir):  # noqa: F811
    s = solved_state()
    s.conv = False
    f = fresh(short_dir, ".vel")
    ref = Ref755()
    state_to_ref(s, ref)
    ref["FNAME"] = str(f)
    ref.call("REFSAVVEL")
    assert not f.exists()
    with pytest.raises(ValueError):
        savvel(s, f)


@st.composite
def profiles(draw):
    n = draw(st.integers(0, 110))
    r = np.sort(np.array(draw(st.lists(finite(0.05, 2.0), min_size=n, max_size=n, unique=True))))
    u = np.array(draw(st.lists(finite(-5.0, 30.0), min_size=n, max_size=n)))
    v = np.array(draw(st.lists(finite(-20.0, 20.0), min_size=n, max_size=n)))
    style = draw(st.sampled_from(["list", "split", "extra", "blank"]))
    return r, u, v, style


def write_profile(path, r, u, v, style):
    r, u, v = ([float(x) for x in a] for a in (r, u, v))
    lines = []
    for k in range(len(r)):
        if style == "split":
            lines += [repr(r[k]), f"{u[k]!r} {v[k]!r}"]
        elif style == "extra":
            lines.append(f"{r[k]!r} {u[k]!r} {v[k]!r} 99.0")
        else:
            lines.append(write_list(r[k], u[k], v[k]))
            if style == "blank":
                lines.append("")
    path.write_text("\n".join(lines) + ("\n" if lines else ""))


@SETTINGS
@given(s=perturbed_states(), prof=profiles(), uwt=st.sampled_from([0.0, 0.5, 1.0, 0.73]),
       vwt=st.sampled_from([0.0, 1.0, -1.0]))
def test_vget(ref_module, short_dir, s, prof, uwt, vwt):  # noqa: F811
    """GETVEL: weights, splines, more than IX records, fewer than 2 (not read)."""
    f = fresh(short_dir, ".vel")
    write_profile(f, *prof)
    state_to_ref(s, ref_module)
    ref_module["FNAME"] = str(f)
    ref_module.call("REFGETVEL", uwt, vwt)
    ok = getvel(s, f, uwt, vwt)
    assert ok == (min(len(prof[0]), 100) >= 2)
    assert_states_match(s, ref_module)


@pytest.mark.parametrize("text", ["0.1 1.0 0.0\n0.2 x 0.0\n0.3 1.0 0.0\n", "", "0.1 2.0\n"])
def test_vget_bad_files(short_dir, text):  # noqa: F811
    f = fresh(short_dir, ".vel")
    f.write_text(text)
    s = solved_state()
    ref = Ref755()
    state_to_ref(s, ref)
    ref["FNAME"] = str(f)
    ref.call("REFGETVEL", 1.0, 1.0)
    assert not getvel(s, f, 1.0, 1.0)
    assert_states_match(s, ref)


@EXPENSIVE
@given(theta=st.sampled_from([15.0, 22.0, 27.0]), uwt=finite(0.3, 1.0), vwt=st.sampled_from([1.0, -1.0]),
       rpm=finite(800.0, 1400.0))
def test_two_rotor_coupling_by_hand(short_dir, theta, uwt, vwt, rpm):  # noqa: F811
    """The XROTOR procedure for a rotor behind another: solve the front one, VPUT, VGET it
    into the rear one (here the same blade), solve. Both steps match the reference."""
    ref, front = rotol_both(theta)
    ref.call("REFRPM", rpm)
    front.solve_rpm(rpm)
    assert_states_match(front.s, ref)
    f = fresh(short_dir, ".vel")
    ref["FNAME"] = str(f)
    ref.call("REFSAVVEL")
    if not front.s.conv:                                   # SAVVEL refuses: nothing written
        assert not f.exists()
        with pytest.raises(ValueError):
            front.save_slipstream(fresh(short_dir, ".vel"))
        return
    ref.call("REFGETVEL", uwt, vwt)
    front.save_slipstream(fresh(short_dir, ".vel"))       # same text, checked in test_vput
    p = front.slipstream_profile()
    front.set_slipstream(p["r"], p["u"], p["v"], uwt, vwt)
    ref.call("REFRPM", rpm)
    front.solve_rpm(rpm)
    assert_states_match(front.s, ref)


def test_load_slipstream_splines(short_dir):  # noqa: F811
    """7.55 LOAD leaves the slipstream splines stale; the fix splines them, which makes the
    loaded rotor behave as it did when the profile was imposed with VGET."""
    s = solved_state()
    r = np.linspace(0.1, s.rad, 8)
    from pyxrotor.slipstream import set_profile
    set_profile(s, r, 5.0 + 3.0 * np.sin(r), 2.0 * r, 1.0, 1.0)
    from pyxrotor.oper import aper
    aper(s, 4, 2, True)
    f = fresh(short_dir)
    save(s, f)
    res = {}
    for fix in (False, True):
        ref, rotor = Ref755(), Rotor(fix_slipstream_splines=fix)
        ref["FNAME"] = str(f)
        ref.call("REFLOAD")
        load(rotor.s, f)
        if not fix:
            assert_states_match(rotor.s, ref)                # 7.55: stale (zero) slopes
            assert (rotor.s.uaddr[: rotor.s.nadd] == 0.0).all()
        else:
            # the file holds 5 significant digits, so the slopes agree to about that
            np.testing.assert_allclose(rotor.s.uaddr[:8], s.uaddr[:8], rtol=2e-3, atol=1e-3)
        res[fix] = rotor.performance()["thrust"]
    assert res[True] != res[False]
