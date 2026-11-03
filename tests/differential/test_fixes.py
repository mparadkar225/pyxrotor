"""Bug-fix options: each against the reference with the fix applied there (patched shim or
prepared inputs), and each shown to change what 7.55 does.

    Fixes.linear_cl         REFDESGENX(..., LFIXCL=.TRUE., CL1, CL2)
    Fixes.section_fallback  unmodified 7.55 routines with the invalid IAERO entries filled in
                            as SETIAERO would (the fixed fallback does exactly that, then
                            proceeds as for a valid index)
    Fixes.design_sections   REFDESGEN(..., LFIXIA=.TRUE.)   (tests in test_design.py)
"""
import numpy as np
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from pyxrotor.aero import AeroSection, getalf, getclcdcm, setiaero
from pyxrotor.constants import PI
from pyxrotor.design import desgen
from pyxrotor.oper import aper, sections, tpq
from pyxrotor.rotor import Rotor
from rotors import perturbed_states, solved_state
from statesync import assert_states_match, state_to_ref
from strategies import alphas, finite, flight, iaero_to_fortran, section_sets, speeds
from test_aero import call_getclcdcm, same, set_stations
from test_design import design_points

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)
from conftest import EXPENSIVE  # noqa: E402


# ------------------------------------------------------------------- linear design CL
@EXPENSIVE
@given(s=design_points(clean=True), cl1=finite(0.2, 1.2), cl2=finite(0.1, 1.2),
       fixia=st.booleans())
def test_linear_cl_fix(ref_module, s, cl1, cl2, fixia):
    """DESGEN with the linear CL re-evaluated on the design grid (and the stall check there),
    including tip CLs past section 1's limits."""
    s.cldes_linear = (cl1, cl2)
    state_to_ref(s, ref_module)
    ref_module.call("REFDESGENX", True, fixia, True, cl1, cl2)
    desgen(s, betz_ok=True, refresh_sections=fixia, fix_linear_cl=True)
    assert_states_match(s, ref_module)


def rotol_sections():
    from test_design import ROTOL_SECTIONS
    return ROTOL_SECTIONS


def test_linear_cl_first_design():
    """A first design with a linear CL: 7.55 gives every station the root CL (all XI are 0
    when SETCLD runs); the fix gives the requested distribution on the design's stations."""
    out = {}
    for fix in (False, True):
        r = Rotor(fix_linear_cl=fix)
        r.atmosphere(0.0)
        r.set_sections(rotol_sections())
        r.s.raddes, r.s.r0des, r.s.rwdes, r.s.veldes, r.s.nblds = 1.6, 0.25, 0.0, 60.0, 3
        r.s.rpmdes, r.s.advdes, r.s.dest, r.s.desp, r.s.pddes = 1200.0, 0.0, False, True, 6.0e5
        r.s.name = "linear CL"
        assert r.redesign(cl_root_tip=(0.7, 0.4))
        xi, cl = r.s.xi[: r.s.ii], r.s.cldes[: r.s.ii]          # the design CL
        out[fix] = r.performance()["efficiency"]
        if fix:
            np.testing.assert_allclose(cl, 0.7 * (1 - xi) + 0.4 * xi, rtol=0, atol=1e-15)
        else:
            assert (cl == 0.7).all()
    assert out[True] != out[False]


def test_linear_cl_after_load_is_not_reapplied(tmp_path):
    """LOAD sets CLDES from the solved CL: the fix must not re-impose an older linear CL."""
    r = Rotor(fix_linear_cl=True)
    r.s.cldes_linear = (0.9, 0.1)
    from pyxrotor.restart import load
    from pathlib import Path
    load(r.s, Path(__file__).resolve().parents[1] / "data" / "p51d_prop.xr")
    assert r.s.cldes_linear is None


# ------------------------------------------------------------------- section fallback
def fill_like_setiaero(iaero, xi, secs, stations):
    out = iaero.copy()
    ref_idx = setiaero(xi, secs)
    for i in stations:
        out[i] = ref_idx[i]
    return out


@SETTINGS
@given(secs=section_sets(), alf=alphas, w=speeds, fl=flight, bad=st.sampled_from([-1, 5, 20]),
       x_over=st.booleans(), data=st.data())
def test_getclcdcm_fixed_fallback(ref_module, secs, alf, w, fl, bad, x_over, data):
    """The fixed fallback equals 7.55 with that station's IAERO set as SETIAERO does,
    including stations outboard of the last section (7.55 falls back to section 1 there)."""
    assume(bad < 0 or bad >= len(secs))
    ref = ref_module
    ii = 30
    xi = set_stations(ref, secs, ii, 0.15)
    ref["VEL"], ref["VSO"] = fl
    is_ = data.draw(st.integers(0, ii - 1))
    iaero = setiaero(xi, secs)
    iaero[is_] = bad
    fixed = fill_like_setiaero(iaero, xi, secs, [is_])
    ref["IAERO"][:ii] = iaero_to_fortran(fixed)
    want = call_getclcdcm(ref, is_, alf, w, 0.0)
    got_iaero = iaero.copy()
    got = getclcdcm(is_, alf, w, 0.0, xi=xi, iaero=got_iaero, sections=secs, vel=fl[0],
                    vso=fl[1], fix_fallback=True)
    for name, value in want.items():
        same(getattr(got, name), value, name)
    np.testing.assert_array_equal(got_iaero, fixed)


def test_fallback_7_55_vs_fixed_cases():
    """The two 7.55 behaviors the fix changes, on three sections at r/R 0, 0.5, 0.8."""
    secs = [AeroSection(xisect=0.0, cdmin=0.02), AeroSection(xisect=0.5, cdmin=0.01),
            AeroSection(xisect=0.8, cdmin=0.006)]
    xi = np.array([0.6, 0.9])
    for i, want_755, want_fix in [(0, 1, 1), (1, 0, 2)]:
        for fix, want in ((False, want_755), (True, want_fix)):
            ia = np.array([-1, -1])
            getclcdcm(i, 0.05, 1.0, 0.0, xi=xi, iaero=ia, sections=secs, vel=50.0, vso=340.0,
                      fix_fallback=fix)
            assert ia[i] == want, (i, fix)
    # between sections the stored index is the inner one in both, but 7.55 *used* the outer
    # one (section 3 at r/R 0.8, extrapolating) for the first lookup:
    a755 = getclcdcm(0, 0.05, 1.0, 0.0, xi=xi, iaero=np.array([-1, -1]), sections=secs, vel=50.0, vso=340.0)
    afix = getclcdcm(0, 0.05, 1.0, 0.0, xi=xi, iaero=np.array([-1, -1]), sections=secs, vel=50.0, vso=340.0,
                     fix_fallback=True)
    assert a755.cdrag < afix.cdrag      # 7.55: section 3 alone; fixed: 2-3 interpolated


@EXPENSIVE
@given(s=perturbed_states(duct=False, added_inflow=False), data=st.data())
def test_section_fallback_in_solver(ref_module, s, data):
    """TPQ and a full APER solve with some stations' section index invalid (e.g. after
    sections were deleted): fixed port == 7.55 with those indices filled in beforehand."""
    ii = s.ii
    secs = sections(s)
    bad_st = data.draw(st.lists(st.integers(0, ii - 1), min_size=1, max_size=ii, unique=True))
    for i in bad_st:
        s.iaero[i] = data.draw(st.sampled_from([-1, len(secs), len(secs) + 3]))
    fixed = s.copy()
    fixed.iaero[:ii] = fill_like_setiaero(s.iaero[:ii], s.xi[:ii], secs, bad_st)
    s.fixes.section_fallback = True
    state_to_ref(fixed, ref_module)
    if data.draw(st.booleans()):
        ref_module.call("TPQ", 1)
        tpq(s, 1)
    else:
        s.nitera = fixed.nitera = 40
        ref_module.call("APER", 4, 2, True)
        aper(s, 4, 2, True)
    assert_states_match(s, ref_module)


def test_fix_all_flag():
    r = Rotor(fix_all=True)
    f = r.s.fixes
    assert f.design_sections and f.linear_cl and f.section_fallback
    assert Rotor().s.fixes == type(f)()           # all off by default
