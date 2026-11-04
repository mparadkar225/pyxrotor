"""Section aerodynamics (xaero.f): CLCDCM, GETCLCDCM, GETALF, SETIAERO, SORTAR.

Every output, including the CLMAX/CLMIN that CLCDCM writes back and the IAERO that
GETCLCDCM can rewrite, must match XROTOR 7.55 bit for bit.
"""
import math

import numpy as np
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from harness import assert_close
from pyxrotor.aero import AeroSection, clcdcm, getalf, getclcdcm, setiaero, sort_sections
from strategies import (alphas, finite, flight, iaero_from_fortran, iaero_to_fortran,
                        load_sections, reynolds, section_sets, sections_params, speeds)

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)        # example count comes from the Hypothesis profile

CLCDCM_OUT = ["clift", "cl_alf", "cl_w", "stallf", "cdrag", "cd_alf", "cd_w", "cd_rey",
              "cmom", "cm_al", "cm_w"]
STATION_OUT = ["clift", "cl_alf", "cl_w", "clmax", "clmin", "dcl_stall", "stallf",
               "cdrag", "cd_alf", "cd_w", "cd_rey", "cmom", "cm_al", "cm_w"]


def same(got, want, name):
    """Bitwise agreement, treating NaN == NaN (the Fortran produces NaN in the same places)."""
    if isinstance(want, (bool, np.bool_)) or isinstance(got, bool):
        assert bool(got) == bool(want), name
        return
    got, want = float(got), float(want)
    if math.isnan(want):
        assert math.isnan(got), f"{name}: got {got!r}, want nan"
    else:
        assert_close(got, want, rtol=0, name=name)


# ------------------------------------------------------------------------------- CLCDCM
@SETTINGS
@given(alf=alphas, w=speeds, rey=reynolds, fl=flight, sec=sections_params())
def test_clcdcm(ref_module, alf, w, rey, fl, sec):
    ref = ref_module
    vel, vso = fl
    ref["VEL"], ref["VSO"] = vel, vso
    out = ref.call("CLCDCM", alf, w, rey,
                   0.0, 0.0, 0.0, False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   sec.a0, sec.clmax, sec.clmin, sec.dclda, sec.dclda_stall, sec.dcl_stall,
                   sec.cdmin, sec.cldmin, sec.dcdcl2, sec.cmcon, sec.mcrit, sec.reref, sec.rexp)
    want = dict(zip(CLCDCM_OUT, [o.item() for o in out[3:14]]))
    want["clmax"], want["clmin"] = out[15].item(), out[16].item()

    got = clcdcm(alf, w, rey, vel, vso, sec)
    for name, value in want.items():
        same(getattr(got, name), value, name)


NON_FINITE = st.sampled_from([math.nan, math.inf, -math.inf, 1e300, -1e300])


@SETTINGS
@given(alf=st.one_of(alphas, NON_FINITE), w=st.one_of(speeds, NON_FINITE), rey=reynolds,
       fl=flight, sec=sections_params(), bad_field=st.sampled_from([None, "mcrit", "clmax", "clmin",
       "a0", "dcl_stall", "cldmin"]), bad=NON_FINITE)
def test_clcdcm_non_finite(ref_module, alf, w, rey, fl, sec, bad_field, bad):
    """Diverged-iteration inputs: pins MAX/MIN NaN handling at every call site, and inf
    arithmetic, inside the compiled CLCDCM itself (a standalone probe behaved differently)."""
    if bad_field:
        setattr(sec, bad_field, bad)
    test_clcdcm.hypothesis.inner_test(ref_module, alf, w, rey, fl, sec)


@pytest.mark.parametrize("rey", [5e-324, 1e-320])
def test_clcdcm_reynolds_underflow(ref_module, rey):
    """REY/REREF underflows to 0 with a negative REXP: the Re correction pow(0, REXP) is a
    pole, +inf in the Fortran (fmath.pow returned nan until the compiled CLCDCM's
    comparison with the Python found it)."""
    sec = AeroSection(a0=0.0, clmax=1.0, clmin=0.0, dclda=2.0, dclda_stall=0.0, dcl_stall=0.0,
                      cdmin=0.03125, cldmin=0.0, dcdcl2=0.0, cmcon=0.0, reref=1e4, rexp=-1.0,
                      mcrit=1.0)
    test_clcdcm.hypothesis.inner_test(ref_module, 0.0, 0.0, rey, (120.0, 340.0), sec)


def test_clcdcm_physics():
    """Sanity checks independent of the Fortran: linear range, stall, Mach limiting."""
    sec = AeroSection()
    low = clcdcm(0.05, 0.3, 0.0, 100.0, 340.0, sec)          # Mach ~0.09
    # Linear Prandtl-Glauert lift, to within the smooth stall limiter's small effect
    # (about 1e-4 here, from the exponential tail of the negative-stall term).
    assert low.clift == pytest.approx(sec.dclda * 0.05 / math.sqrt(1 - (30 / 340) ** 2), rel=2e-4)
    assert not low.stallf
    assert low.cdrag == pytest.approx(sec.cdmin + sec.dcdcl2 * (low.clift - sec.cldmin) ** 2, rel=1e-3)

    stalled = clcdcm(0.6, 0.3, 0.0, 100.0, 340.0, sec)
    assert stalled.stallf and stalled.clift < sec.clmax + 0.5 and stalled.cdrag > 10 * low.cdrag

    sonic = clcdcm(0.05, 1.2, 0.0, 300.0, 340.0, sec)        # local Mach ~1.06
    assert sonic.mach_limited and not low.mach_limited


# ------------------------------------------------------------------------------ SETIAERO
@SETTINGS
@given(secs=section_sets(), ii=st.integers(3, 99), xi0=finite(0.05, 0.5))
def test_setiaero(ref_module, secs, ii, xi0):
    ref = ref_module
    load_sections(ref, secs)
    ref["II"], ref["XI0"], ref["XITIP"], ref["IXSPAC"] = ii, xi0, 1.0, 2
    ref.call("SETX")
    ref.call("SETIAERO")
    got = setiaero(ref.get("XI")[:ii], secs)
    assert np.array_equal(got, iaero_from_fortran(ref.get("IAERO")[:ii]))


# -------------------------------------------------------------------------------- SORTAR
@SETTINGS
@given(xs=st.lists(st.sampled_from([0.0, 0.2, 0.5, 0.7, 1.0]), min_size=1, max_size=8))
def test_sort_sections(ref_module, xs):
    """Includes duplicate radii, to pin down SORTAR's tie order."""
    ref = ref_module
    secs = [AeroSection(xisect=x, cdmin=float(k)) for k, x in enumerate(xs)]   # cdmin tags identity
    load_sections(ref, secs)
    n = len(secs)
    s, w = ref.call("SORTAR", n, np.array([sec.xisect for sec in secs]),
                    np.asfortranarray(ref.get("AERODATA")[:, :n]), 14)[1:3]
    got = sort_sections(secs)
    assert [sec.xisect for sec in got] == list(s)
    assert [sec.cdmin for sec in got] == list(w[6])


# ----------------------------------------------------------------------------- GETCLCDCM
def set_stations(ref, secs, ii, xi0):
    load_sections(ref, secs)
    ref["II"], ref["XI0"], ref["XITIP"], ref["IXSPAC"] = ii, xi0, 1.0, 2
    ref.call("SETX")
    return ref.get("XI")[:ii]


def call_getclcdcm(ref, is_, alf, w, rey):
    out = ref.call("GETCLCDCM", is_ + 1, alf, w, rey,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return dict(zip(STATION_OUT, [o.item() for o in out[4:18]]))


@SETTINGS
@given(secs=section_sets(), alf=alphas, w=speeds, rey=reynolds, fl=flight, data=st.data())
def test_getclcdcm(ref_module, secs, alf, w, rey, fl, data):
    """With IAERO as SETIAERO leaves it, and with stale indices pointing at any section."""
    ref = ref_module
    ii = 30
    xi = set_stations(ref, secs, ii, 0.15)
    ref["VEL"], ref["VSO"] = fl
    is_ = data.draw(st.integers(0, ii - 1))
    if data.draw(st.booleans()):
        iaero = setiaero(xi, secs)
    else:                                                  # stale but valid indices
        iaero = np.array(data.draw(st.lists(st.integers(0, len(secs) - 1), min_size=ii, max_size=ii)))
    ref["IAERO"][:ii] = iaero_to_fortran(iaero)

    want = call_getclcdcm(ref, is_, alf, w, rey)
    got = getclcdcm(is_, alf, w, rey, xi=xi, iaero=iaero.copy(), sections=secs, vel=fl[0], vso=fl[1])
    for name, value in want.items():
        same(getattr(got, name), value, name)


@SETTINGS
@given(secs=section_sets(), alf=alphas, w=speeds, fl=flight, bad=st.sampled_from([-1, 5, 20]),
       data=st.data())
def test_getclcdcm_invalid_index(ref_module, secs, alf, w, fl, bad, data):
    """The fallback search for an invalid IAERO (e.g. INIT's 0), including what it writes back."""
    ref = ref_module
    ii = 30
    xi = set_stations(ref, secs, ii, 0.15)
    ref["VEL"], ref["VSO"] = fl
    is_ = data.draw(st.integers(0, ii - 1))
    assume(bad < 0 or bad >= len(secs))
    iaero = setiaero(xi, secs)
    iaero[is_] = bad
    ref["IAERO"][:ii] = iaero_to_fortran(iaero)

    want = call_getclcdcm(ref, is_, alf, w, 0.0)
    got_iaero = iaero.copy()
    got = getclcdcm(is_, alf, w, 0.0, xi=xi, iaero=got_iaero, sections=secs, vel=fl[0], vso=fl[1])
    for name, value in want.items():
        same(getattr(got, name), value, name)
    assert got_iaero[is_] == iaero_from_fortran(ref.get("IAERO")[is_])


# -------------------------------------------------------------------------------- GETALF
@SETTINGS
@given(secs=section_sets(), cl=finite(-1.8, 2.2), w=finite(0.2, 4.0), fl=flight, data=st.data())
def test_getalf(ref_module, secs, cl, w, fl, data):
    """Both converged inversions and targets beyond stall, where the Newton iteration fails."""
    ref = ref_module
    ii = 30
    xi = set_stations(ref, secs, ii, 0.15)
    ref["VEL"], ref["VSO"] = fl
    iaero = setiaero(xi, secs)
    ref["IAERO"][:ii] = iaero_to_fortran(iaero)
    is_ = data.draw(st.integers(0, ii - 1))

    out = ref.call("GETALF", is_ + 1, cl, w, 0.0, 0.0, 0.0, False)
    got = getalf(is_, cl, w, xi=xi, iaero=iaero.copy(), sections=secs, vel=fl[0], vso=fl[1])
    for name, k in (("alf", 3), ("alf_cl", 4), ("alf_w", 5), ("stallf", 6)):
        same(getattr(got, name), out[k].item(), name)


def test_getalf_reports_failure():
    secs = [AeroSection()]
    xi = np.linspace(0.2, 1.0, 10)
    iaero = setiaero(xi, secs)
    ok = getalf(5, 0.5, 1.0, xi=xi, iaero=iaero, sections=secs, vel=50.0, vso=340.0)
    assert ok.converged and abs(clcdcm(ok.alf, 1.0, 0.0, 50.0, 340.0, secs[0]).clift - 0.5) < 1e-6
    bad = getalf(5, 3.0, 1.0, xi=xi, iaero=iaero, sections=secs, vel=50.0, vso=340.0)
    assert not bad.converged
