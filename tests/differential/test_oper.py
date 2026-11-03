"""TPQ, APINIT, APITER, APER: the complete compared State must match after each call."""
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor.oper import apinit, aper, apiter, tpq
from rotors import perturbed_states
from statesync import assert_states_match, state_to_ref
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)
from conftest import EXPENSIVE  # noqa: E402


@SETTINGS
@given(s=perturbed_states(), itype=st.sampled_from([1, 2, 3]))
def test_tpq(ref_module, s, itype):
    state_to_ref(s, ref_module)
    ref_module.call("TPQ", itype)
    tpq(s, itype)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=perturbed_states(), itype=st.sampled_from([2, 3]), data=st.data())
def test_tpq_design_negative_circulation(ref_module, s, itype, data):
    """Design modes with some GAM <= 0: TPQ(2) borrows a neighbour's chord there."""
    ii = s.ii
    for i in data.draw(st.lists(st.integers(0, ii - 1), min_size=1, max_size=4, unique=True)):
        s.gam[i] = -abs(s.gam[i]) * data.draw(finite(0.0, 1.0))
    state_to_ref(s, ref_module)
    ref_module.call("TPQ", itype)
    tpq(s, itype)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=perturbed_states(), adv_scale=finite(0.5, 1.8))
def test_apinit(ref_module, s, adv_scale):
    s.adv *= adv_scale
    state_to_ref(s, ref_module)
    ref_module.call("APINIT")
    apinit(s)
    assert_states_match(s, ref_module)


@st.composite
def operating_cases(draw, allow_vortex=False):
    """A state plus (ISPEC, ICON) for APITER: targets near the current solution, options varied."""
    s = draw(perturbed_states())
    rng = np.random.default_rng(draw(st.integers(0, 2**32 - 1)))
    ii = s.ii
    s.gam[:ii] *= 1.0 + draw(finite(0.0, 0.3)) * rng.standard_normal(ii)
    s.beta[:ii] += draw(finite(-0.05, 0.05))
    s.beta0[:ii] = s.beta[:ii]
    s.free = draw(st.booleans())
    s.fast = draw(st.booleans())
    s.vrtx = allow_vortex and not s.fast and draw(st.booleans())
    s.nitera = draw(st.sampled_from([0, 1, 3, 40, 40, 40]))   # mostly full solves
    s.conv = False
    ispec = draw(st.sampled_from([1, 2, 3, 4, 5]))
    icon = draw(st.sampled_from([1, 2]))
    dim = s.rho * s.vel**2 * s.rad**2
    s.tspec = s.ttot * dim * draw(finite(0.7, 1.3))
    s.qspec = s.qtot * dim * s.rad * draw(finite(0.7, 1.3))
    s.pspec = s.ptot * dim * s.vel * draw(finite(0.7, 1.3))
    if ispec == 5 and draw(st.booleans()):
        # engine power-rpm line around the current rpm
        rpm = s.vel / (s.rad * s.adv * 3.141592654 / 30.0)
        n = draw(st.integers(3, 6))
        rpms = np.sort(rpm * rng.uniform(0.5, 1.5, n)) + np.arange(n)
        pwr = s.pspec * rng.uniform(0.6, 1.4, n)
        s.lpwrvar, s.npwrvar = True, n
        s.rpmvar[:n], s.pwrvar[:n] = rpms, pwr
        s.xpwrvar[:n] = rng.normal(0.0, abs(s.pspec) / max(rpm, 1.0), n)
    return s, ispec, icon


@EXPENSIVE
@given(case=operating_cases())
def test_apiter(ref_module, case):
    s, ispec, icon = case
    state_to_ref(s, ref_module)
    ref_module.call("APITER", ispec, icon)
    apiter(s, ispec, icon)
    assert_states_match(s, ref_module)


@EXPENSIVE
@given(case=operating_cases(), linit=st.booleans())
def test_aper(ref_module, case, linit):
    s, ispec, icon = case
    state_to_ref(s, ref_module)
    ref_module.call("APER", ispec, icon, linit)
    aper(s, ispec, icon, linit)
    assert_states_match(s, ref_module)


@settings(deadline=None, max_examples=3)
@given(case=operating_cases(allow_vortex=True))
def test_apiter_vortex_wake(ref_module, case):
    """With the discrete vortex wake (slow: a few examples only)."""
    s, ispec, icon = case
    s.vrtx, s.fast, s.duct = True, False, False
    s.nitera = min(s.nitera, 3)
    s.adw = max(s.adw, 0.15)
    state_to_ref(s, ref_module)
    ref_module.call("APITER", ispec, icon)
    apiter(s, ispec, icon)
    assert_states_match(s, ref_module)


@pytest.mark.parametrize("bad", ["all_chord", "one_chord", "beta", "adv_zero"])
def test_apinit_non_finite(ref, bad):
    """APINIT's MAX(VHSQ,-0.25) and MIN/MAX(DADW, k*ADW) return their *second* argument for a
    NaN first one (MAXSD/MINSD source operand; see tools/minmax_sites.py), e.g. after a
    diverged design leaves NaN geometry."""
    from rotors import solved_state
    s = solved_state()
    if bad == "all_chord":
        s.ch[:s.ii] = np.nan
    elif bad == "one_chord":
        s.ch[s.ii // 2] = np.nan
    elif bad == "beta":
        s.beta[:3] = np.nan
    else:
        s.adv = 0.0
    state_to_ref(s, ref)
    ref.call("APINIT")
    apinit(s)
    assert_states_match(s, ref)
