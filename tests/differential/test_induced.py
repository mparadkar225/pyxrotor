"""VCALC, GRADMO, HELICO, VORSEGVEL, VRTXCO against the reference (complete State)."""
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from harness import assert_close
from pyxrotor.induced import gradmo, helico, vcalc, vorsegvel, vrtxco
from rotors import perturbed_states, solved_state
from statesync import assert_states_match, state_to_ref
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)


@st.composite
def wake_inputs(draw):
    """A converged state with perturbed GAM, ADW, blade count and duct flag."""
    s = draw(perturbed_states(added_inflow=False))
    rng = np.random.default_rng(draw(st.integers(0, 2**32 - 1)))
    s.gam[:s.ii] *= 1.0 + draw(finite(0.0, 0.5)) * rng.standard_normal(s.ii)
    s.nblds = draw(st.integers(1, 6))
    s.adw = abs(s.adw) * draw(finite(0.5, 2.0))
    return s


@SETTINGS
@given(s=wake_inputs())
def test_vcalc(ref_module, s):
    state_to_ref(s, ref_module)
    ref_module.call("VCALC")
    vcalc(s)
    assert_states_match(s, ref_module)


def call_model(ref, name, s):
    ref.call(name, 100, s.ii, s.nblds, s.duct, s.rake, ref["XI"], ref["XV"], ref["GAM"], s.adw,
             ref["VIND_GAM"], ref["VIND_ADW"])


@SETTINGS
@given(s=wake_inputs())
def test_gradmo(ref_module, s):
    state_to_ref(s, ref_module)
    call_model(ref_module, "GRADMO", s)
    gradmo(s)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=wake_inputs())
def test_helico(ref_module, s):
    state_to_ref(s, ref_module)
    call_model(ref_module, "HELICO", s)
    helico(s)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=wake_inputs(), stale=finite(-1.0, 1.0))
def test_helico_reads_gam_beyond_ii(ref_module, s, stale):
    """HELICO reads GAM(I) past the blade stations (I up to IIMAX-1); those values matter."""
    s.gam[s.ii:] = stale
    state_to_ref(s, ref_module)
    call_model(ref_module, "HELICO", s)
    helico(s)
    assert_states_match(s, ref_module)


@SETTINGS
@given(m=st.integers(1, 20), data=st.data())
def test_vorsegvel(ref_module, m, data):
    pts = st.one_of(finite(-3.0, 3.0), st.just(0.0))
    a = data.draw(arrays(np.float64, (3, m), elements=pts))
    b = data.draw(arrays(np.float64, (3, m), elements=pts))
    uvw, ua, ub = vorsegvel(a, b, full=True)
    for k in range(m):
        out = ref_module.call("VORSEGVEL", a[:, k].copy(), b[:, k].copy(), np.zeros(3),
                              np.zeros((3, 3), order="F"), np.zeros((3, 3), order="F"))
        assert_close(uvw[:, k], out[2], rtol=0, name="UVW")
        assert_close(ua[:, :, k], out[3], rtol=0, name="UVW_A")
        assert_close(ub[:, :, k], out[4], rtol=0, name="UVW_B")


@settings(deadline=None, max_examples=6)
@given(s=wake_inputs(), rake=st.sampled_from([0.0, 0.1]))
def test_vrtxco(ref_module, s, rake):
    """Expensive (about 10^6 segments), so few examples even under the deep profile."""
    s.rake = rake
    s.adw = max(s.adw, 0.15)                 # keep the helix (50/ADW) a manageable length
    state_to_ref(s, ref_module)
    call_model(ref_module, "VRTXCO", s)
    vrtxco(s)
    assert_states_match(s, ref_module)
