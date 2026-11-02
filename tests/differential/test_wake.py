"""UVADD, CSCALC, XWINIT, SETXW: outputs and the complete compared State must match."""
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from harness import assert_close
from pyxrotor.wake import Velocities, cscalc, setxw, uvadd, xwinit
from rotors import perturbed_states
from statesync import assert_states_match, state_to_ref
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)


@SETTINGS
@given(s=perturbed_states(), frac=finite(-0.1, 1.3))
def test_uvadd(ref_module, s, frac):
    state_to_ref(s, ref_module)
    xiw = frac
    _, wa, wt = ref_module.call("UVADD", xiw, 0.0, 0.0)
    got = uvadd(s, xiw)
    assert_close(got[0], wa.item(), rtol=0, name="WA")
    assert_close(got[1], wt.item(), rtol=0, name="WT")


@SETTINGS
@given(s=perturbed_states(), data=st.data())
def test_cscalc(ref_module, s, data):
    i = data.draw(st.integers(0, s.ii - 1))
    state_to_ref(s, ref_module)
    out = ref_module.call("CSCALC", i + 1, *([0.0] * 22))
    got = cscalc(s, i)
    for name, o in zip(Velocities._fields, out[1:]):
        assert_close(getattr(got, name), o.item(), rtol=0, name=name)


@SETTINGS
@given(s=perturbed_states())
def test_xwinit(ref_module, s):
    state_to_ref(s, ref_module)
    ref_module.call("XWINIT")
    xwinit(s)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=perturbed_states())
def test_setxw(ref_module, s):
    state_to_ref(s, ref_module)
    ref_module.call("SETXW")
    setxw(s)
    assert_states_match(s, ref_module)
