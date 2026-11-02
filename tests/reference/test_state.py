"""pyxrotor.State mirrors XROTOR.INC exactly, and round-trips through the reference."""
import pytest

from pyxrotor import _layout
from pyxrotor.state import State
from statesync import assert_states_match, state_from_ref, state_to_ref

pytestmark = pytest.mark.reference


def test_layout_is_current():
    """src/pyxrotor/_layout.py must be regenerated (tools/gen_layout.py) if XROTOR.INC changes."""
    import tools.gen_layout as gen
    assert (gen.OUT.read_text() == gen.render())


def test_state_covers_every_common_variable(ref):
    s = State()
    assert set(State.SPEC) == {v.name for b in ref.blocks.values() for v in b.vars}
    for name, (kind, bounds) in State.SPEC.items():
        value = getattr(s, name.lower())
        if bounds:
            assert value.shape == ref[name].shape, name


def test_round_trip_through_reference(ref):
    """Start-up state -> State -> a fresh reference instance -> identical."""
    from reference.ref755 import Ref755
    s = state_from_ref(ref)
    other = Ref755(start=False)
    state_to_ref(s, other)
    assert_states_match(state_from_ref(other), ref)
    assert s.ii == 30 and s.iaero[0] == -1          # INIT leaves IAERO = 0, i.e. unset
