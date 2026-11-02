"""SETX: the ported radial grid must match the reference bit for bit."""
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from harness import assert_close
from pyxrotor.geometry import setx

pytestmark = pytest.mark.differential


@settings(max_examples=300, deadline=None)
@given(
    ii=st.integers(min_value=3, max_value=99),          # XV(II+1) must fit in XV(IX)
    xi0=st.floats(min_value=0.01, max_value=0.6),
    xitip=st.floats(min_value=0.8, max_value=1.0),
    ixspac=st.sampled_from([1, 2]),
)
def test_setx_matches_reference(ref_module, ii, xi0, xitip, ixspac):
    ref = ref_module
    ref["II"], ref["XI0"], ref["XITIP"], ref["IXSPAC"] = ii, xi0, xitip, ixspac
    ref.call("SETX")

    grid = setx(ii, xi0, xitip, ixspac)
    assert grid.dt == ref.get("DT")
    assert_close(grid.t, ref.get("T")[:ii], rtol=0, name="T")
    assert_close(grid.xi, ref.get("XI")[:ii], rtol=0, name="XI")
    assert_close(grid.dxi, ref.get("DXI")[:ii], rtol=0, name="DXI")
    assert_close(grid.xv, ref.get("XV")[: ii + 1], rtol=0, name="XV")
