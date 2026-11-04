"""pyxrotor.fmath against the C library's (the Fortran's) IEEE results."""
import ctypes
import math

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor import fmath as fm

pytestmark = pytest.mark.differential
_libm = ctypes.CDLL("libm.so.6")
_libm.pow.restype = ctypes.c_double
_libm.pow.argtypes = [ctypes.c_double, ctypes.c_double]

SPECIAL = st.sampled_from([0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 3.0, -3.0, 1e-300,
                           5e-324, -5e-324, 1e300, -1e300, math.inf, -math.inf, math.nan])
NUM = st.one_of(SPECIAL, st.floats(allow_nan=True, allow_infinity=True))


def same(a, b):
    return (math.isnan(a) and math.isnan(b)) or (a == b and math.copysign(1, a) == math.copysign(1, b))


@settings(deadline=None, max_examples=2000)
@given(x=NUM, y=NUM)
def test_pow_matches_libm(x, y):
    assert same(fm.pow(x, y), _libm.pow(x, y)), (x, y, fm.pow(x, y), _libm.pow(x, y))


def test_pow_poles():
    assert fm.pow(0.0, -1.0) == math.inf
    assert fm.pow(-0.0, -1.0) == -math.inf
    assert fm.pow(-0.0, -2.0) == math.inf
    assert fm.pow(0.0, -0.5) == math.inf
    assert math.isnan(fm.pow(-2.0, 0.5))
