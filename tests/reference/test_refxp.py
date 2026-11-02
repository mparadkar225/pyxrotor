"""Checks on the xrotor-python reference build.

xrotor-python does not currently reproduce its own bundled test case: Newton's iteration
diverges on the NASA Report 640 propeller at every optimization level tried (-O0, -O1),
with and without the double-precision interface patch, so the fault is upstream rather
than in our build. Until that is fixed it is a structural reference only, not an oracle.
The xfail is strict, so this test starts failing if the problem goes away.
"""
import ctypes

import pytest

from reference import refxp

pytestmark = pytest.mark.reference


def test_interface_patched_to_double(xrotor_python):
    XRotor, _ = xrotor_python
    import xrotor.model as model
    assert model.Performance()._rpm.__class__ is ctypes.c_double


@pytest.mark.xfail(strict=True, reason="xrotor-python's own test case diverges (upstream issue)")
def test_bundled_case_solve_for_thrust(xrotor_python):
    XRotor, Case = xrotor_python
    xr = XRotor()
    xr.print = False
    xr.case = Case.from_dict(refxp.test_case())
    rms = xr.operate(thrust=500)
    perf = xr.performance
    # Expected values from xrotor-python's xrotor/test.py
    assert abs(rms) < 1.0e-7
    assert perf.thrust == pytest.approx(500, abs=0.5)
    assert perf.rpm == pytest.approx(2019, abs=0.5)
    assert perf.efficiency == pytest.approx(0.596, abs=5e-4)
