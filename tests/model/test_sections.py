"""Section estimates (pyxrotor.sections): the plain thickness estimate reproduces the
formulas the measurement comparisons were built on, and the calibrated critical Mach number
follows its stated rule. Whether the calibration fits the data is tested in
tests/measurement/test_calibrated_sections.py."""
import numpy as np
import pytest

from pyxrotor import sections
from pyxrotor.constants import PI

pytestmark = pytest.mark.model


@pytest.mark.parametrize("t", [0.02, 0.044, 0.06, 0.072, 0.15, 0.25, 0.4])
def test_thickness_estimate(t):
    """mcrit_shift=0: Mcrit = max(0.85 - 1.2 t, 0.5) exactly, for both families; a number
    is added to it exactly (as the comparisons added their shifts)."""
    base = max(0.85 - 1.2 * t, 0.5)
    assert sections.clark_y(t, 0.5, 0.0).mcrit == base
    assert sections.naca16(t, 0.5, 0.4, 0.0).mcrit == base
    assert sections.naca16(t, 0.5, 0.4, 0.15).mcrit == base + 0.15


def test_clark_y_values():
    """The Clark Y estimate as validation/naca747 defined it (R&M 2357's sections)."""
    s = sections.clark_y(0.10, 0.6, 0.0)
    assert s.a0 == -3.7 * PI / 180.0 and s.dclda == 6.3 - 4.0 * (0.10 - 0.066)
    assert (s.clmax, s.clmin, s.cldmin, s.cmcon) == (1.3, -0.5, 0.35, -0.1)
    assert s.cdmin == 0.0065 + 0.25 * 0.10 * 0.10 and s.dcdcl2 == 0.007 + 0.01 * (0.10 - 0.066)
    assert sections.clark_y(0.25, 0.3, 0.0).clmax == 1.1
    assert sections.clark_y(0.25, 0.3, 0.0).a0 == -6.0 * PI / 180.0


def test_naca16_values():
    s = sections.naca16(0.05, 0.8, 0.5, 0.0)
    assert s.a0 == -0.5 / (2.0 * PI) and s.dclda == 6.0
    assert (s.clmax, s.clmin, s.cldmin, s.cmcon) == (0.5 + 0.85, 0.5 - 0.85, 0.5, -0.125)
    assert sections.naca16(0.15, 0.3, 0.5, 0.0).clmax == 0.5 + 1.1


def test_calibrated_shift_rule():
    """+0.15 at 5% and thicker, +0.20 at 4.5% and thinner, linear between; it is the
    default (mcrit_shift=None)."""
    f = sections.calibrated_shift
    assert f(0.05) == f(0.06) == f(0.3) == 0.15
    assert f(0.045) == f(0.02) == f(0.0) == 0.20
    assert abs(f(0.0475) - 0.175) < 1e-12
    for t in (0.03, 0.047, 0.065):
        assert sections.mcrit(t) == max(0.85 - 1.2 * t, 0.5) + f(t)
        assert sections.naca16(t, 0.9, 0.3).mcrit == sections.mcrit(t)
        assert sections.clark_y(t, 0.9).mcrit == sections.mcrit(t)


def test_calibrated_mcrit_falls_with_thickness():
    """No jump up with thickness: a thicker section never gets a higher Mcrit."""
    t = np.linspace(0.0, 0.4, 4001)
    m = np.array([sections.mcrit(x) for x in t])
    assert np.all(np.diff(m) <= 1e-15)
