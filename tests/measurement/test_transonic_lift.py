"""The transonic lift loss (pyxrotor.transonic) on the P-47 flight tests
(validation/transonic_lift): with the calibrated drag (Mcrit + 0.15), the 2-D lift
divergence (shift 0) loses 0.09-0.12 too much at tip Mach 1.08-1.2. The same +0.15 shift
on the lift divergence keeps every band within 0.08 of the low-speed bias, with at most 3
unconverged runs per propeller. On the plain Mcrit the lift loss moves the error by under
0.1 (the drag rise dominates)."""
import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.measurement
spec = importlib.util.spec_from_file_location(
    "transonic_lift_compare",
    Path(__file__).resolve().parents[2] / "validation" / "transonic_lift" / "compare.py")
cmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmp)


@pytest.mark.parametrize("prop", ["l177", "l594", "l7d10"])
def test_p47_lift_loss(prop):
    shifted = cmp.p47(prop, "calibrated, lift loss +0.15")
    assert all(abs(v) < 0.08 for v in shifted["band_minus_bias"].values() if v is not None)
    assert shifted["unconverged"] <= 3
    raw = cmp.p47(prop, "calibrated, lift loss +0")["band_minus_bias"]["1.08-1.2"]
    assert raw < -0.085
    plain = cmp.p47(prop, "plain")["band_minus_bias"]
    plain_ll = cmp.p47(prop, "plain, lift loss +0")["band_minus_bias"]
    for k, v in plain.items():
        if v is not None and plain_ll[k] is not None:
            assert abs(plain_ll[k] - v) < 0.1, k
