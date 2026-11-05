"""XROTOR against the YP-39 high-speed flight test of a Hamilton Standard 3155-6 three-blade
propeller (validation/naca_yp39; Vogeley 1942, NTRS 20150011182), helical tip Mach number
0.95-1.06. What these tests pin, comparing with the test's "working part" efficiency (the
model has no round shanks):
  * 7.55 with Clark Y Mcrit = 0.85 - 1.2 t loses 0.1-0.9 too much, more as the tip Mach
    number rises;
  * Mcrit + 0.15, the P-47 correction for Clark Y-like sections, follows the test within
    0.04 in every tip Mach band, and converges on every run;
  * the test does lose to compressibility: without it the model is 0.03-0.05 high.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
HERE = Path(__file__).resolve().parents[2] / "validation" / "naca_yp39"
_spec = importlib.util.spec_from_file_location("naca_yp39_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

BANDS = [(0.94, 0.99), (0.99, 1.035), (1.035, 1.07)]


@pytest.fixture(scope="module")
def runs():
    out = {}
    for name, kw in {"7.55": {}, "shift15": {"mcrit_shift": 0.15},
                     "incompressible": {"mach_scale": 0.01}}.items():
        res = []
        for row in model.runs():
            o = model.run(row, **kw)
            res.append((float(row["Mt"]), o["converged"], o["eta"] - float(row["eta_working"])))
        out[name] = res
    return out


def _band_means(res):
    means = []
    for a, b in BANDS:
        d = [e for mt, ok, e in res if a <= mt < b and ok]
        means.append(np.mean(d) if d else np.nan)
    return means


def test_755_loses_far_too_much(runs):
    m = _band_means(runs["7.55"])
    assert m[0] < -0.08 and m[1] < -0.2 and m[2] < -0.5, m
    assert m[0] > m[1] > m[2]


def test_mcrit_plus_015_follows_the_test(runs):
    assert all(ok for _, ok, _ in runs["shift15"])
    m = _band_means(runs["shift15"])
    assert all(abs(x) < 0.04 for x in m), m


def test_compressibility_loss_is_real(runs):
    m = _band_means(runs["incompressible"])
    assert all(0.025 < x < 0.06 for x in m), m
