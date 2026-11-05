"""XROTOR's compressibility model against flight tests of a Curtiss No. 714-1C2-12 four-blade
propeller on a Republic P-47C (NACA ACR L4L07 / WR L-177, Vogeley 1944; validation/naca_p47).

62 runs: climbs at tip Mach numbers 0.78-0.95 and high-speed runs at airplane Mach up to
0.71, tip Mach up to 1.10. Each is modeled at the test's J, power coefficient, rpm, density
and Mach number (POWE at fixed rpm). What these tests pin:
  * the transcribed tables are self-consistent;
  * without compressibility the model is 0.02-0.05 high in efficiency at every tip Mach
    (the Clark Y section estimate, as on R&M 2357);
  * with 7.55's compressibility model and Mcrit = 0.85 - 1.2 t the losses start too early
    and grow far too fast: 0.09 too low by tip Mach 0.92-0.96, 0.15 at 1.0, and above 1.02
    efficiency 0.3-0.4 too low or no converged solution;
  * raising every section's Mcrit by 0.15 brings the losses into line: model minus test is
    then the same +0.02-0.04 bias at every tip Mach up to 1.10.
"""
import csv
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
HERE = Path(__file__).resolve().parents[2] / "validation" / "naca_p47"
_spec = importlib.util.spec_from_file_location("naca_p47_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

BANDS = [(0.0, 0.85), (0.85, 0.90), (0.90, 0.96), (0.96, 1.02), (1.02, 1.2)]


@pytest.fixture(scope="module")
def results():
    out = {}
    for name, kw in {"7.55": {}, "shift": {"mcrit_shift": 0.15},
                     "incompressible": {"mach_scale": 0.01}}.items():
        res = []
        for row in model.runs():
            m = model.run(row, **kw)
            m["eta_test"], m["Mt_test"] = float(row["eta"]), float(row["Mt"])
            res.append(m)
        out[name] = res
    return out


def band_errors(res):
    """Mean (model - test) efficiency per tip-Mach band (converged runs), and the number of
    unconverged runs per band."""
    means, failed = [], []
    for a, b in BANDS:
        rs = [r for r in res if a <= r["Mt_test"] < b]
        d = [r["eta"] - r["eta_test"] for r in rs if r["converged"]]
        means.append(float(np.mean(d)) if d else float("nan"))
        failed.append(sum(not r["converged"] for r in rs))
    return np.array(means), failed


@pytest.mark.parametrize("name", ["l177_table1.csv", "l177_table2.csv"])
def test_tables_are_consistent(name):
    """eta = CT J / CP and Mt = M sqrt(1 + (pi/J)^2), to the tables' rounding (J has 2
    decimals): catches transcription errors."""
    for row in model.table(name):
        J, cp, ct, eta, M, mt = (float(row[k]) for k in ("J", "CP", "CT", "eta", "M", "Mt"))
        assert abs(ct * J / cp - eta) < 0.012, row["run"]
        assert abs(M * math.sqrt(1.0 + (math.pi / J) ** 2) - mt) < 0.013, row["run"]


def test_low_tip_mach_matches(results):
    """Below tip Mach 0.85 compressibility hardly matters: 7.55 and the incompressible model
    agree with each other and with the test (bias +0.00 to +0.05)."""
    for name in ("7.55", "incompressible"):
        d = [r["eta"] - r["eta_test"] for r in results[name] if r["Mt_test"] < 0.85]
        assert all(r["converged"] for r in results[name] if r["Mt_test"] < 0.85)
        assert -0.04 < np.mean(d) < 0.05, name


def test_incompressible_bias_is_flat(results):
    """Without compressibility the model is high by a roughly constant amount below tip
    Mach 0.96, the section-model bias, and the test's losses then show above that."""
    m, _ = band_errors(results["incompressible"])
    assert np.all((m[:3] > 0.01) & (m[:3] < 0.07)), m
    assert m[3] > 0.07 and m[4] > 0.15, m


def test_755_overpredicts_compressibility_losses(results):
    """XROTOR 7.55's model with Mcrit = 0.85 - 1.2 t: losses too early and far too large."""
    m, failed = band_errors(results["7.55"])
    assert m[2] < -0.06 and m[3] < -0.10, m
    assert failed[4] >= 6 or m[4] < -0.3, (m, failed)


def test_mcrit_shift_reproduces_the_losses(results):
    """With every section's Mcrit raised by 0.15 the model minus test is the incompressible
    bias (+0.01 to +0.05) at every tip Mach up to 1.10, and every run converges."""
    m, failed = band_errors(results["shift"])
    assert not any(failed), failed
    assert np.all((m > 0.0) & (m < 0.055)), m
