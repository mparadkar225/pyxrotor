"""XROTOR's compressibility model against Langley flight tests of three propellers on the
Republic P-47 (validation/naca_p47): the Curtiss No. 714-1C2-12 four-blade (NACA ACR L4L07 /
WR L-177, Vogeley 1944), the Hamilton Standard No. 6507A-2 three-blade (WR L-594, Gardner
1945) and the Curtiss No. 838-1C2-18 three-blade with trailing-edge extensions (RM L7D10,
Gardner 1947). The first tests below are the Curtiss 714's; test_all_three_* cover all.

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


# Printed values that disagree with their own row (kept as printed; the model uses M).
MISPRINTS = {("l7d10_table3.csv", "22-15"): "Mt 1.090, but M 0.618 and J 2.200 give 1.077"}


@pytest.mark.parametrize("name", [t for prop in model.PROPS.values() for t in prop[1]])
def test_tables_are_consistent(name):
    """eta = CT J / CP and Mt = M sqrt(1 + (pi/J)^2), to the tables' rounding (J has 2
    decimals in L-177, 3 in the others): catches transcription errors."""
    for row in model.table(name):
        J, cp, ct, eta, M, mt = (float(row[k]) for k in ("J", "CP", "CT", "eta", "M", "Mt"))
        assert abs(ct * J / cp - eta) < 0.012, row["run"]
        if (name, row["run"]) not in MISPRINTS:
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


ALL_BANDS = [(0.80, 0.90), (0.90, 0.96), (0.96, 1.02), (1.02, 1.08), (1.08, 1.2)]


@pytest.fixture(scope="module")
def all_results():
    out = {}
    for prop in model.PROPS:
        for name, kw in {"7.55": {}, "shift": {"mcrit_shift": 0.15}}.items():
            res = []
            for row in model.runs(prop):
                m = model.run(row, **kw)
                m["eta_test"], m["Mt_test"] = float(row["eta"]), float(row["Mt"])
                res.append(m)
            out[prop, name] = res
    return out


def _mean(res, a, b):
    d = [r["eta"] - r["eta_test"] for r in res if a <= r["Mt_test"] < b and r["converged"]]
    return float(np.mean(d)) if d else float("nan")


@pytest.mark.parametrize("prop", list(model.PROPS))
def test_all_three_755_overpredicts_losses(all_results, prop):
    """For every propeller, 7.55's model with Mcrit = 0.85 - 1.2 t is 0.1 or more too low in
    efficiency by tip Mach 0.96-1.02, and far worse above."""
    res = all_results[prop, "7.55"]
    assert _mean(res, 0.96, 1.02) < -0.10
    assert _mean(res, 1.02, 1.2) < -0.25


@pytest.mark.parametrize("prop", list(model.PROPS))
def test_all_three_mcrit_shift_reproduces_the_losses(all_results, prop):
    """For every propeller, with Mcrit + 0.15 the error at each tip Mach is the propeller's
    own low-speed bias (its mean error below tip Mach 0.90) within +-0.04, up to tip Mach
    1.17, and at least 85% of the runs converge (the rest stall in a limit cycle near the
    local-Mach cap)."""
    res = all_results[prop, "shift"]
    bias = _mean(res, 0.0, 0.90)
    assert 0.0 < bias < 0.07, bias
    for a, b in ALL_BANDS:
        m = _mean(res, a, b)
        if not np.isnan(m):
            assert abs(m - bias) < 0.04, (a, b, m, bias)
    assert sum(r["converged"] for r in res) >= 0.85 * len(res)
