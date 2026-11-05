"""XROTOR against four four-blade NACA 16-series propellers tested on a YP-47M in the NACA
Cleveland altitude wind tunnel, free-stream Mach 0.3-0.5 (validation/naca_yp47m). The test
efficiencies have no wall correction and one average installation drag coefficient, so one
installation offset per propeller (a constant DCT/J^2, fitted at tip Mach < 0.8 before
breakdown) is removed first. What these tests pin:
  * the offsets behave like a drag area: DCT/J^2 between 0.001 and 0.004;
  * the traced Aeroproducts and Curtiss 732 blades reproduce their reports' activity factors;
  * at helical tip Mach 0.99-1.03 (Curtiss 836), 7.55 with Mcrit = 0.85 - 1.2 t loses 0.13-0.18
    too much; Mcrit + 0.15 is within 0.03, Mcrit + 0.20 and no compressibility are too high
    at tip Mach 1.03;
  * at tip Mach 0.84-0.87 the settings differ by less than 0.06 and Mcrit + 0.15 is within
    0.04 on every curve;
  * past the breakdown at high CP and low J (blade stall), the model with its stall CL of
    design CL + 0.85 misses the drop: it is 0.15-0.25 high on the Aeroproducts curves.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
HERE = Path(__file__).resolve().parents[2] / "validation" / "naca_yp47m"
_spec = importlib.util.spec_from_file_location("naca_yp47m_compare", HERE / "compare.py")
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)
model = compare.model

SETTINGS = {"7.55": {}, "mcrit+0.15": {"mcrit_shift": 0.15}, "mcrit+0.20": {"mcrit_shift": 0.20},
            "incompressible": {"mach_scale": 0.01}}


@pytest.fixture(scope="module")
def results():
    rows = []
    for p in model.points():
        row = {"propeller": p["propeller"], "M0": float(p["M0"]), "J": float(p["J"]),
               "CP": float(p["CP"]), "eta_test": float(p["eta"])}
        row["Mt"] = row["M0"] * np.sqrt(1.0 + (np.pi / row["J"]) ** 2)
        for name, kw in SETTINGS.items():
            if name != "incompressible" and row["Mt"] < 0.80:
                row[name] = {"usable": False, "eta": np.nan}
                continue
            o = model.run(row["propeller"], row["M0"], row["J"], row["CP"], **kw)
            row[name] = {"usable": o["usable"], "eta": o["eta"]}
        rows.append(row)
    cv = compare.curves(rows)
    off = compare.offsets(rows)
    for r in rows:
        r["offset"] = off[r["propeller"]] * r["J"] ** 3 / r["CP"]
    return cv, off


def _mean(rows, name, past=False):
    v = [r[name]["eta"] - r["offset"] - r["eta_test"] for r in rows
         if r["past_breakdown"] == past and r[name]["usable"]]
    return np.mean(v) if v else np.nan


def test_installation_offsets(results):
    _, off = results
    assert set(off) == set(model.GEO)
    assert all(0.001 < f < 0.004 for f in off.values()), off


@pytest.mark.parametrize("prop", ["aeroproducts_H20C-162", "curtiss_732-1C2-0"])
def test_traced_activity_factor(prop):
    g = model.GEO[prop]
    assert abs(g["activity_factor_traced"] / g["activity_factor"] - 1.0) < 0.01


@pytest.mark.parametrize("M0, J", [(0.3, 1.0), (0.5, 1.75)])
def test_sonic_tip_wants_mcrit_plus_015(results, M0, J):
    cv, _ = results
    rows = cv[("curtiss_836-14C2-18R1", M0, J)]
    assert _mean(rows, "7.55") < -0.1
    assert abs(_mean(rows, "mcrit+0.15")) < 0.03
    if J == 1.75:                               # tip Mach 1.03
        assert _mean(rows, "mcrit+0.20") > 0.04
        assert _mean(rows, "incompressible") > 0.08


def test_tip_mach_085_settings_close(results):
    cv, _ = results
    n = 0
    for key, rows in cv.items():
        if 0.83 < rows[0]["Mt"] < 0.88:
            a, b = _mean(rows, "7.55"), _mean(rows, "mcrit+0.15")
            assert abs(b) < 0.04, key
            assert -0.06 < a - b < 0.0, key
            n += 1
    assert n == 5


def test_stall_breakdown_missed(results):
    cv, _ = results
    for J in (1.75, 2.15):
        rows = cv[("aeroproducts_H20C-162", 0.4, J)]
        assert _mean(rows, "incompressible", past=True) > 0.15, J
