"""The calibrated critical Mach number of pyxrotor.sections (mcrit_shift=None: + 0.15 at
5% thickness and above, + 0.20 at 4.5% and below) against every transonic dataset in
validation/, with pyxrotor.robust's power solve where the test fixed the power:
  * P-47 flight tests (Curtiss 714, HS 6507A-2, Curtiss 838; tips 5-7%): at tip Mach 0.9-1.17
    every band within 0.04 of the propeller's own low-speed bias; at most 1 run in 112
    without a usable solution;
  * YP-39 flight test (Clark Y, 3155-6): within 0.035 of the working-part efficiency in every
    tip-Mach band, tip Mach 0.95-1.06;
  * YP-47M tunnel (Curtiss 836, tip Mach 0.99-1.03): within 0.035 after the installation
    offset, before the breakdown; at J 1.75 without the lightest-loaded point (CP 0.19,
    where the test is 0.17 lower and matches 7.55's Mcrit instead; validation/naca_yp47m);
  * RM L51L28 static tests (16-series, 4.4% at 0.94R): from Mt 0.8 to 1.0 the power at the
    test's thrust rises within 10% of the test's on the cambered blades and on the
    symmetric blade at 8 deg; at 12 deg the symmetric blade's rises 36% too much (+ 0.20
    alone gives 26%);
  * NASA SR-2 propfan (2% tips): within 0.025 of its Mach 0.6 bias up to Mach 0.85.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
VAL = Path(__file__).resolve().parents[2] / "validation"


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, VAL / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.mark.parametrize("prop", ["l177", "l594", "l7d10"])
def test_p47(prop):
    p47 = _load("cal_p47_model", "naca_p47/model.py")
    res = [(o["Mt"], o["converged"], o["eta"] - float(row["eta"]))
           for row in p47.runs(prop) for o in [p47.run(row, mcrit_shift=None)]]
    bias = np.mean([e for mt, ok, e in res if mt < 0.9 and ok])
    assert 0.0 < bias < 0.06, bias
    for a, b in [(0.9, 0.96), (0.96, 1.02), (1.02, 1.08), (1.08, 1.2)]:
        e = [e for mt, ok, e in res if a <= mt < b and ok]
        assert abs(np.mean(e) - bias) < 0.04, (a, b, np.mean(e), bias)
    assert sum(1 for _, ok, _ in res if not ok) <= 1


def test_yp39():
    yp39 = _load("cal_yp39_model", "naca_yp39/model.py")
    d = [(o["Mt"], o["converged"], o["eta"] - float(row["eta_working"]))
         for row in yp39.runs() for o in [yp39.run(row, mcrit_shift=None)]]
    assert all(ok for _, ok, _ in d)
    for a, b in [(0.94, 0.99), (0.99, 1.035), (1.035, 1.07)]:
        assert abs(np.mean([e for mt, _, e in d if a <= mt < b])) < 0.035, (a, b)


def test_yp47m_sonic_tips():
    cmp = _load("cal_yp47m_compare", "naca_yp47m/compare.py")
    model = cmp.model
    rows = []
    for p in model.points():
        row = {"propeller": p["propeller"], "M0": float(p["M0"]), "J": float(p["J"]),
               "CP": float(p["CP"]), "eta_test": float(p["eta"])}
        row["Mt"] = row["M0"] * np.sqrt(1.0 + (np.pi / row["J"]) ** 2)
        o = model.run(row["propeller"], row["M0"], row["J"], row["CP"], mach_scale=0.01)
        row["incompressible"] = {"usable": o["usable"], "eta": o["eta"]}
        if row["propeller"].startswith("curtiss_836") and row["Mt"] > 0.98:
            o = model.run(row["propeller"], row["M0"], row["J"], row["CP"], mcrit_shift=None)
            row["cal"] = {"usable": o["usable"], "eta": o["eta"]}
        rows.append(row)
    cv = cmp.curves(rows)
    off = cmp.offsets(rows)
    for key in [("curtiss_836-14C2-18R1", 0.3, 1.0), ("curtiss_836-14C2-18R1", 0.5, 1.75)]:
        e = [r["cal"]["eta"] - off[key[0]] * r["J"] ** 3 / r["CP"] - r["eta_test"]
             for r in cv[key] if not r["past_breakdown"] and r["cal"]["usable"]
             and not (r["J"] == 1.75 and r["CP"] < 0.25)]
        assert len(e) >= 2 and abs(np.mean(e)) < 0.035, (key, e)


def test_l51l28():
    cmp = _load("cal_l51l28_compare", "naca_l51l28/compare.py")
    cmp.SETTINGS = {"cal": {"mcrit_shift": None}}
    cmp.MATCHED_MT = (0.8, 1.0)
    m = cmp.at_test_thrust()
    for cld in (0.0, 0.3, 0.5):
        for beta in (4.0, 8.0, 12.0):
            d = cmp.drift(m, cld, beta, "cal")
            if (cld, beta) == (0.0, 4.0):
                continue                      # no solution reaches the test's thrust at Mt 1.0
            if (cld, beta) == (0.0, 12.0):
                assert 1.2 < d[1.0] < 1.45, d
            else:
                assert abs(d[1.0] - 1.0) < 0.1, (cld, beta, d)


def test_sr2():
    pf = _load("cal_propfan_model", "nasa_propfan/model.py")
    rows = [r for r in pf.table("design_point_eta.csv")
            if r["prop"] == "SR-2" and r["source"] == "TM-79124 fig. 14"]
    err = {}
    for r in rows:
        o = pf.run("SR-2", float(r["M"]), 3.06, cp=1.7, mcrit_shift=None)
        assert o["usable"], r["M"]
        err[float(r["M"])] = o["eta"] - float(r["eta"]) / 100.0
    bias = err[0.6]
    assert all(abs(e - bias) < 0.025 for e in err.values()), err
