"""Rerun the transonic propeller comparisons with sections fitted from NASA CR-176883's data
banks (pyxrotor.databank.clark_y_section / naca16_section) in place of pyxrotor.sections'
thickness estimates, and write results.json.

Every dataset's own model is used unchanged except for its section builder. Settings:
  est+0     sections' estimate, plain Mcrit 0.85 - 1.2 t (what 7.55 users would type)
  est+cal   sections' estimate, calibrated Mcrit (+0.15 / +0.20)
  bank+S    the bank's fitted section, its fitted Mcrit + S

    PYTHONPATH=src python validation/cr176883_sections/compare.py [dataset ...]
"""
from __future__ import annotations

import importlib.util
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pyxrotor import databank, sections

HERE = Path(__file__).resolve().parent
VAL = HERE.parent
BANK_SHIFTS = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25)
SETTINGS = ["est+0", "est+cal"] + [f"bank+{s:.2f}" for s in BANK_SHIFTS]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, VAL / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _shift(setting):
    """(builders use the bank?, mcrit_shift to pass)"""
    if setting == "est+0":
        return False, 0.0
    if setting == "est+cal":
        return False, None
    return True, float(setting.split("+")[1])


def _patch(m, bank):
    """Swap a model module's section builder for the bank's."""
    if not bank:
        return
    if hasattr(m, "sections_"):
        m.sections_ = SimpleNamespace(clark_y=databank.clark_y_section,
                                      naca16=databank.naca16_section, mcrit=sections.mcrit)
    if hasattr(m, "naca16"):
        m.naca16 = databank.naca16_section
    if hasattr(m, "section16"):
        m.section16 = databank.naca16_section


def _bands(res, bands, bias=0.0):
    out = {}
    for a, b in bands:
        e = [e for mt, ok, e in res if a <= mt < b and ok]
        out[f"{a}-{b}"] = None if not e else round(float(np.mean(e)) - bias, 4)
    return out


def p47(prop, setting):
    bank, shift = _shift(setting)
    m = _load(f"cr_p47_{prop}_{setting}", "naca_p47/model.py")
    _patch(m, bank)
    res = [(o["Mt"], o["converged"], o["eta"] - float(row["eta"]))
           for row in m.runs(prop) for o in [m.run(row, mcrit_shift=shift)]]
    bias = float(np.mean([e for mt, ok, e in res if mt < 0.9 and ok]))
    return {"bias": round(bias, 4),
            "band_minus_bias": _bands(res, [(0.9, 0.96), (0.96, 1.02), (1.02, 1.08), (1.08, 1.2)], bias),
            "unconverged": sum(1 for _, ok, _ in res if not ok), "runs": len(res)}


def yp39(setting):
    bank, shift = _shift(setting)
    m = _load(f"cr_yp39_{setting}", "naca_yp39/model.py")
    _patch(m, bank)
    res = [(o["Mt"], o["converged"], o["eta"] - float(row["eta_working"]))
           for row in m.runs() for o in [m.run(row, mcrit_shift=shift)]]
    return {"band": _bands(res, [(0.94, 0.99), (0.99, 1.035), (1.035, 1.07)]),
            "unconverged": sum(1 for _, ok, _ in res if not ok)}


def yp47m(setting):
    """Curtiss 836 at tip Mach > 0.98 (the sonic tips), after each propeller's installation
    offset as in test_calibrated_sections; the offset is refitted for the bank sections from
    their own incompressible solves (the offset is a drag area, not compressibility)."""
    bank, shift = _shift(setting)
    cmp = _load(f"cr_yp47m_{setting}", "naca_yp47m/compare.py")
    model = cmp.model
    _patch(model, bank)
    rows = []
    for p in model.points():
        row = {"propeller": p["propeller"], "M0": float(p["M0"]), "J": float(p["J"]),
               "CP": float(p["CP"]), "eta_test": float(p["eta"])}
        row["Mt"] = row["M0"] * np.sqrt(1.0 + (np.pi / row["J"]) ** 2)
        o = model.run(row["propeller"], row["M0"], row["J"], row["CP"], mach_scale=0.01)
        row["incompressible"] = {"usable": o["usable"], "eta": o["eta"]}
        if row["propeller"].startswith("curtiss_836") and row["Mt"] > 0.98:
            o = model.run(row["propeller"], row["M0"], row["J"], row["CP"], mcrit_shift=shift)
            row["cal"] = {"usable": o["usable"], "eta": o["eta"]}
        rows.append(row)
    cv = cmp.curves(rows)
    off = cmp.offsets(rows)
    out = {"offsets": {k: round(float(v), 5) for k, v in off.items()}}
    for key in [("curtiss_836-14C2-18R1", 0.3, 1.0), ("curtiss_836-14C2-18R1", 0.5, 1.75)]:
        e = [r["cal"]["eta"] - off[key[0]] * r["J"] ** 3 / r["CP"] - r["eta_test"]
             for r in cv[key] if not r["past_breakdown"] and r["cal"]["usable"]
             and not (r["J"] == 1.75 and r["CP"] < 0.25)]
        out[f"M0 {key[1]} J {key[2]}"] = None if not e else round(float(np.mean(e)), 4)
    return out


def l51l28(setting):
    bank, shift = _shift(setting)
    cmp = _load(f"cr_l51l28_{setting}", "naca_l51l28/compare.py")
    _patch(cmp.model, bank)
    cmp.SETTINGS = {"s": {"mcrit_shift": shift}}
    cmp.MATCHED_MT = (0.4, 0.8, 1.0)
    mt = cmp.at_test_thrust()
    out = {}
    for cld in (0.0, 0.3, 0.5):
        for beta in (4.0, 8.0, 12.0):
            rs = {r["Mt"]: r for r in mt if r["design_cl"] == cld and r["beta75"] == beta}
            low = rs.get(0.4)
            d = cmp.drift(mt, cld, beta, "s")
            out[f"cld {cld} beta {beta}"] = {
                "CP_ratio_Mt0.4": None if not low or low["s"] is None
                else round(low["s"]["CP"] / low["CP_test"], 3),
                "drift_Mt1.0": None if 1.0 not in d else round(d[1.0], 3)}
    return out


def sr2(setting):
    bank, shift = _shift(setting)
    pf = _load(f"cr_pf_{setting}", "nasa_propfan/model.py")
    _patch(pf, bank)
    rows = [r for r in pf.table("design_point_eta.csv")
            if r["prop"] == "SR-2" and r["source"] == "TM-79124 fig. 14"]
    err = {}
    for r in rows:
        o = pf.run("SR-2", float(r["M"]), 3.06, cp=1.7, mcrit_shift=shift)
        err[float(r["M"])] = None if not o["usable"] else o["eta"] - float(r["eta"]) / 100.0
    bias = err[0.6]
    return {"bias_M0.6": None if bias is None else round(bias, 4),
            "minus_bias": {str(k): None if v is None or bias is None else round(v - bias, 4)
                           for k, v in sorted(err.items())}}


JOBS = {"p47_l177": lambda s: p47("l177", s), "p47_l594": lambda s: p47("l594", s),
        "p47_l7d10": lambda s: p47("l7d10", s), "yp39": yp39, "yp47m": yp47m,
        "l51l28": l51l28, "sr2": sr2}


def _job(args):
    name, setting = args
    try:
        return name, setting, JOBS[name](setting)
    except Exception as e:                       # record, keep the other jobs
        return name, setting, {"error": repr(e)}


def main(names):
    out_path = HERE / "results.json"
    results = json.loads(out_path.read_text()) if out_path.exists() else {}
    # results already in results.json are kept (delete the file to redo everything)
    work = [(n, s) for n in names for s in SETTINGS if s not in results.get(n, {})]
    with ProcessPoolExecutor() as ex:
        for name, setting, r in ex.map(_job, work):
            results.setdefault(name, {})[setting] = r
            print(name, setting, json.dumps(r), flush=True)
            out_path.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:] or list(JOBS))
