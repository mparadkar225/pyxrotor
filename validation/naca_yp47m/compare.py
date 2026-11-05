"""Compare the YP-47M propeller model with the altitude-wind-tunnel efficiency curves;
writes results.json and prints the tables of README.md.

    PYTHONPATH=src python validation/naca_yp47m/compare.py

The test efficiencies are not corrected for tunnel-wall constriction and use one average
installation drag coefficient; the reports call them trends only. Away from stall and at
low tip Mach number the model's efficiency is above the test's by an amount that behaves
like a constant drag area: (eta_model - eta_test) CP / J^3 = DCT / J^2 is nearly constant
for each propeller. So one installation offset f per propeller is fitted on those points,
and eta_model - f J^3 / CP is compared everywhere else.
"""
from __future__ import annotations

import importlib.util
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location(
    "naca_yp47m_model", Path(__file__).resolve().parent / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

SETTINGS = {"7.55": dict(mcrit_shift=0.0), "mcrit+0.15": dict(mcrit_shift=0.15),
            "mcrit+0.20": dict(mcrit_shift=0.20), "incompressible": dict(mach_scale=0.01),
            "mcrit+0.15, stall cld+0.5": dict(mcrit_shift=0.15, stall_margin=0.5)}
FIT_MT = 0.80          # offset fitted on curves below this helical tip Mach number
BREAKDOWN = 0.08       # a point is past breakdown once the test is this far below its curve's peak


def curves(results):
    out = defaultdict(list)
    for r in results:
        out[(r["propeller"], r["M0"], r["J"])].append(r)
    for rows in out.values():
        rows.sort(key=lambda r: r["CP"])
        peak = 0.0
        past = False
        for r in rows:
            peak = max(peak, r["eta_test"])
            past = past or r["eta_test"] < peak - BREAKDOWN
            r["past_breakdown"] = past
    return out


def offsets(results):
    """Per propeller: median of (eta_incompressible - eta_test) CP / J^3 over the fit points."""
    f = defaultdict(list)
    for r in results:
        o = r["incompressible"]
        if r["Mt"] < FIT_MT and not r["past_breakdown"] and o["usable"]:
            f[r["propeller"]].append((o["eta"] - r["eta_test"]) * r["CP"] / r["J"] ** 3)
    return {p: float(np.median(v)) for p, v in f.items()}


def residuals(results, off):
    for r in results:
        corr = off[r["propeller"]] * r["J"] ** 3 / r["CP"]
        r["installation_offset"] = corr
        for name in SETTINGS:
            o = r[name]
            o["residual"] = o["eta"] - corr - r["eta_test"] if o["usable"] else None


def main():
    results = []
    for p in model.points():
        row = {"propeller": p["propeller"], "M0": float(p["M0"]), "J": float(p["J"]),
               "CP": float(p["CP"]), "eta_test": float(p["eta"]), "figure": p["figure"]}
        for name, kw in SETTINGS.items():
            o = model.run(row["propeller"], row["M0"], row["J"], row["CP"], **kw)
            row[name] = {k: o[k] for k in ("eta", "beta75", "CT", "method", "usable")}
            if name == "7.55":
                row["Mt"] = o["Mt"]          # helical tip Mach number at the test's M0 and J
        results.append(row)
    cv = curves(results)
    off = offsets(results)
    residuals(results, off)
    Path(model.HERE / "results.json").write_text(json.dumps(
        {"offsets_dCT_over_J2": off, "points": results}, indent=1))
    report(cv, off)


def fmt(v):
    return "—" if v is None or not np.isfinite(v) else f"{v:+.3f}"


def report(cv, off):
    print("installation offset f = DCT/J^2:", {k: round(v, 5) for k, v in off.items()})
    names = list(SETTINGS)
    print("\nMean residual (model - offset - test) per curve, before / past breakdown [n]")
    print("| propeller | M0 | J | Mt |", " | ".join(names), "|")
    for key, rows in sorted(cv.items(), key=lambda kv: -kv[1][0]["Mt"]):
        cells = []
        for name in names:
            parts = []
            for past in (False, True):
                v = [r[name]["residual"] for r in rows if r["past_breakdown"] == past]
                if not v:
                    parts.append("")
                    continue
                ok = [x for x in v if x is not None]
                parts.append(f"{np.mean(ok):+.3f} [{len(ok)}/{len(v)}]" if ok else f"— [0/{len(v)}]")
            cells.append(" / ".join(p for p in parts if p))
        print(f"| {key[0]} | {key[1]} | {key[2]} | {rows[0]['Mt']:.2f} |", " | ".join(cells), "|")


if __name__ == "__main__":
    main()
