"""Compare the RM L6L27a model with figures 27 and 28; writes results.json and prints the
tables of README.md.

    PYTHONPATH=src python validation/naca_l6l27a/compare.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("naca_l6l27a_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

N_RPS = model.RPM / 60.0
SETTINGS = {"calibrated": {}, "7.55": {"mcrit_shift": 0.0},
            "calibrated, no shank model": {"shank": None}}


def j_of_mph(mph):
    return mph * 0.44704 / (N_RPS * model.D)


def fig28():
    out = []
    for r in model.table("fig28_difference.csv"):
        J = j_of_mph(float(r["airspeed_mph"]))
        rec = {"airspeed_mph": float(r["airspeed_mph"]), "J": J,
               "test_percent": float(r["difference_percent"])}
        for name, kw in SETTINGS.items():
            a = model.run(model.P03, J, 0.2, **{k: v for k, v in kw.items() if k != "shank"})
            b = model.run(model.P03R, J, 0.2, **kw)
            rec[name] = {"eta_03": a["eta"], "eta_03R": b["eta"],
                         "percent": 100.0 * (a["eta"] - b["eta"]) / a["eta"],
                         "shank_power_fraction": b["shank_power_fraction"],
                         "usable": a["usable"] and b["usable"]}
        rec["Mt"] = a["Mt"]
        out.append(rec)
    return out


def fig27():
    out = []
    for r in model.table("fig27_efficiency.csv"):
        J, cp = float(r["J"]), float(r["CP"])
        o = model.run(r["propeller"], J, cp)
        out.append({"propeller": r["propeller"], "CP": cp, "J": J, "eta_test": float(r["eta"]),
                    "eta": o["eta"], "usable": o["usable"]})
    return out


def main():
    a, b = fig28(), fig27()
    (HERE / "results.json").write_text(json.dumps({"fig28": a, "fig27": b}, indent=1))
    print("Figure 28: (eta_03 - eta_03R) / eta_03, percent, CP 0.2, 1140 rpm")
    print("| mph | J | Mt | test |", " | ".join(SETTINGS), "|")
    for r in a:
        print(f"| {r['airspeed_mph']:.0f} | {r['J']:.2f} | {r['Mt']:.2f} | {r['test_percent']:.1f} |",
              " | ".join(f"{r[k]['percent']:.1f}" for k in SETTINGS), "|")
    print("\nFigure 27: model - test efficiency (calibrated Mcrit, shank model)")
    for prop in (model.P03, model.P03R):
        d = [x["eta"] - x["eta_test"] for x in b if x["propeller"] == prop and x["usable"]]
        print(prop, f"mean {np.mean(d):+.3f}, range {min(d):+.3f} to {max(d):+.3f}, n {len(d)}")


if __name__ == "__main__":
    main()
