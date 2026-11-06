"""Compare the SR-2 and SR-3 models with the NASA Lewis 8x6 ft tunnel data; writes
results.json and prints the summary tables of README.md.

    PYTHONPATH=src python validation/nasa_propfan/compare.py

  A. design point, J = 3.06, CP = 1.7 (1.695), against Mach number: SR-2 from TM-79124
     fig. 14, SR-3 from CR-3505 Table I (and fig. 14 as a check). Net efficiency, blade angle.
  B. TM-79124 figs. 12-13: net efficiency against J at 70-100% of design power loading
     (CP = 0.05933 f J^3) at Mach 0.6 and 0.8.
  C. SR-3 at a fixed blade angle (CR-3505 p. 19): CP falls 9% and 25% from Mach 0.70 to
     0.80 and 0.85 at JREF = 3.1, beta = 60.5 deg, and CT by 11% and 29%. The model is set
     to the blade angle that gives the test's CP at Mach 0.70 (it needs about 5 deg less, as
     did Hamilton Standard's own lifting-line programs) and held there.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("nasa_propfan_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

SETTINGS = {
    "7.55": {}, "mcrit+0.10": {"mcrit_shift": 0.10}, "mcrit+0.15": {"mcrit_shift": 0.15},
    "mcrit+0.20": {"mcrit_shift": 0.20},
    "sweep Mcrit/cos": {"sweep": "mcrit"}, "cosine": {"sweep": "cosine"},
    "cosine+0.10": {"sweep": "cosine", "mcrit_shift": 0.10},
    "cosine+0.15": {"sweep": "cosine", "mcrit_shift": 0.15},
    "core sweep": {"sweep": "core"}, "core sweep+0.05": {"sweep": "core", "mcrit_shift": 0.05},
    "core sweep+0.10": {"sweep": "core", "mcrit_shift": 0.10},
    "core sweep, incompressible": {"sweep": "core", "mach_scale": 0.01},
    "incompressible": {"mach_scale": 0.01},
}
SWEEP_ONLY = {"sweep Mcrit/cos", "cosine", "cosine+0.10", "cosine+0.15", "core sweep",
              "core sweep+0.05", "core sweep+0.10", "core sweep, incompressible"}  # SR-3 only


def settings(prop):
    return {k: v for k, v in SETTINGS.items() if prop == "SR-3" or k not in SWEEP_ONLY}


def slim(o):
    return {k: o[k] for k in ("usable", "converged", "CP", "CT", "eta", "beta75") if k in o}


def design_point():
    out = []
    for row in model.table("design_point_eta.csv"):
        if row["basis"] != "Dref":
            continue
        prop, M, J, cp = row["prop"], float(row["M"]), float(row["J"]), float(row["CP"])
        rec = {"prop": prop, "source": row["source"], "M": M, "J": J, "CP": cp,
               "eta_test": float(row["eta"]) / 100.0}
        for name, kw in settings(prop).items():
            rec[name] = slim(model.run(prop, M, J, cp=cp, **kw))
        out.append(rec)
    return out


def power_loading():
    out = []
    for row in model.table("power_loading_eta.csv"):
        prop, M, J, cp = row["prop"], float(row["M"]), float(row["J"]), float(row["CP"])
        rec = {"prop": prop, "M": M, "J": J, "CP": cp, "percent": int(row["percent_design_loading"]),
               "eta_test": float(row["eta"]) / 100.0}
        for name, kw in settings(prop).items():
            rec[name] = slim(model.run(prop, M, J, cp=cp, **kw))
        out.append(rec)
    return out


def fixed_pitch():
    """CR-3505: SR-3 at JREF 3.1, beta 60.5 deg. For each setting the blade angle giving the
    test's CP 1.64 at Mach 0.70, then CP and CT at 0.80 and 0.85 relative to 0.70."""
    test = {"CP": {0.70: 1.64, 0.80: 1.49, 0.85: 1.23}, "CT_ratio": {0.80: 0.89, 0.85: 0.71}}
    out = {"test": test}
    for name, kw in settings("SR-3").items():
        ref = model.run("SR-3", 0.70, 3.1, cp=1.64, **kw)
        if not ref["usable"]:
            out[name] = None
            continue
        b = ref["beta75"]
        pts = {0.70: slim(model.run("SR-3", 0.70, 3.1, beta75=b, **kw))}
        for M in (0.75, 0.80, 0.85):
            pts[M] = slim(model.run("SR-3", M, 3.1, beta75=b, **kw))
        out[name] = {"beta75": b, "beta_offset": b - 60.5,
                     "points": {str(k): v for k, v in pts.items()}}
    return out


def loss(rows, name, M0=0.6):
    """eta(M) - eta(M0) of the model and of the test, by Mach number."""
    by = {r["M"]: r for r in rows}
    if M0 not in by or not by[M0][name]["usable"]:
        return {}
    res = {}
    for M, r in sorted(by.items()):
        if r[name]["usable"]:
            res[M] = (r[name]["eta"] - by[M0][name]["eta"], r["eta_test"] - by[M0]["eta_test"])
    return res


def main():
    A = design_point()
    B = power_loading()
    C = fixed_pitch()
    (HERE / "results.json").write_text(json.dumps({"design_point": A, "power_loading": B,
                                                   "fixed_pitch": C}, indent=0))
    for prop, src in (("SR-2", "TM-79124 fig. 14"), ("SR-3", "CR-3505 Table I row 1")):
        rows = [r for r in A if r["prop"] == prop and r["source"] == src]
        print(f"\n{prop} at J 3.06, CP 1.7 ({src}): net efficiency, test and model "
              "(model - test), blade angle")
        print("        " + "  ".join(f"M {r['M']:.2f}" for r in rows))
        print("test    " + "  ".join(f"{r['eta_test']:.3f} " for r in rows))
        for name in settings(prop):
            print(f"{name:16s}" + "  ".join(
                (f"{r[name]['eta'] - r['eta_test']:+.3f}" if r[name]["usable"] else "  -   ") +
                (f" b{r[name]['beta75']:.1f}" if r[name]["usable"] else "      ") for r in rows))
    print("\nTM-79124 figs 12-13, mean (model - test) net efficiency [usable/points]")
    for prop in ("SR-2", "SR-3"):
        for M in (0.6, 0.8):
            rows = [r for r in B if r["prop"] == prop and r["M"] == M]
            line = f"{prop} M {M}: "
            for name in settings(prop):
                d = [r[name]["eta"] - r["eta_test"] for r in rows if r[name]["usable"]]
                line += f"{name} {np.mean(d):+.3f} [{len(d)}/{len(rows)}]  " if d else f"{name} - "
            print(line)
    print("\nSR-3 at fixed pitch (JREF 3.1): CP and CT relative to Mach 0.70")
    t = C["test"]
    print(f"test: CP {t['CP'][0.80] / t['CP'][0.70]:.3f} {t['CP'][0.85] / t['CP'][0.70]:.3f}  "
          f"CT {t['CT_ratio'][0.80]:.2f} {t['CT_ratio'][0.85]:.2f}  (beta 60.5)")
    for name, v in C.items():
        if name == "test" or v is None:
            continue
        p = v["points"]
        if not p["0.7"]["usable"]:
            continue
        f = lambda M, q: (p[M][q] / p["0.7"][q]) if p[M]["usable"] else float("nan")
        print(f"{name:16s} beta {v['beta75']:.2f} ({v['beta_offset']:+.2f}): CP {f('0.8', 'CP'):.3f} "
              f"{f('0.85', 'CP'):.3f}  CT {f('0.8', 'CT'):.3f} {f('0.85', 'CT'):.3f}")


if __name__ == "__main__":
    main()
