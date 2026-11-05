"""Compare the static models with RM L51L28's faired curves (figures 6-8), for several
critical-Mach settings, in two ways:

  * at the test's blade angle (results.json): CT and CP;
  * at the test's thrust (results_matched.json): the model's blade angle is swept and its
    CP interpolated at the test's CT. This takes out the blade twist under load, which the
    report says is part of the rise of CT and CP with tip Mach number and which the model
    doesn't have. The ratio model/test CP at matched thrust, normalized by its value at
    Mt = 0.8, is the drift printed: 1 means the model's power rise with tip Mach number
    matches the test's.

    PYTHONPATH=src python validation/naca_l51l28/compare.py
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("naca_l51l28_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

SETTINGS = {"7.55": {}, "mcrit+0.10": {"mcrit_shift": 0.10}, "mcrit+0.15": {"mcrit_shift": 0.15},
            "mcrit+0.20": {"mcrit_shift": 0.20}, "incompressible": {"mach_scale": 0.01}}
MATCHED_MT = (0.4, 0.6, 0.8, 0.9, 0.95, 1.0)


def at_test_angle():
    out = []
    for row in model.tests():
        cld, beta, mt = float(row["design_cl"]), float(row["beta75"]), float(row["Mt"])
        rec = {"design_cl": cld, "beta75": beta, "Mt": mt, "CT_test": float(row["CT"]),
               "CP_test": float(row["CP"])}
        for name, kw in SETTINGS.items():
            m = model.static(cld, beta, mt, **kw)
            rec[name] = {k: m[k] for k in ("converged", "CT", "CP")}
        out.append(rec)
    return out


def polar(cld, mt, kw):
    """(beta, CT, CP) of converged solutions, beta -2 to 20 deg, up to the CT maximum."""
    pts = []
    for b in np.arange(-2.0, 20.01, 1.0):
        o = model.static(cld, float(b), mt, **kw)
        if o["converged"]:
            pts.append((float(b), o["CT"], o["CP"]))
    p = np.array(pts)
    k = 1
    while k < len(p) and p[k, 1] > p[k - 1, 1]:
        k += 1
    return p[:k]


def at_test_thrust():
    tests = model.tests()
    out = []
    for cld in (0.0, 0.3, 0.5):
        for mt in MATCHED_MT:
            pol = {name: polar(cld, mt, kw) for name, kw in SETTINGS.items()}
            for t in tests:
                if float(t["design_cl"]) != cld or abs(float(t["Mt"]) - mt) > 1e-9:
                    continue
                ct, cp = float(t["CT"]), float(t["CP"])
                rec = {"design_cl": cld, "beta75": float(t["beta75"]), "Mt": mt, "CT_test": ct,
                       "CP_test": cp}
                for name, p in pol.items():
                    if len(p) >= 2 and p[0, 1] <= ct <= p[-1, 1]:
                        rec[name] = {"CP": float(np.interp(ct, p[:, 1], p[:, 2])),
                                     "beta75": float(np.interp(ct, p[:, 1], p[:, 0]))}
                    else:
                        rec[name] = None           # no converged solution at that thrust
                out.append(rec)
    return out


def drift(matched, cld, beta, name):
    """{Mt: (model CP / test CP at the test's CT) / (the same at Mt 0.8)}."""
    rs = {r["Mt"]: r for r in matched if r["design_cl"] == cld and r["beta75"] == beta}
    if 0.8 not in rs or rs[0.8][name] is None:
        return {}
    ref = rs[0.8][name]["CP"] / rs[0.8]["CP_test"]
    return {mt: r[name]["CP"] / r["CP_test"] / ref for mt, r in rs.items() if r[name] is not None}


def main():
    fixed = at_test_angle()
    (HERE / "results.json").write_text(json.dumps(fixed, indent=0))
    matched = at_test_thrust()
    (HERE / "results_matched.json").write_text(json.dumps(matched, indent=0))
    print("model CP / test CP at the test's CT (Mt 0.4, 0.6, 0.8), and its drift from Mt 0.8 "
          "(Mt 0.9, 0.95, 1.0); - = no converged solution reaches the test's CT")
    for cld in (0.0, 0.3, 0.5):
        for beta in (0.0, 4.0, 8.0, 12.0):
            if not drift(matched, cld, beta, "incompressible"):
                continue
            line = f"cld {cld} beta {beta:4.1f}"
            rs = {r["Mt"]: r for r in matched if r["design_cl"] == cld and r["beta75"] == beta}
            for name in SETTINGS:
                d = drift(matched, cld, beta, name)
                lo = " ".join(f"{rs[mt][name]['CP'] / rs[mt]['CP_test']:.2f}" if mt in rs and rs[mt][name] else "  - "
                              for mt in (0.4, 0.6, 0.8))
                hi = " ".join(f"{d[mt]:.2f}" if mt in d else "  - " for mt in (0.9, 0.95, 1.0))
                line += f" | {name}: {lo} / {hi}"
            print(line)


if __name__ == "__main__":
    main()
