"""Compare the YP-39 model with the flight runs; writes results.json and prints the table
of README.md.

    PYTHONPATH=src python validation/naca_yp39/compare.py
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("naca_yp39_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)

SETTINGS = {"7.55": {}, "mcrit+0.10": {"mcrit_shift": 0.10}, "mcrit+0.15": {"mcrit_shift": 0.15},
            "mcrit+0.20": {"mcrit_shift": 0.20}, "incompressible": {"mach_scale": 0.01}}
BANDS = [(0.94, 0.99), (0.99, 1.035), (1.035, 1.07)]


def main():
    out = []
    for row in model.runs():
        rec = {k: float(v) for k, v in row.items()}
        for name, kw in SETTINGS.items():
            o = model.run(row, **kw)
            rec[name] = {k: o[k] for k in ("converged", "method", "eta", "beta75", "CT")}
        out.append(rec)
    (HERE / "results.json").write_text(json.dumps(out, indent=0))
    print("mean (model - test working-part efficiency) by tip Mach [converged/runs]")
    for name in SETTINGS:
        line = f"{name:15s}"
        for a, b in BANDS:
            rs = [r for r in out if a <= r["Mt"] < b]
            d = [r[name]["eta"] - r["eta_working"] for r in rs if r[name]["converged"]]
            line += f"  {a:.2f}-{b:.2f}: {np.mean(d):+.3f} [{len(d)}/{len(rs)}]"
        print(line)
    print("test entire propeller minus working part (shank drag):",
          " ".join(f"{r['eta_entire'] - r['eta_working']:+.3f}" for r in out))


if __name__ == "__main__":
    main()
