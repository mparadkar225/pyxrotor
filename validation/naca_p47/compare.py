"""Compare the model with the L-177 flight data by tip-Mach band, for several critical-Mach
settings; writes results.json (every run, every setting) and prints the band summary.

    PYTHONPATH=src python validation/naca_p47/compare.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from model import run, runs  # noqa: E402

CASES = {"7.55 model": {}, "Mcrit+0.05": {"mcrit_shift": 0.05}, "Mcrit+0.10": {"mcrit_shift": 0.10},
         "Mcrit+0.15": {"mcrit_shift": 0.15}, "Mcrit+0.20": {"mcrit_shift": 0.20},
         "no compressibility": {"mach_scale": 0.01}}
BANDS = [(0.0, 0.80), (0.80, 0.85), (0.85, 0.90), (0.90, 0.96), (0.96, 1.02), (1.02, 1.2)]


def main():
    out = {}
    for name, kw in CASES.items():
        res = []
        for row in runs():
            m = run(row, **kw)
            m.update(eta_test=float(row["eta"]), CT_test=float(row["CT"]),
                     Mt_test=float(row["Mt"]), table=row["table"],
                     beta_test=float(row["beta_deg"]) if row["beta_deg"] else None)
            res.append(m)
        out[name] = res
    (HERE / "results.json").write_text(json.dumps(out, indent=1))
    print("mean (model - test) efficiency by tip-Mach band [converged/runs]")
    print(" " * 20 + "".join(f"{a:.2f}-{b:.2f}      " for a, b in BANDS))
    for name, res in out.items():
        line = f"{name:20s}"
        for a, b in BANDS:
            rs = [r for r in res if a <= r["Mt_test"] < b]
            d = [r["eta"] - r["eta_test"] for r in rs if r["converged"]]
            line += (f"{np.mean(d):+.3f} [{len(d)}/{len(rs)}]  " if d else f"  --   [0/{len(rs)}]  ")
        print(line)


if __name__ == "__main__":
    main()
