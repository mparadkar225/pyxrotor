"""Run the Report 747 models at the digitized test points and write results.json.

    python validation/naca747/compare.py [--step 0.1]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "src")]
import model  # noqa: E402


def job(args):
    kind, beta, J = args
    fn = {"single": model.single, "dual": model.dual, "disk4": model.disk4}[kind]
    try:
        r = fn(beta, J)
    except Exception as e:                       # noqa: BLE001 (record and carry on)
        r = {"converged": False, "error": repr(e)}
    return {"kind": kind, "beta": beta, "J": J, **r}


def points(step):
    data = json.loads((HERE / "tests.json").read_text())
    want = {}
    for c in data["curves"]:
        lo, hi = c["J"][0], c["J"][-1]
        js = np.round(np.arange(np.ceil(lo / step) * step, hi + 1e-9, step), 3)
        for J in js:
            want.setdefault((c["rotation"], c["beta_front"]), set()).add(float(J))
    jobs = []
    for (rot, beta), js in sorted(want.items()):
        for J in sorted(js):
            jobs.append((rot, beta, J))
            if rot == "single":
                jobs.append(("disk4", beta, J))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    jobs = points(a.step)
    with ProcessPoolExecutor(a.workers) as ex:
        out = list(ex.map(job, jobs))
    (HERE / "results.json").write_text(json.dumps(out, indent=1))
    bad = [o for o in out if not o["converged"]]
    print(f"{len(out)} points, {len(bad)} not converged")


if __name__ == "__main__":
    main()
