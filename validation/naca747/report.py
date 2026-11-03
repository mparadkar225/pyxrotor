"""Compare results.json (the models) with tests.json (the digitized tunnel data).

Only J where both CT and CP were digitized are used: the CT traces cover the attached-flow
branch, while the CP traces also run into the stalled region at low V/nD, which a strip
model with a simple post-stall CL/CD is not meant to reproduce.

    python validation/naca747/report.py        # prints the tables, writes summary.json
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
res = json.loads((HERE / "results.json").read_text())
data = json.loads((HERE / "tests.json").read_text())
C = {(c["rotation"], c["quantity"], c["beta_front"]): c for c in data["curves"]}


def meas(rot, q, b, J):
    c = C.get((rot, q, b))
    # rear CP curves (figure 12) have gaps where the dashes merge into the solid line
    if c and c["J"][0] - 1e-9 <= J <= c["J"][-1] + 1e-9 and np.min(np.abs(np.array(c["J"]) - J)) < 0.03:
        return float(np.interp(J, c["J"], c["value"]))
    return None


def table():
    pts = {}
    for r in res:
        rot = "dual" if r["kind"] == "dual" else "single"
        t, p = meas(rot, "CT", r["beta"], r["J"]), meas(rot, "CP", r["beta"], r["J"])
        if t is None or p is None or not r["converged"]:
            continue
        pts.setdefault((r["kind"], r["beta"]), []).append(
            dict(J=r["J"], CT=r["CT"], CP=r["CP"], eta=r["eta"], CTm=t, CPm=p, etam=t * r["J"] / p))
    return pts


def main():
    pts = table()
    out = {}
    print("model / measured - 1 over the attached-flow range (mean, worst); peak efficiency")
    print("kind    beta  n   CT              CP              eta model-meas")
    for (kind, b), ps in sorted(pts.items()):
        ct = np.array([p["CT"] / p["CTm"] - 1 for p in ps])
        cp = np.array([p["CP"] / p["CPm"] - 1 for p in ps])
        de = max(p["eta"] for p in ps) - max(p["etam"] for p in ps)
        w = lambda e: f"{100*e.mean():+5.1f} ({100*e[np.argmax(abs(e))]:+5.1f})"
        print(f"{kind:7s} {b:4.0f} {len(ps):2d}  {w(ct)}   {w(cp)}   {de:+.3f}")
        out[f"{kind}_{b:.0f}"] = dict(points=ps, CT_err=ct.tolist(), CP_err=cp.tolist(), peak_eta_diff=de)
    # interaction: dual relative to single at the same front blade angle and J
    print("\ndual / single at equal beta_front and J (the counter-rotation effect)")
    print("beta  J     CT ratio meas/model   CP ratio meas/model   eta diff meas/model")
    inter = []
    for b in sorted({k[1] for k in pts}):
        s = {p["J"]: p for p in pts.get(("single", b), [])}
        for p in pts.get(("dual", b), []):
            q = s.get(p["J"])
            if q is None:
                continue
            row = dict(beta=b, J=p["J"], ct_m=p["CTm"] / q["CTm"], ct_p=p["CT"] / q["CT"],
                       cp_m=p["CPm"] / q["CPm"], cp_p=p["CP"] / q["CP"],
                       de_m=p["etam"] - q["etam"], de_p=p["eta"] - q["eta"])
            inter.append(row)
            print(f"{b:4.0f} {p['J']:.2f}   {row['ct_m']:.3f} / {row['ct_p']:.3f}         "
                  f"{row['cp_m']:.3f} / {row['cp_p']:.3f}         {row['de_m']:+.3f} / {row['de_p']:+.3f}")
    out["interaction"] = inter
    # front and rear separately (figure 12), dual, over the attached-flow range
    print("\ndual, front and rear CP separately (figure 12), attached flow: model / measured")
    print("beta  J     front            rear             rear/front")
    split = []
    for r in res:
        if r["kind"] != "dual" or not r["converged"]:
            continue
        b, J = r["beta"], r["J"]
        f, rr = meas("dual", "CP_front", b, J), meas("dual", "CP_rear", b, J)
        if f is None or rr is None or meas("dual", "CT", b, J) is None:
            continue
        row = dict(beta=b, J=J, front=r["CP_front"] / f, rear=r["CP_rear"] / rr,
                   ratio_p=r["CP_rear"] / r["CP_front"], ratio_m=rr / f)
        split.append(row)
        print(f"{b:4.0f} {J:.2f}   {r['CP_front']:.4f}/{f:.4f}  {r['CP_rear']:.4f}/{rr:.4f}  "
              f"{row['ratio_p']:.3f}/{row['ratio_m']:.3f}")
    out["split"] = split
    (HERE / "summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
