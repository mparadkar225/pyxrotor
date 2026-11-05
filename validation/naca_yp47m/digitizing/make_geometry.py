"""geometry.json from the blade-form traces (bladeform_traces.json, written by bladeform.py).

Sampled at r/R 0.20-0.95 every 0.025; a trace is held at its end value beyond its range
(np.interp). Additions by eye:
* Curtiss 836 b/D inboard of 0.55R: the trace followed the blade-angle curve where the two
  cross (0.3R), so b/D there is read by eye from figure 1.
* Tips: the traces stop at 0.91-0.96R. The chord at 0.975R and 1.0R is 0.75 and 0.35 of
  the 0.95R chord (the figures' rounded tips), the blade angle is extrapolated along its
  last 0.05R slope.
Run from this folder: python make_geometry.py
"""
import json

import numpy as np

T = json.load(open("bladeform_traces.json"))
# table on p. 3-4 of each report (checked on the page images)
INFO = {
    "aeroproducts_H20C-162": dict(ntrs="20030063991", D_in=12 * 12 + 7, activity_factor=107),
    "curtiss_838-1C2-18R1": dict(ntrs="20030063996", D_in=13 * 12, activity_factor=120),
    "curtiss_732-1C2-0": dict(ntrs="20030064146", D_in=13 * 12, activity_factor=150),
    "curtiss_836-14C2-18R1": dict(ntrs="20030063262", D_in=13 * 12, activity_factor=98),
}
EYE_836_BD = [(0.15, 0.027), (0.20, 0.039), (0.30, 0.050), (0.40, 0.063), (0.50, 0.068)]

R = np.round(np.arange(0.20, 0.9501, 0.025), 4)
out = {}
for name, info in INFO.items():
    tr = {k: np.array(v) for k, v in T[name]["traces"].items()}
    bd = tr["b_over_D"]
    if name.startswith("curtiss_836"):
        keep = bd[:, 0] > 0.55
        bd = np.vstack([np.array(EYE_836_BD), bd[keep]])
    cl = tr.get("cl_design", tr.get("cl_design_by_eye"))
    th = tr["theta_deg"]
    b = np.interp(R, bd[:, 0], bd[:, 1])
    t = np.interp(R, tr["h_over_b"][:, 0], tr["h_over_b"][:, 1])
    c = np.interp(R, cl[:, 0], cl[:, 1])
    a = np.interp(R, th[:, 0], th[:, 1])
    slope = (a[-1] - a[-3]) / (R[-1] - R[-3])
    r = np.concatenate([R, [0.975, 1.0]])
    b = np.concatenate([b, [0.75 * b[-1], 0.35 * b[-1]]])
    t = np.concatenate([t, [t[-1], t[-1]]])
    c = np.concatenate([c, [c[-1], c[-1]]])
    a = np.concatenate([a, [a[-1] + 0.025 * slope, a[-1] + 0.05 * slope]])
    tw = a - np.interp(0.75, r, a)
    x = np.linspace(0.15, 1.0, 2001)
    af = 1e5 / 16 * np.trapezoid(np.interp(x, np.r_[0.15, r], np.r_[b[0], b]) * x ** 3, x)
    out[name] = dict(info, nblds=4, activity_factor_traced=round(float(af), 1),
                     r_over_R=[float(v) for v in r], b_over_D=[round(float(v), 5) for v in b],
                     h_over_b=[round(float(v), 5) for v in t], cl_design=[round(float(v), 4) for v in c],
                     **{"twist_deg_from_0.75R": [round(float(v), 3) for v in tw]})
    print(name, "AF traced", out[name]["activity_factor_traced"], "quoted", info["activity_factor"])
json.dump(out, open("../geometry.json", "w"), indent=1)
