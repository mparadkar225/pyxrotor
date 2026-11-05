"""Write traces_fig6-8.json and static_fig6-8.csv from run_mt2.py's mt_full.json:
the traces thinned to steps of 0.01 in Mt, and the curves sampled every 0.05 in Mt where
both CT and CP were traced.

    python sample.py mt_full.json ..
"""
import csv
import json
import sys

import numpy as np

src, out = sys.argv[1], sys.argv[2].rstrip("/") + "/"
d = json.load(open(src))
FIG = {16: ('6a', '0', 'CT'), 17: ('6b', '0', 'CP'), 18: ('7a', '3', 'CT'), 19: ('7b', '3', 'CP'),
       20: ('8a', '5', 'CT'), 21: ('8b', '5', 'CP')}
# figure 8(a), 12 deg: past Mt 0.89 the curve runs along the CT = 0.105 grid line; read by eye
MANUAL = {(20, 12): [(0.92, 0.1045), (0.95, 0.1050), (0.98, 0.1047), (1.00, 0.1043), (1.03, 0.1034)]}
traces = {}
for p, curves in d.items():
    p = int(p)
    fig, cam, q = FIG[p]
    for beta, v in zip([0, 4, 8, 12, 16], curves):
        v = np.array(v)
        v = v[np.argsort(v[:, 0])]
        if (p, beta) in MANUAL:
            m = np.array(MANUAL[p, beta])
            v = np.vstack([v[v[:, 0] < m[0, 0] - 0.01], m])
        keep = [0]
        for i in range(1, len(v)):
            if v[i, 0] - v[keep[-1], 0] >= 0.01:
                keep.append(i)
        if keep[-1] != len(v) - 1:
            keep.append(len(v) - 1)
        traces.setdefault(cam, {}).setdefault(str(beta), {})[q] = \
            [[round(a, 4), round(b, 5)] for a, b in v[keep]]
json.dump({"source": "NACA RM L51L28 (Wood & Swihart 1952), figures 6-8, traced from 300 dpi page "
           "images (pages 15-20); Mt, coefficient pairs along each faired curve",
           "manual": "figure 8(a), beta 12 deg, Mt > 0.9: read by eye where the curve runs along "
           "the CT = 0.105 grid line",
           "curves": traces}, open(out + 'traces_fig6-8.json', 'w'), indent=0)
rows = []
for cam, bs in traces.items():
    for beta, qs in bs.items():
        for mt in np.round(np.arange(0.30, 1.001, 0.05), 2):
            vals = {}
            for q in ('CT', 'CP'):
                v = np.array(qs[q])
                if v[0, 0] - 0.005 <= mt <= v[-1, 0] + 0.005:
                    vals[q] = float(np.interp(mt, v[:, 0], v[:, 1]))
            if len(vals) == 2:
                rows.append(dict(design_cl=f'0.{cam}', beta75=beta, Mt=f'{mt:.2f}',
                                 CT=f"{vals['CT']:.4f}", CP=f"{vals['CP']:.4f}",
                                 CT_over_CP=f"{vals['CT'] / vals['CP']:.3f}"))
with open(out + 'static_fig6-8.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(rows[0]))
    w.writeheader()
    w.writerows(rows)
