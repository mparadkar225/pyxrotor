import sys, json; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from dig2 import *
from auto import curves

def same(a, b, gap_max=200, tol_ov=4, tol_gap=7):
    ax, bx = sorted(a), sorted(b)
    lo, hi = max(ax[0], bx[0]), min(ax[-1], bx[-1])
    if hi - lo >= 6:
        xs = [x for x in ax if lo <= x <= hi]
        d = [abs(a[x] - np.interp(x, bx, [b[k] for k in bx])) for x in xs]
        return np.median(d) < tol_ov
    if ax[-1] < bx[0]: L, R = a, b
    else: L, R = b, a
    lx, rx = sorted(L), sorted(R)
    gap = rx[0] - lx[-1]
    if gap > gap_max: return False
    pl = np.polyfit(lx[-25:], [L[x] for x in lx[-25:]], 1)
    pr = np.polyfit(rx[:25], [R[x] for x in rx[:25]], 1)
    xm = 0.5 * (lx[-1] + rx[0])
    return abs(np.polyval(pl, xm) - np.polyval(pr, xm)) < tol_gap and abs(pl[0] - pr[0]) < 0.15

def merged(p):
    im, G, cs = curves(p)
    groups = []
    for c in sorted(cs, key=len, reverse=True):
        for g in groups:
            if any(same(c, h) for h in g):
                g.append(c); break
        else: groups.append([c])
    out = []
    for g in groups:
        pts = {}
        for c in g:
            for x, y in c.items(): pts.setdefault(x, []).append(y)
        out.append({x: float(np.median(v)) for x, v in sorted(pts.items())})
    return im, G, out
