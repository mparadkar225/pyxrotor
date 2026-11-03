"""Digitize NACA 747 figure 12: front (solid) and rear (dashed) CP, separated by connected
components (dashes are short pieces, solid curves long ones)."""
import copy
import json
import numpy as np
from scipy import ndimage as nd
from charts import chart, D
from digitize import crossings

c = chart("fig12")
x0, y0, x1, y1 = c.box
d = np.zeros_like(c.dark)
d[y0:y1, x0:x1] = c.img[y0:y1, x0:x1] < 135
lab, n = nd.label(d, structure=np.ones((3, 3)))
ext = np.array([max(s[0].stop - s[0].start, s[1].stop - s[1].start) for s in nd.find_objects(lab)])
solid = np.isin(lab, np.flatnonzero(ext >= 70) + 1)
dash = d & ~solid & np.isin(lab, np.flatnonzero(ext >= 12) + 1)

def sub(mask):
    k = copy.copy(c)
    k.ink = mask
    return k

cs, cd = sub(solid), sub(dash | solid)   # dashed traced on dash+solid: dashes merge into solid at crossings
cd_only = sub(dash)
BETAS = [20, 25, 30, 35, 40, 45, 50]

def seeds(k, y):
    xs = crossings(k, y, min_sep_px=6)
    return xs

if __name__ == "__main__":
    import sys
    y = float(sys.argv[1]) if len(sys.argv) > 1 else 0.06
    print("solid", np.round(seeds(cs, y), 3))
    print("dash ", np.round(seeds(cd_only, y), 3))

def column(k, J, half=3):
    xp = int(round(float(k.to_px(J, 0)[0])))
    col = k.ink[y0:y1, xp - half:xp + half + 1].any(axis=1)
    idx = np.flatnonzero(col) + y0
    if not idx.size:
        return []
    g = np.split(idx, np.flatnonzero(np.diff(idx) > 3) + 1)
    return sorted(round(float(k.to_val(0, x.mean())[1]), 4) for x in g)

SEEDS = {"front": {20: (0.36, 0.0548), 25: (0.36, 0.0894), 30: (0.36, 0.1339), 35: (0.36, 0.1783),
                   40: (0.36, 0.224), 45: (0.36, 0.2706), 50: (0.36, 0.3087)},
         "rear": {20: (0.36, 0.0586), 25: (0.36, 0.0945), 30: (0.36, 0.1462), 35: (0.36, 0.2029),
                  40: (0.36, 0.2533), 45: (1.0, 0.2662), 50: (0.36, 0.3916)}}

def trace_all(window=9, y_min=0.0):
    out = {}
    for rot, k, gap in (("front", cs, 60), ("rear", cd_only, 400)):
        for b, seed in SEEDS[rot].items():
            x, y = k.trace(seed, (0.25, 4.0), window=window, max_gap=gap)
            m = y >= y_min
            out[(rot, b)] = (x[m], y[m])
    return out


def _runs(line, lo, hi):
    lo, hi = max(int(lo), 0), int(hi) + 1
    idx = np.flatnonzero(line[lo:hi])
    if not idx.size:
        return []
    return [g.mean() + lo for g in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)]


def trace_rear(seed, J_range, window=7, look=14, step=2):
    """Follow a dashed curve: dash pixels first; through a dash gap coast (if dashes resume
    within `look` px ahead); where the dash has merged into the solid curve (touching it),
    follow the solid pixels nearest the prediction."""
    sx, sy = (float(v) for v in c.to_px(*seed))
    lo, hi = sorted(float(c.to_px(v, 0)[0]) for v in J_range)
    out = []
    for direction in (1, -1):
        x, yv, slope, n = sx, sy, 0.0, 0
        while lo <= x <= hi:
            pred = yv + slope * step * direction if n else yv
            p = int(round(x))
            dr = [r for r in _runs(dash[:, p], pred - window, pred + window)]
            hit = None
            if dr:
                hit = min(dr, key=lambda r: abs(r - pred))
            else:
                ahead = any(_runs(dash[:, int(round(x + k * direction))], pred + slope * k * direction - window,
                                  pred + slope * k * direction + window) for k in range(2, look, 2))
                if not ahead:
                    sr = _runs(solid[:, p], pred - window, pred + window)
                    if sr:
                        hit = min(sr, key=lambda r: abs(r - pred))
            if hit is not None:
                if n:
                    new = (hit - yv) / (step * direction)
                    slope = 0.6 * slope + 0.4 * new if n > 2 else new
                yv = hit
                n += 1
                out.append((x, yv))
            else:
                yv = pred
            x += step * direction
    a = np.array(sorted(set(out)))
    xv, yvv = c.to_val(a[:, 0], a[:, 1])
    return xv, yvv


def fronts(y_min=0.035):
    t = trace_all()
    out = {}
    for b in BETAS:
        x, y = t[("front", b)]
        m = y >= y_min
        x, y = x[m], y[m]
        o = np.argsort(x); x, y = x[o], y[o]
        # drop a hop onto the next curve: past J = 1 the curve only descends
        for i in range(1, len(x)):
            if x[i] > 1.0 and y[i] > y[i - 1] + 0.004:
                x, y = x[:i], y[:i]
                break
        out[b] = (x, y)
    return out


def rear_points(total, front, tol=0.012, step=0.01):
    """Dashed-pixel heights near (total - front) at each J; None where no dash is visible."""
    J = np.arange(max(total[0][0], front[0][0]), min(total[0][-1], front[0][-1]) + 1e-9, step)
    pts = []
    for j in J:
        est = np.interp(j, *total) - np.interp(j, *front)
        xp = int(round(float(c.to_px(j, 0)[0])))
        col = dash[y0:y1, xp - 1:xp + 2].any(axis=1)
        idx = np.flatnonzero(col) + y0
        if not idx.size:
            continue
        g = np.split(idx, np.flatnonzero(np.diff(idx) > 2) + 1)
        vals = [float(c.to_val(0, q.mean())[1]) for q in g]
        v = min(vals, key=lambda v: abs(v - est))
        if abs(v - est) < tol:
            pts.append((float(j), v, est))
    return np.array(pts)
