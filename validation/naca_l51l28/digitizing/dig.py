"""Digitize gridded NACA charts: grid detection, off-grid run finding, curve following."""
import numpy as np, json, sys
from PIL import Image, ImageDraw

def load(p):
    return np.array(Image.open(p).convert('L')).astype(int)

def grid(im, thr=130, frac=0.5):
    d = im < thr
    rows = np.where(d.mean(axis=1) > frac * d.mean(axis=1).max())[0]
    cols = np.where(d.mean(axis=0) > frac * d.mean(axis=0).max())[0]
    def grp(a):
        g = []
        for v in a:
            if g and v - g[-1][-1] <= 3: g[-1].append(v)
            else: g.append([v])
        return [(float(np.mean(x)), len(x)) for x in g]
    return grp(rows), grp(cols)

def lines(im, axis, lo, hi, thr=130, minfrac=0.6):
    """Grid lines within the frame: rows (axis=1) whose dark fraction over [lo,hi] columns."""
    d = im < thr
    if axis == 1: f = d[:, lo:hi].mean(axis=1)
    else: f = d[lo:hi, :].mean(axis=0)
    idx = np.where(f > minfrac)[0]
    g = []
    for v in idx:
        if g and v - g[-1][-1] <= 3: g[-1].append(v)
        else: g.append([v])
    return np.array([np.mean(x) for x in g]), [ (min(x), max(x)) for x in g]

def runs(col, thr=130):
    ys = np.where(col < thr)[0]
    g = []
    for y in ys:
        if g and y - g[-1][-1] <= 2: g[-1].append(y)
        else: g.append([y])
    return g

def follow(im, seed, xlo, xhi, gridrows, gridcols, step=3, tol=6, nfit=15, maxmiss=12, maxlen=9, thr=130, gpad=4):
    out = {}
    for d in (+1, -1):
        xs, ys = [seed[0]], [float(seed[1])]
        x = seed[0]; miss = 0
        while True:
            x += d * step
            if x < xlo or x > xhi: break
            if np.min(np.abs(gridcols - x)) <= gpad:   # vertical grid line: skip column
                continue
            k = min(nfit, len(xs))
            if k >= 4:
                p = np.polyfit(xs[-k:], ys[-k:], 1 if k < 8 else 2); yp = np.polyval(p, x)
            else: yp = ys[-1]
            lo = int(yp - 3 * tol)
            g = runs(im[lo:int(yp + 3 * tol) + 1, x], thr)
            c = [lo + np.mean(a) for a in g if len(a) <= maxlen]
            c = [v for v in c if np.min(np.abs(gridrows - v)) > gpad]
            c = [v for v in c if abs(v - yp) <= (tol if len(xs) >= 4 else 2 * tol)]
            if not c:
                miss += 1
                if miss > maxmiss: break
                continue
            y = min(c, key=lambda v: abs(v - yp)); miss = 0
            xs.append(x); ys.append(y); out[x] = y
    out[seed[0]] = float(seed[1])
    return dict(sorted(out.items()))

def frame(im, nrows, ncols, thr=130):
    """The chart's grid: nrows rows and ncols columns, evenly spaced (to a few px)."""
    out = []
    for axis, n in ((1, nrows), (0, ncols)):
        d = (im < thr)
        f = d.mean(axis=axis)
        cand = []
        idx = np.where(f > 0.45 * f.max())[0]
        g = []
        for v in idx:
            if g and v - g[-1][-1] <= 3: g[-1].append(v)
            else: g.append([v])
        cand = np.array([np.mean(x) for x in g])
        best = None
        for i in range(len(cand)):
            for j in range(i + 1, len(cand)):
                s = (cand[j] - cand[i])
                for kk in range(1, n):
                    sp = s / kk
                    if sp < 40: break
                    pred = cand[i] + sp * np.arange(n)
                    m = [np.min(np.abs(cand - p)) for p in pred]
                    sc = sum(v < 4 for v in m)
                    if best is None or sc > best[0]:
                        best = (sc, [cand[np.argmin(np.abs(cand - p))] if mm < 4 else p for p, mm in zip(pred, m)])
        out.append(np.array(best[1]))
    return out  # rows (top->bottom), cols (left->right)
