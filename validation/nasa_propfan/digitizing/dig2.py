"""Gridded-chart digitizing with skewed grids: each grid line is fitted as a straight line,
pixel -> value by interpolating between the lines at the point (piecewise-linear in both)."""
import numpy as np, json
from PIL import Image, ImageDraw
from dig import load, frame, runs

class Grid:
    def __init__(self, im, nrows, ncols, xv, yv, thr=130):
        r, c = frame(im, nrows, ncols, thr)
        self.xv, self.yv = xv, yv      # (first, step) values of col 0 / row 0 (top)
        d = im < thr
        def centroid(prof, x0):
            w = np.arange(int(x0) - 12, int(x0) + 13); pr = prof[w]
            return (w * pr).sum() / pr.sum()
        V = []
        for x0 in c:   # vertical lines, x = a + b y, from 6 bands
            pts = []
            for y0 in np.linspace(r[0] + 15, r[-1] - 115, 8):
                prof = d[int(y0):int(y0) + 100, :].mean(axis=0)
                pts.append((y0 + 50, centroid(prof, x0)))
            pts = np.array(pts); V.append(np.polyfit(pts[:, 0], pts[:, 1], 1))
        H = []
        for y0 in r:
            pts = []
            for x0 in np.linspace(c[0] + 15, c[-1] - 115, 8):
                prof = d[:, int(x0):int(x0) + 100].mean(axis=1)
                pts.append((x0 + 50, centroid(prof, y0)))
            pts = np.array(pts); H.append(np.polyfit(pts[:, 0], pts[:, 1], 1))
        self.V, self.H = np.array(V), np.array(H)
        self.resid = None
    def xs_at(self, y): return self.V[:, 0] * y + self.V[:, 1]
    def ys_at(self, x): return self.H[:, 0] * x + self.H[:, 1]
    def value(self, x, y):
        k = np.interp(x, self.xs_at(y), np.arange(len(self.V)))
        j = np.interp(y, self.ys_at(x), np.arange(len(self.H)))
        return self.xv[0] + k * self.xv[1], self.yv[0] + j * self.yv[1]
    def pixel(self, vx, vy):
        k = (vx - self.xv[0]) / self.xv[1]; j = (vy - self.yv[0]) / self.yv[1]
        x, y = np.interp(k, np.arange(len(self.V)), self.xs_at(self.ys_at(1000).mean())), 0
        for _ in range(4):
            y = np.interp(j, np.arange(len(self.H)), self.ys_at(x))
            x = np.interp(k, np.arange(len(self.V)), self.xs_at(y))
        return x, y
    def on_vline(self, x, y, pad): return np.min(np.abs(self.xs_at(y) - x)) <= pad
    def on_hline(self, x, y, pad): return np.min(np.abs(self.ys_at(x) - y)) <= pad

def follow(im, G, seed, xlo, xhi, step=3, tol=5, nfit=18, maxmiss=15, maxlen=9, thr=140, gpad=4):
    out = {}
    for d in (+1, -1):
        xs, ys = [seed[0]], [float(seed[1])]
        x = seed[0]; miss = 0
        while True:
            x += d * step
            if x < xlo or x > xhi: break
            k = min(nfit, len(xs))
            if k >= 4:
                p = np.polyfit(xs[-k:], ys[-k:], 1 if k < 10 else 2); yp = np.polyval(p, x)
            else: yp = ys[-1]
            if G.on_vline(x, yp, gpad + 1): continue
            lo = int(yp - 3 * tol)
            g = runs(im[lo:int(yp + 3 * tol) + 1, x], thr)
            c = [lo + np.mean(a) for a in g if len(a) <= maxlen]
            c = [v for v in c if not G.on_hline(x, v, gpad)]
            c = [v for v in c if abs(v - yp) <= (tol if len(xs) >= 4 else 2 * tol)]
            if not c:
                miss += 1
                if miss > maxmiss: break
                continue
            y = min(c, key=lambda v: abs(v - yp)); miss = 0
            xs.append(x); ys.append(y); out[x] = y
    out[seed[0]] = float(seed[1])
    return dict(sorted(out.items()))

def extend(im, G, pts, direction, xlo, xhi, step=3, tol=4, nfit=20, maxmiss=10, thr=140, gpad=4, max_on=45):
    """Continue a traced curve (dict x->y) past its end, also through stretches where it
    runs along a horizontal grid line (then the line's own centre is taken)."""
    xs = sorted(pts); 
    if direction < 0: xs = xs[::-1]
    X = list(xs[-nfit:] if direction > 0 else xs[-nfit:]); X = xs[-nfit:]
    Y = [pts[x] for x in X]
    x = xs[-1]; miss = 0; on = 0; new = {}; last_off = xs[-1]
    while True:
        x += direction * step
        if x < xlo or x > xhi: break
        k = min(nfit, len(X)); p = np.polyfit(X[-k:], Y[-k:], 2 if k >= 10 else 1); yp = np.polyval(p, x)
        if G.on_vline(x, yp, gpad + 1): continue
        lo = int(yp - 3 * tol); g = runs(im[lo:int(yp + 3 * tol) + 1, x], thr)
        c = [(lo + np.mean(a), G.on_hline(x, lo + np.mean(a), gpad)) for a in g if len(a) <= 14]
        c = [(v, o) for v, o in c if abs(v - yp) <= tol]
        if not c:
            miss += 1
            if miss > maxmiss: break
            continue
        offg = [v for v, o in c if not o]
        if offg: y = min(offg, key=lambda v: abs(v - yp)); on = 0; last_off = x
        else:
            y = min([v for v, o in c], key=lambda v: abs(v - yp)); on += 1
            if on > max_on: break
        miss = 0; X.append(x); Y.append(y); new[x] = y
    if new:
        lo_ = min(last_off, xs[-1]) if direction > 0 else max(last_off, xs[-1])
        new = {k: v for k, v in new.items() if (k <= max(last_off, xs[-1]) if direction > 0 else k >= min(last_off, xs[-1]))}
    # drop a trailing on-grid stretch that ended the trace (it may be the grid line alone)
    return new
