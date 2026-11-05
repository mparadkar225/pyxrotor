"""TM-79124 figures 12-13: net efficiency against J at 70-100% of design power loading,
gridless line charts. Per panel: the left axis line and its ticks, the bottom axis and its
ticks give the calibration; the four curves are followed from seeds where they are apart,
ordered top to bottom as 70, 80, 90, 100% (their order at high J)."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from PIL import Image, ImageDraw
from scipy.signal import find_peaks
from dig2 import follow, runs
from fig16 import NoGrid

def calib(im, box, yvals, xvals):
    """box = (x0, y0, x1, y1) around the panel's axes. Returns (axis x, y ticks, x ticks)."""
    x0, y0, x1, y1 = box
    d = im < 128
    c = d[y0:y1, x0:x0 + 200].sum(axis=0); ax = x0 + int(np.argmax(c))
    r = d[y0:y1, ax + 40:x1].sum(axis=1); bot = y0 + int(np.argmax(r))
    t = d[y0:bot + 3, ax + 3:ax + 25].sum(axis=1)
    yk, _ = find_peaks(t, height=12, distance=25); yk = y0 + yk
    t = d[bot - 30:bot - 3, ax + 20:x1].sum(axis=0)
    xk, _ = find_peaks(t, height=15, distance=25); xk = ax + 20 + xk
    yk = yk.tolist()
    if not yk or abs(yk[-1] - bot) > 4: yk.append(bot)       # the axis is the lowest tick
    yk = yk[-len(yvals):]
    # x ticks: keep those on an even spacing from the axis (labels near the axis add peaks)
    sp = (xk[-1] - ax) / (len(xvals) - 1)
    xs = [ax]
    for k in range(1, len(xvals)):
        near = [v for v in xk if abs(v - (ax + k * sp)) < 0.2 * sp]
        xs.append(min(near, key=lambda v: abs(v - (ax + k * sp))) if near else ax + k * sp)
    return ax, bot, np.array(yk, float), np.array(xs, float)

def make_maps(yk, xk, yvals, xvals):
    yv = np.array(sorted(yvals, reverse=True), float)
    py = np.polyfit(yk, yv, 1); px = np.polyfit(xk, xvals, 1)
    return (lambda X, Y: (float(np.polyval(px, X)), float(np.polyval(py, Y))),
            lambda J, e: (float((J - px[1]) / px[0]), float((e - py[1]) / py[0])),
            np.abs(np.polyval(py, yk) - yv).max(), np.abs(np.polyval(px, xk) - xvals).max())

def trace(im, box, yvals, xvals, seed_J, seeds_eta, xlim, tol=3, name='p'):
    ax, bot, yk, xk = calib(im, box, yvals, xvals)
    val, pix, ey, ex = make_maps(yk, xk, yvals, xvals)
    out = {}
    o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
    colors = [(255, 0, 0), (0, 150, 0), (0, 0, 255), (200, 0, 200)]
    for lab, eta, col in zip((70, 80, 90, 100), seeds_eta, colors):
        X, Y = pix(seed_J, eta); X = int(round(X))
        g = runs(im[int(Y) - 15:int(Y) + 16, X], 128)
        Y = min([int(Y) - 15 + np.mean(a) for a in g], key=lambda v: abs(v - Y))
        tr = follow(im, NoGrid(), (X, Y), int(pix(xlim[0], 0)[0]), int(pix(xlim[1], 0)[0]),
                    step=2, tol=tol, nfit=25, maxlen=10, thr=128, maxmiss=6)
        out[lab] = [val(x, y) for x, y in sorted(tr.items())]
        for x, y in tr.items(): d.ellipse([x - 1, y - 1, x + 1, y + 1], fill=col)
    x0, y0, x1, y1 = box
    o.crop((x0, y0, x1, y1)).save(f'ov_{name}.png')
    return out, (ey, ex, yk.tolist(), xk.tolist())
