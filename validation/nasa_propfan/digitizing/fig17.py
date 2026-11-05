"""CR-3505 figure 17: SR-3 aerodynamic sweep (solid line), traced on the rotated page with
the grid removed by a 4x4 opening. Calibration: bilinear between the frame corners, from
the fitted grid lines (the chart is skewed 18 px across and 8-9 px down)."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as nd
from dig2 import follow, runs
from fig16 import NoGrid

C = {(0, 0): (600.0, 427.6), (1, 0): (2928.0, 406.4), (0, 1): (610.0, 1888.5), (1, 1): (2919.0, 1867.4)}
X0, X1, Y0, Y1 = 0.2, 1.0, 60.0, -40.0     # frame values: left, right, top, bottom

def pix(xv, yv):
    u = (xv - X0) / (X1 - X0); v = (yv - Y0) / (Y1 - Y0)
    p = [(1 - u) * (1 - v) * np.array(C[0, 0]) + u * (1 - v) * np.array(C[1, 0]) +
         (1 - u) * v * np.array(C[0, 1]) + u * v * np.array(C[1, 1])][0]
    return float(p[0]), float(p[1])

def val(X, Y):
    u, v = 0.5, 0.5
    for _ in range(30):
        x, y = pix(X0 + u * (X1 - X0), Y0 + v * (Y1 - Y0))
        e = 1e-4
        xu, yu = pix(X0 + (u + e) * (X1 - X0), Y0 + v * (Y1 - Y0))
        xv, yv = pix(X0 + u * (X1 - X0), Y0 + (v + e) * (Y1 - Y0))
        Jm = np.array([[(xu - x) / e, (xv - x) / e], [(yu - y) / e, (yv - y) / e]])
        du, dv = np.linalg.solve(Jm, [X - x, Y - y]); u += du; v += dv
    return X0 + u * (X1 - X0), Y0 + v * (Y1 - Y0)

if __name__ == "__main__":
    im = np.array(Image.open(sys.argv[1]).convert('L').rotate(-90, expand=True)).astype(int)
    img = np.where(nd.binary_opening(im < 120, np.ones((4, 4))), 0, 255)
    # The dense grid partly survives the opening, so instead of following column by column,
    # take the connected ink component through a seed on the curve, and per column the
    # middle of its short vertical runs, then drop points off a local quadratic.
    ink = img < 128
    lab, n = nd.label(ink, np.ones((3, 3)))
    X, Y = pix(0.812, 42.4); X = int(X)
    ys = np.where(ink[int(Y) - 20:int(Y) + 21, X])[0] + int(Y) - 20
    k = lab[int(ys[np.argmin(abs(ys - Y))]), X]
    comp = lab == k
    tr = {}
    for x in range(int(pix(0.2, 0)[0]) + 4, int(pix(1.0, 0)[0]) - 4, 3):
        yy = np.where(comp[:, x])[0]
        if not len(yy): continue
        g = np.split(yy, np.where(np.diff(yy) > 1)[0] + 1)
        g = [a for a in g if len(a) <= 10]
        if len(g) == 1: tr[x] = float(np.mean(g[0]))
    xs = np.array(sorted(tr)); yv = np.array([tr[x] for x in xs])
    keep = []
    for i, x in enumerate(xs):
        m = (abs(xs - x) < 90) & (xs != x)
        p = np.polyfit(xs[m], yv[m], 2); keep.append(abs(np.polyval(p, x) - yv[i]) < 3)
    tr = {int(x): float(y) for x, y, kk in zip(xs, yv, keep) if kk}
    pts = [list(map(float, val(x, y))) for x, y in sorted(tr.items())]
    json.dump({"aero_sweep_deg": pts}, open('fig17_trace.json', 'w'))
    o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
    for x, y in tr.items(): d.ellipse([x - 2, y - 2, x + 2, y + 2], fill=(255, 0, 0))
    o.crop((550, 350, 3000, 1950)).resize((1225, 800)).save('ov_fig17.png')
    for r in (0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0):
        print(r, round(float(np.interp(r, *np.array(pts).T)), 2))
    print(pts[0][0], pts[-1][0])
