"""CR-3505 figure 16 (SR-3 blade characteristics): trace b/D, t/b, CLD, dbeta on the page
image rotated upright, with the grid removed by a 4x4 binary opening (the grid lines are
1-2 px, the curves 4-6 px)."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as nd
from dig2 import follow, runs

class NoGrid:
    def on_vline(self, *a): return False
    def on_hline(self, *a): return False

def val(X, Y):
    """Pixel -> (r/R, left-axis value) for the rotated page: rows tilt 9 px / 1400 px,
    columns 8 px / 850 px (fitted to the solid grid lines)."""
    Yc = Y - 0.00643 * (X - 2000.0)
    Xc = X + 0.0094 * (Y - 1200.0)
    return 0.2 + (Xc - 1066.8) / (2915.8 - 1066.8) * 0.8, (1800.5 - Yc) / (1800.5 - 647.5) * 0.20

def pix(xv, yv):
    X, Y = 2000.0, 1200.0
    for _ in range(5):
        Xc = 1066.8 + (xv - 0.2) / 0.8 * (2915.8 - 1066.8)
        Yc = 1800.5 - yv / 0.20 * (1800.5 - 647.5)
        X = Xc - 0.0094 * (Y - 1200.0); Y = Yc + 0.00643 * (X - 2000.0)
    return X, Y

if __name__ == "__main__":
    src = sys.argv[1]
    im = np.array(Image.open(src).convert('L').rotate(-90, expand=True)).astype(int)
    ink = nd.binary_opening(im < 120, np.ones((4, 4)))
    img = np.where(ink, 0, 255)
    # (seed r/R, seed value, r/R range, tolerance px); dbeta and t/b cross at 0.905, so their
    # outer parts are traced from their own seeds
    seeds = {"b_over_D": [(0.55, 0.205, 0.2, 1.0, 5)], "CLD": [(0.80, 0.127, 0.2, 1.0, 5)],
             "dbeta": [(0.70, 0.047, 0.2, 0.90, 3), (0.97, 0.011, 0.915, 1.0, 3)],
             "t_over_b": [(0.60, 0.034, 0.2, 0.89, 3), (0.97, 0.018, 0.92, 1.0, 3)]}
    out = {}
    o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
    cols = {"b_over_D": (255, 0, 0), "CLD": (0, 150, 0), "dbeta": (0, 0, 255), "t_over_b": (200, 0, 200)}
    x0, _ = pix(0.2, 0.1); x1, _ = pix(1.0, 0.1)
    for k, parts in seeds.items():
      tr = {}
      for xv, yv, lo, hi, tol in parts:
        X, Y = pix(xv, yv)
        X = int(round(X))
        g = runs(img[int(Y) - 25:int(Y) + 26, X], 128)
        c = [int(Y) - 25 + np.mean(a) for a in g]
        Y = min(c, key=lambda v: abs(v - Y))
        tr.update(follow(img, NoGrid(), (X, Y), max(int(pix(lo, yv)[0]), int(x0) + 3),
                         min(int(pix(hi, yv)[0]), int(x1) - 3), step=3, tol=tol, nfit=30,
                         maxlen=14, thr=128, maxmiss=6))
      if True:
        out[k] = [list(map(float, val(x, y))) for x, y in sorted(tr.items())]
        for x, y in tr.items(): d.ellipse([x - 2, y - 2, x + 2, y + 2], fill=cols[k])
        print(k, len(tr), round(out[k][0][0], 3), round(out[k][-1][0], 3))
    json.dump(out, open('fig16_traces.json', 'w'))
    o.crop((1000, 350, 3000, 1900)).resize((1000, 775)).save('ov_fig16.png')


def steep_root(src, out_json='fig16_traces.json'):
    """t/b inboard of 0.31R is nearly vertical: trace it row by row (r/R as a function of
    t/b) between t/b 0.10 and 0.235, the ink run nearest the previous row's."""
    im = np.array(Image.open(src).convert('L').rotate(-90, expand=True)).astype(int)
    ink = nd.binary_opening(im < 120, np.ones((4, 4)))
    pts = []
    X = pix(0.31, 0.096)[0]
    for yv in np.arange(0.098, 0.236, 0.002):
        Xp, Y = pix(0.27, yv)
        row = ink[int(round(Y)), :]
        xs = np.where(row[int(X) - 40:int(X) + 10])[0] + int(X) - 40
        if not len(xs): break
        g = np.split(xs, np.where(np.diff(xs) > 2)[0] + 1)
        c = [np.mean(a) for a in g]
        X = min(c, key=lambda v: abs(v - X))
        pts.append(val(X, Y))
    d = json.load(open(out_json))
    d["t_over_b_root"] = [list(map(float, p)) for p in pts]
    json.dump(d, open(out_json, 'w'))
    return pts
