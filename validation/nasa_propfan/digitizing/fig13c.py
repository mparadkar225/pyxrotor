"""TM X-73612 figure 13(c) (SR-2 blade characteristics), gridless: trace b/D, dbeta, CLD,
t/b. Calibration from the t/b axis ticks (0 at y 1178, 0.22 at 652) and the r/R ticks
(0.4-1.0 at x 1528-1811; 0.2 is the t/b axis line at 1434). The other axes share the
t/b scale: b/D = 0.04 + 0.5 y, CLD = -0.3 + 5 y, dbeta = -12 + 200 y (deg)."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from PIL import Image, ImageDraw
from dig2 import follow, runs
from fig16 import NoGrid

def val(X, Y): return 0.2 + (X - 1434.0) / (1811.0 - 1434.0) * 0.8, (1178.0 - Y) / (1178.0 - 652.0) * 0.22
def pix(x, y): return 1434.0 + (x - 0.2) / 0.8 * (1811.0 - 1434.0), 1178.0 - y / 0.22 * (1178.0 - 652.0)

if __name__ == "__main__":
    im = np.array(Image.open(sys.argv[1]).convert('L')).astype(int)
    seeds = {"b_over_D_axis": [(0.6, 0.216, 0.21, 0.999, 4)],
             "dbeta_axis": [(0.45, 0.112, 0.21, 0.999, 3)],
             "CLD_axis": [(0.8, 0.081, 0.21, 0.999, 4)],
             "t_over_b": [(0.6, 0.033, 0.21, 0.999, 3)]}
    out = {}
    o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
    cols = [(255, 0, 0), (0, 0, 255), (0, 150, 0), (200, 0, 200)]
    for (k, parts), col in zip(seeds.items(), cols):
        tr = {}
        for xv, yv, lo, hi, tol in parts:
            X, Y = pix(xv, yv); X = int(round(X))
            g = runs(im[int(Y) - 15:int(Y) + 16, X], 128)
            Y = min([int(Y) - 15 + np.mean(a) for a in g], key=lambda v: abs(v - Y))
            tr.update(follow(im, NoGrid(), (X, Y), int(pix(lo, 0)[0]), int(pix(hi, 0)[0]),
                             step=2, tol=tol, nfit=25, maxlen=18, thr=128, maxmiss=8))
        out[k] = [list(map(float, val(x, y))) for x, y in sorted(tr.items())]
        for x, y in tr.items(): d.ellipse([x - 1, y - 1, x + 1, y + 1], fill=col)
        print(k, len(tr), round(out[k][0][0], 3), round(out[k][-1][0], 3))
    # steep parts, row by row: the b/D drop at the tip, t/b inboard of 0.37R
    def rows(x0, y_from, y_to, dy, win=(-30, 30)):
        pts, X = [], pix(x0, y_from)[0]
        for yv in np.arange(y_from, y_to, dy):
            Y = int(round(pix(0.5, yv)[1]))
            xs = np.where(im[Y, int(X) + win[0]:int(X) + win[1]] < 128)[0] + int(X) + win[0]
            if not len(xs): break
            g = np.split(xs, np.where(np.diff(xs) > 2)[0] + 1)
            X = min([np.mean(a) for a in g if len(a) < 14] or [X], key=lambda v: abs(v - X))
            pts.append(val(X, Y)); d.ellipse([X - 1, Y - 1, X + 1, Y + 1], fill=(255, 150, 0))
        return [list(map(float, p)) for p in pts]
    out["b_over_D_axis_tip"] = rows(0.925, 0.198, 0.04, -0.002)
    out["t_over_b_root"] = rows(0.366, 0.063, 0.215, 0.002, (-25, 8))
    json.dump(out, open('fig13c_traces.json', 'w'))
    o.crop((1300, 600, 1850, 1250)).resize((1100, 1300)).save('ov_fig13c.png')
