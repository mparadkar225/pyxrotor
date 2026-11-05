"""YP-47M reports' figure 1 (blade-form curves): trace b/D, h/b, blade angle and design
section lift coefficient on a coarse grid. Curves are followed column by column,
skipping grid columns and runs on grid rows (dig2.follow), then extended through stretches
that run along a grid line (dig2.extend). Calibration: the detected grid rows and columns,
interpolated piecewise linearly; each curve's scale is value = a + b * k, k the row count
above the bottom frame line, read from the figure's axis labels."""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as nd
from scipy.signal import find_peaks
from dig2 import follow, runs, extend


class LineGrid:
    """Grid lines at fixed pixel positions (the pages here are skewed by < 2 px)."""
    def __init__(self, rows, cols):
        self.rows, self.cols = np.asarray(rows), np.asarray(cols)
    def xs_at(self, y): return self.cols
    def ys_at(self, x): return self.rows
    def on_vline(self, x, y, pad): return np.min(np.abs(self.cols - x)) <= pad
    def on_hline(self, x, y, pad): return np.min(np.abs(self.rows - y)) <= pad

FIGS = {
 "aeroproducts_H20C-162": dict(page="20030063991/pages/p18-18.png", bottom=2200, cols=(577, 1780),
     scales={"b_over_D": (0, 0.01), "h_over_b": (0, 0.05), "theta_deg": (0, 5.0), "cl_design": (-0.7, 0.1)},
     seeds={"b_over_D": [(0.65, 0.0834)], "h_over_b": [(0.65, 0.0695)], "theta_deg": [(0.65, 23.8)]},
     manual={"cl_design": [(0.29, 0.0), (0.30, 0.05), (0.31, 0.30), (0.33, 0.45), (0.36, 0.50), (1.0, 0.50)]}),
 "curtiss_838-1C2-18R1": dict(page="20030063996/pages/p14-14.png", bottom=2402, cols=(643, 1877),
     scales={"b_over_D": (0, 0.01), "h_over_b": (0, 0.1), "theta_deg": (0, 5.0), "cl_design": (-0.8, 0.1)},
     seeds={"b_over_D": [(0.75, 0.0754), (0.25, 0.0429)], "h_over_b": [(0.75, 0.054), (0.25, 0.312)],
            "theta_deg": [(0.75, 26.25), (0.35, 46.1)]},
     manual={"cl_design": [(0.35, 0.75), (0.40, 0.67), (0.45, 0.61), (0.48, 0.60), (1.0, 0.60)]}),
 "curtiss_732-1C2-0": dict(page="20030064146/pages/p9-09.png", bottom=2107, cols=(695, 1928),
     scales={"b_over_D": (0.02, 0.01), "h_over_b": (0, 0.1), "theta_deg": (0, 5.0), "cl_design": (-0.5, 0.1)},
     seeds={"b_over_D": [(0.45, 0.1065), (0.55, 0.1123), (0.85, 0.106, 0.8, 1.0), (0.95, 0.0848, 0.8, 1.0)], "h_over_b": [(0.45, 0.08)],
            "theta_deg": [(0.45, 40.95), (0.85, 23.8)], "cl_design": [(0.75, 0.622), (0.25, 0.438)]}),
 "hamilton_6507A-2": dict(page="20030065899/pages/p17-17.png", bottom=2211, cols=(587, 1834),
     scales={"b_over_D": (0, 0.01), "h_over_b": (0, 0.05), "theta_deg": (0, 5.0), "cl_design": (-0.5, 0.1)},
     seeds={"b_over_D": [(0.85, 0.0709)], "h_over_b": [(0.65, 0.078), (0.85, 0.0655)],
            "theta_deg": [(0.65, 23.5)], "cl_design": [(0.65, 0.461)]}),
 "curtiss_836-14C2-18R1": dict(page="20030063262/pages/p17-17.png", bottom=2119, cols=(665, 1889),
     scales={"b_over_D": (0, 0.01), "h_over_b": (0, 0.1), "theta_deg": (0, 10.0), "cl_design": (-0.4, 0.1)},
     seeds={"b_over_D": [(0.75, 0.0617, 0.45, 1.0), (0.35, 0.0466, 0.14, 0.48)], "h_over_b": [(0.75, 0.065), (0.35, 0.142)],
            "theta_deg": [(0.75, 26.2)]},
     manual={"cl_design": [(0.42, 0.60), (0.46, 0.53), (0.50, 0.50), (1.0, 0.50)]}),
}


class _NoGrid:
    def on_vline(self, *a): return False
    def on_hline(self, *a): return False


NOGRID = _NoGrid()


def calib(im, bottom, cols):
    d = im < 140
    r = d[:, cols[0] + 20:cols[1] - 20].mean(axis=1)
    pr, _ = find_peaks(r, height=0.3, distance=60)
    pr = pr[pr <= bottom + 5]
    rows = np.array(sorted(pr))[::-1].astype(float)          # bottom first: k = 0, 1, ...
    c = d[int(rows[-1]) + 20:int(rows[0]) - 20, :].mean(axis=0)
    pc, _ = find_peaks(c, height=0.3, distance=60)
    pc = pc[(pc >= cols[0] - 5) & (pc <= cols[1] + 5)].astype(float)
    return rows, pc


def trace_fig(name, cfg, out_png):
    im = np.array(Image.open(cfg["page"]).convert("L")).astype(int)
    rows, cols = calib(im, cfg["bottom"], cfg["cols"])
    ks = np.arange(len(rows)); xs = np.linspace(0, 1, len(cols))
    kpix = lambda k: float(np.interp(k, ks, rows))
    pixk = lambda y: float(np.interp(-y, -rows, ks))
    xpix = lambda x: float(np.interp(x, xs, cols))
    pixx = lambda X: float(np.interp(X, cols, xs))
    G = LineGrid(rows, cols)
    img = np.where(nd.binary_opening(im < 130, np.ones((4, 4))), 0, 255)
    o = Image.fromarray(im.astype(np.uint8)).convert("RGB"); dr = ImageDraw.Draw(o)
    out = {}
    colors = {"b_over_D": (255, 0, 0), "h_over_b": (200, 0, 200), "theta_deg": (0, 0, 255), "cl_design": (0, 150, 0)}
    for cname, seeds in cfg["seeds"].items():
        a, b = cfg["scales"][cname]
        tr = {}
        for sd in seeds:
            x, v = sd[0], sd[1]
            lo, hi = (sd[2], sd[3]) if len(sd) > 2 else (0.0, 1.0)
            X = int(round(xpix(x))); Y = kpix((v - a) / b)
            g = runs(img[int(Y) - 30:int(Y) + 31, X], 128)
            if not g:
                print(name, cname, "no ink at seed"); continue
            Y = min([int(Y) - 30 + np.mean(q) for q in g], key=lambda q: abs(q - Y))
            t = follow(img, G, (X, Y), max(int(cols[0]) + 3, int(xpix(lo))),
                       min(int(cols[-1]) - 3, int(xpix(hi))), step=3, tol=4, nfit=25, maxlen=16,
                       thr=128, maxmiss=40, gpad=4)
            tr.update(t)
        # drop stretches that run along a grid row for more than 60 px (the tracer left the
        # curve for a grid line the opening did not remove)
        xs_ = sorted(tr); on = [np.min(np.abs(rows - tr[x])) < 3.0 for x in xs_]
        i = 0
        while i < len(xs_):
            if on[i]:
                j = i
                while j + 1 < len(xs_) and on[j + 1]: j += 1
                if xs_[j] - xs_[i] > 60:
                    for x in xs_[i:j + 1]: tr.pop(x)
                i = j + 1
            else:
                i += 1
        out[cname] = [(pixx(X), a + b * pixk(Y)) for X, Y in sorted(tr.items())]
        for X, Y in tr.items(): dr.ellipse([X - 2, Y - 2, X + 2, Y + 2], fill=colors[cname])
        if tr:
            print(name, cname, len(tr), round(out[cname][0][0], 3), round(out[cname][-1][0], 3))
    for cname, pts in cfg.get("manual", {}).items():          # read by eye: on a grid row
        out[cname + "_by_eye"] = pts
        a, b = cfg["scales"][cname]
        for x, v in pts:
            X, Y = xpix(x), kpix((v - a) / b); dr.rectangle([X - 5, Y - 5, X + 5, Y + 5], outline=(255, 140, 0), width=3)
    o.crop((int(cols[0]) - 50, int(rows[-1]) - 50, int(cols[-1]) + 50, int(rows[0]) + 50)).resize(
        ((int(cols[-1] - cols[0]) + 100) // 2, (int(rows[0] - rows[-1]) + 100) // 2)).save(out_png)
    return out, {"rows": rows.tolist(), "cols": cols.tolist()}


if __name__ == "__main__":
    res = {}
    for name, cfg in FIGS.items():
        tr, cal = trace_fig(name, cfg, f"ov_bf_{name}.png")
        res[name] = {"traces": tr, "calibration": cal}
    json.dump(res, open("bladeform_traces.json", "w"))
