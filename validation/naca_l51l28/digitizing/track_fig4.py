import sys, json
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter, minimum_filter

def track(img, seed, xlo, xhi, win=10, step=2, slope_n=6):
    """slope_n: points used to extrapolate (more = holds direction across grid lines)."""
    """Follow a dark curve from seed (x, y) to xlo and xhi, darkest pixel near the
    extrapolated position, on a lightly smoothed image."""
    out = {seed[0]: seed[1]}
    for d in (+1, -1):
        xs, ys = [seed[0]], [float(seed[1])]
        x = seed[0]
        while True:
            x += d * step
            if x < xlo or x > xhi:
                break
            if len(xs) >= 3:
                k = min(slope_n, len(xs))
                p = np.polyfit(xs[-k:], ys[-k:], 1)
                yp = np.polyval(p, x)
            else:
                yp = ys[-1]
            w = win if len(xs) >= 3 else 3 * win      # wider until a slope is known
            lo, hi = int(yp - w), int(yp + w) + 1
            col = img[lo:hi, x]
            y = lo + int(np.argmin(col))
            if col.min() > 140:     # lost the ink
                break
            xs.append(x); ys.append(float(y)); out[x] = float(y)
    return dict(sorted(out.items()))

if __name__ == "__main__":
    cfg = json.load(open(sys.argv[1]))
    im = np.array(Image.open(cfg["image"]).convert("L")).astype(float)
    if cfg.get("remove_grid"):
        # Grid lines run across the whole plot: rows (columns) whose mean over the plot box
        # is darker than their neighbours'. Replace each such line by the mean of the rows
        # (columns) just outside it, which keeps the curves crossing it.
        from scipy.signal import find_peaks
        x0, y0, x1, y1 = cfg["plot_box"]
        box = im[y0:y1, x0:x1]
        for axis in (1, 0):
            prof = box.mean(axis=axis)
            pk, props = find_peaks(-prof, prominence=cfg.get("grid_prominence", 3), width=(0, 6))
            for p, w in zip(pk, props["widths"]):
                h = int(np.ceil(w / 2)) + 1
                if axis == 1:
                    a, b = y0 + p - h, y0 + p + h
                    im[a:b + 1, x0:x1] = 0.5 * (im[a - 1, x0:x1] + im[b + 1, x0:x1])[None, :]
                else:
                    a, b = x0 + p - h, x0 + p + h
                    im[y0:y1, a:b + 1] = 0.5 * (im[y0:y1, a - 1] + im[y0:y1, b + 1])[:, None]
        Image.fromarray(np.clip(im, 0, 255).astype("uint8")).save(cfg["overlay"].replace(".png", "_nogrid.png"))
    sm = gaussian_filter(im, 1.2)
    res = {}
    ov = Image.open(cfg["image"]).convert("RGB")
    dr = ImageDraw.Draw(ov)
    for name, c in cfg["curves"].items():
        t = track(sm, tuple(c["seed"]), c["xlo"], c["xhi"], c.get("win", 10), slope_n=c.get("slope_n", 6))
        res[name] = t
        pts = list(t.items())
        dr.line(pts, fill=(255, 0, 0), width=3)
    ov.save(cfg["overlay"])
    json.dump({k: {str(x): y for x, y in v.items()} for k, v in res.items()}, open(cfg["out"], "w"))
