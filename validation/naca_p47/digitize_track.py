import sys, json
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter, minimum_filter

def track(img, seed, xlo, xhi, win=10, step=2, slope_n=6):
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
            lo, hi = int(yp - win), int(yp + win) + 1
            col = img[lo:hi, x]
            y = lo + int(np.argmin(col))
            if col.min() > 140:     # lost the ink
                break
            xs.append(x); ys.append(float(y)); out[x] = float(y)
    return dict(sorted(out.items()))

if __name__ == "__main__":
    cfg = json.load(open(sys.argv[1]))
    im = np.array(Image.open(cfg["image"]).convert("L")).astype(float)
    # suppress thin grid lines: minimum over a 3x3 then smooth slightly
    sm = gaussian_filter(im, 1.2)
    res = {}
    ov = Image.open(cfg["image"]).convert("RGB")
    dr = ImageDraw.Draw(ov)
    for name, c in cfg["curves"].items():
        t = track(sm, tuple(c["seed"]), c["xlo"], c["xhi"], c.get("win", 10))
        res[name] = t
        pts = list(t.items())
        dr.line(pts, fill=(255, 0, 0), width=3)
    ov.save(cfg["overlay"])
    json.dump({k: {str(x): y for x, y in v.items()} for k, v in res.items()}, open(cfg["out"], "w"))
