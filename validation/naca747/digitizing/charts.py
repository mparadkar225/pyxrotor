"""Chart boxes and calibrations for NACA 747 figures 8-13 (see ../README.md and tests.json notes:
figures 9c, 13a, 13b were later shifted by one grid step, 9c also in J; 13c uses fig13c in fig12.py)."""
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
# 300 dpi page renders: pdftoppm -r 300 -png tr747.pdf hi  ->  hi-08.png ... (pages 8-11 used)
D = os.environ.get("NACA747_PAGES", ".")
from digitize import Chart
SPEC = {   # page, box, x calib, y calib
 "fig8":   (8,  (470, 1790, 2430, 2655), (0.0, 0.1), (0.27, -0.01)),
 "fig9a":  (9,  (370, 570, 1150, 1350),  (0.0, 0.1), (0.0, 0.01, "last")),
 "fig9b":  (9,  (1215, 570, 2325, 1350), (0.0, 0.1), (0.0, 0.02, "last")),
 "fig9c":  (9,  (385, 1890, 2325, 2730), (0.0, 0.1), (0.0, 0.05, "last")),
 "fig11":  (10, (465, 1188, 2385, 2094), (0.0, 0.1), (0.0, 0.01, "last")),
 "fig12":  (10, (480, 2280, 2385, 3066), (0.0, 0.1), (0.0, 0.02, "last")),
 "fig13a": (11, (456, 711, 1176, 1605),  (0.0, 0.1), (0.0, 0.01, "last")),
 "fig13b": (11, (1254, 699, 2310, 1605), (0.0, 0.1), (0.0, 0.02, "last")),
 "fig13c": (11, (417, 1866, 2334, 2520), (0.0, 0.1), (0.0, 0.05, "last")),
}
def chart(name):
    page, box, xc, yc = SPEC[name]
    c = Chart.from_image(f"{D}/hi-{page:02d}.png", box=box, grid=230, grid_frac=0.3)
    c.calibrate(x=xc, y=[-v if i == 1 else v for i, v in enumerate(yc)] if len(yc) == 3 else yc)
    return c


def fig13c_50():
    """Dual 50 deg total CP. The scan is stretched by up to 4 px here, so values come from a
    piecewise-linear interpolation between the detected grid lines (one spurious line at
    y = 1880 dropped)."""
    import numpy as np
    from digitize import descending_branch, snap
    c = Chart.from_image(f"{D}/hi-11.png", box=(417, 1866, 2334, 2520), ink=160, grid=230, grid_frac=0.3)
    keep = [i for i, h in enumerate(c.hlines) if abs(h - 1880) > 4]
    c.calibrate(x=(0.0, 0.1), y=(1.0, -0.05), ylines=keep)
    sx, sy = snap(c, 2.55, 0.325, tol=0.2)
    x, y = c.trace((sx, sy), (1.8, 3.2))
    x, y = descending_branch(x, y, sx, y_min=0.06)
    h = np.array([v for i, v in enumerate(c.hlines) if i in keep])
    _, yp = c.to_px(x, y)
    return x, 1.0 - 0.05 * np.interp(yp, h, np.round((h - h[0]) / 32.1))
