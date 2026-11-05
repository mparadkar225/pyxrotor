"""RM L6L27a figure 28 (efficiency difference and helical tip Mach number against airspeed):
trace both curves column by column from seeds at 225 mph, skipping grid lines.
Run from a folder holding 20050019371/pages/p69-69.png (naca-archive debug-figs, 300 dpi).
Grid (found from long dark runs): 200 mph at x 669, 450 mph at 1865.5; 0 at y 2184, 14 at 501.
"""
import numpy as np
from PIL import Image

im = np.array(Image.open("20050019371/pages/p69-69.png").convert("L")).astype(int)
X0, SX = 669, (1865.5 - 669) / 250
Y0, SY = 2184, (2184 - 501) / 14
gx = np.array([669.0, 788.0, 908.5, 1027.5, 1148.5, 1267.5, 1388.5, 1507.5, 1627.0, 1747.5, 1865.5])
gy = np.array([Y0 - k * SY for k in range(15)])


def clusters(x):
    col = im[480:2190, x - 1:x + 2].min(1)
    ys = [i + 480 for i in range(len(col)) if col[i] < 140 and np.min(np.abs(gy - (i + 480))) > 4]
    cl = []
    for y in ys:
        if cl and y - cl[-1][-1] <= 2:
            cl[-1].append(y)
        else:
            cl.append([y])
    return [np.mean(c) for c in cl]


res = {}
for name, seed in (("diff", 2.8), ("mt", 2.1)):
    prev, out = Y0 - seed * SY, {}
    for mph in range(226, 451):
        x = int(round(X0 + (mph - 200) * SX))
        if np.min(np.abs(gx - x)) < 5:
            continue
        c = clusters(x)
        if c:
            y = min(c, key=lambda v: abs(v - prev))
            if abs(y - prev) < 45:
                prev, out[mph] = y, (Y0 - y) / SY
    k = sorted(out)
    res[name] = (np.array(k), np.array([out[m] for m in k]))
print("airspeed_mph,difference_percent,helical_tip_mach")
for m in range(225, 451, 25):
    print(f"{m},{np.interp(m, *res['diff']):.2f},{0.4 + 0.1 * np.interp(m, *res['mt']):.3f}")
