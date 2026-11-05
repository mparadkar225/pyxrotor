"""TM X-73612 figure 11: local axial Mach number near the blade leading edge at M0 = 0.8
(potential flow about the 0.35 D nacelle). Ticks: r/R 0.2-1.0 at x 502-1289; M 0.8, 0.7, 0.6
at y 658, 782, 910 (0.6 is the axis line)."""
import sys, json
import numpy as np
from PIL import Image
im = np.array(Image.open(sys.argv[1]).convert('L')).astype(int)
pts = []
Y = 870.0
for X in range(518, 1286, 4):
    col = np.where(im[int(Y) - 12:int(Y) + 13, X] < 128)[0] + int(Y) - 12
    if not len(col): continue
    g = np.split(col, np.where(np.diff(col) > 2)[0] + 1)
    Y = min([np.mean(a) for a in g], key=lambda v: abs(v - Y))
    pts.append((0.2 + (X - 502) / (1289 - 502) * 0.8, 0.8 - (Y - 658) / (910 - 658) * 0.2))
json.dump({"local_axial_mach_M0_0.8": pts}, open('fig11_trace.json', 'w'))
for r in (0.21, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
    print(r, round(float(np.interp(r, *np.array(pts).T)), 4))
