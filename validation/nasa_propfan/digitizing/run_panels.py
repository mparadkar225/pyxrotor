import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import trace
from dig2 import load
P = {  # name: page, box, y tick values, x tick values, seed J, seed eta (70..100%), J range, prop, M
 '12a': ('p19-19', (780, 450, 1800, 1100), [82, 80, 78, 76], [2.6, 2.8, 3.0, 3.2, 3.4, 3.6], 3.05, [80.57, 80.02, 79.17, 78.11], (2.6, 3.6), 'SR-2', 0.6),
 '12e': ('p20-20', (780, 1330, 1800, 1830), [84, 82, 80, 78, 76], [2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8], 3.2, [82.51, 81.83, 80.86, 79.5], (2.2, 3.8), 'SR-3', 0.6),
 '13a': ('p21-21', (780, 560, 1800, 1100), [80, 78, 76, 74], [2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8], 3.1, [77.06, 76.74, 76.13, 75.38], (2.6, 3.8), 'SR-2', 0.8),
 '13e': ('p22-22', (560, 800, 1500, 1260), [82, 80, 78, 76], [2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8], 3.3, [80.1, 79.65, 78.76, 77.54], (2.6, 3.8), 'SR-3', 0.8),
}
out = {}
for k, (pg, box, yv, xv, sj, se, xl, prop, M) in P.items():
    im = load(f'{sys.argv[1]}/p{pg.split("-")[0][1:]}-{pg.split("-")[1]}.png')
    tr, cal = trace(im, box, yv, xv, sj, se, xl, name=k)
    out[k] = {"prop": prop, "M": M, "curves": {str(a): b for a, b in tr.items()}, "calibration_max_error": cal[:2]}
    print(k, {a: (round(b[0][0], 2), round(b[-1][0], 2), len(b)) for a, b in tr.items()}, cal[:2])
json.dump(out, open('fig12_13_traces.json', 'w'))
