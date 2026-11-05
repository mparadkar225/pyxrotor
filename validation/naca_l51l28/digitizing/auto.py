import sys, json; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from dig2 import *
from scipy.ndimage import label, find_objects
from cand2 import PAGES, grid_for

def gridmask(im, G, pad=4):
    h, w = im.shape; yy, xx = np.mgrid[0:h, 0:w]
    m = np.zeros(im.shape, bool)
    for a, b in G.V: m |= np.abs(xx - (a * yy + b)) <= pad
    for a, b in G.H: m |= np.abs(yy - (a * xx + b)) <= pad
    return m

def curves(p, minw=120):
    im, G = grid_for(p)
    x0, x1 = int(G.xs_at(1500)[0]) + 6, int(G.xs_at(1500)[-1]) - 6
    y0, y1 = int(G.ys_at(1300)[0]) + 6, int(G.ys_at(1300)[-1]) - 6
    ink = (im < 140) & ~gridmask(im, G)
    ink[:, :x0] = ink[:, x1:] = False; ink[:y0] = ink[y1:] = False
    lab, n = label(ink, np.ones((3, 3)))
    segs = []
    for i, sl in enumerate(find_objects(lab)):
        w = sl[1].stop - sl[1].start; h = sl[0].stop - sl[0].start
        if w >= 50 and h < 0.6 * w + 40:
            ys, xs = np.where(lab[sl] == i + 1)
            xm = int(np.median(xs)) + sl[1].start
            ym = float(np.mean(ys[xs + sl[1].start == xm])) + sl[0].start
            segs.append((w, xm, ym))
    segs.sort(reverse=True)
    out = []
    for w, xm, ym in segs:
        if any(abs(np.interp(xm, list(map(int, c)), list(c.values())) - ym) < 6 and min(map(int, c)) <= xm <= max(map(int, c)) for c in out):
            continue
        c = follow(im, G, (xm, ym), x0, x1)
        c = {int(k): v for k, v in c.items()}
        if max(c) - min(c) >= minw: out.append(c)
    return im, G, out

if __name__ == "__main__":
    res = {}
    for p in map(int, sys.argv[1:]):
        im, G, cs = curves(p)
        o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
        cols = [(255,0,0),(0,160,0),(0,0,255),(200,0,200),(0,170,170),(230,120,0),(120,60,0),(0,0,0)]
        info = []
        for i, c in enumerate(cs):
            xs = sorted(c); col = cols[i % len(cols)]
            for x in xs: d.ellipse([x-3, c[x]-3, x+3, c[x]+3], outline=col)
            v0 = G.value(xs[0], c[xs[0]]); v1 = G.value(xs[-1], c[xs[-1]])
            d.text((xs[-1] + 8, c[xs[-1]] - 10), str(i), fill=col)
            info.append((i, round(v0[0], 3), round(v0[1], 4), round(v1[0], 3), round(v1[1], 4), len(xs)))
        print(p, info)
        o.resize((o.width // 2, o.height // 2)).crop((150, 150, 1150, 1400)).save(f'auto{p}.png')
        res[p] = [{str(k): v for k, v in c.items()} for c in cs]
    json.dump(res, open('auto_' + '_'.join(sys.argv[1:]) + '.json', 'w'))
