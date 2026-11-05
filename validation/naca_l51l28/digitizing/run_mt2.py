import sys, json; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from merge import *
from run_mt import digitize, BETAS
def full(p):
    im, G, cv = digitize(p)
    x0, x1 = int(G.xs_at(1500)[0]) + 6, int(G.xs_at(1500)[-1]) - 6
    out = []
    for v in cv:
        pts = {}
        for mt, y in v:
            x, yy = G.pixel(mt, y); pts[int(round(x))] = yy
        for d in (+1, -1):
            pts.update(extend(im, G, pts, d, x0, x1))
        xs = sorted(pts)
        out.append(np.array([G.value(x, pts[x]) for x in xs]))
    return im, G, out
if __name__ == "__main__":
    allv = {}
    for p in map(int, sys.argv[1:]):
        im, G, cv = full(p)
        o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
        cols = [(255,0,0),(0,160,0),(0,0,255),(200,0,200),(0,170,170)]
        for i, v in enumerate(cv):
            for mt, y in v[::2]:
                x, yy = G.pixel(mt, y); d.ellipse([x-3, yy-3, x+3, yy+3], fill=cols[i % 5])
        print(p, len(cv), [(round(v[0,0],3), round(v[-1,0],3)) for v in cv])
        o.resize((o.width // 2, o.height // 2)).crop((150, 150, 1150, 1450)).save(f'mt{p}.png')
        allv[p] = [v.tolist() for v in cv]
    json.dump(allv, open('mt_full.json', 'w'))
