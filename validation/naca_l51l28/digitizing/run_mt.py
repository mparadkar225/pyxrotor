import sys, json; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from merge import *
BETAS = [0, 4, 8, 12, 16]
def digitize(p, minpts=60):
    im, G, cs = merged(p)
    cv = []
    for c in cs:
        xs = sorted(c)
        vals = np.array([G.value(x, c[x]) for x in xs])
        if len(xs) < minpts: continue
        if vals[0, 0] > 0.9 and vals[:, 1].max() < 0.009: continue   # NACA logo
        cv.append(vals)
    def at45(v):
        k = np.argsort(np.abs(v[:, 0] - 0.45))[:20]
        pp = np.polyfit(v[k, 0], v[k, 1], 1); return np.polyval(pp, 0.45)
    cv.sort(key=at45)
    return im, G, cv
if __name__ == "__main__":
    allv = {}
    for p in map(int, sys.argv[1:]):
        im, G, cv = digitize(p)
        o = Image.fromarray(im.astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(o)
        cols = [(255,0,0),(0,160,0),(0,0,255),(200,0,200),(0,170,170),(230,120,0),(0,0,0)]
        for i, v in enumerate(cv):
            for mt, y in v[::2]:
                x, yy = G.pixel(mt, y); d.ellipse([x-3, yy-3, x+3, yy+3], fill=cols[i % 7])
            x, yy = G.pixel(*v[-1]); d.text((x + 10, yy - 20), f'{BETAS[i] if i < 5 else "?"}', fill=cols[i % 7])
        print(p, len(cv), [(round(v[0,0],2), round(v[-1,0],2)) for v in cv])
        o.resize((o.width // 2, o.height // 2)).crop((150, 150, 1150, 1450)).save(f'mt{p}.png')
        allv[p] = [v.tolist() for v in cv]
    json.dump(allv, open('mt_' + '_'.join(sys.argv[1:]) + '.json', 'w'))
