import sys; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from dig2 import *
PAGES = {16: (26,), 17: (26,), 18: (26,), 19: (26,), 20: (26,), 21: (15,)}
def grid_for(p):
    im = load(f'19930086559/pages/p{p}-{p}.png'); n = PAGES[p][0]
    return im, Grid(im, n, 19, (0.2, 0.05), (0.005 * (n - 1), -0.005))
if __name__ == "__main__":
    for p in PAGES:
        im, G = grid_for(p)
        for mt in map(float, sys.argv[1].split(',')):
            x, _ = G.pixel(mt, 0.06); x = int(round(x))
            if G.on_vline(x, 1500, 6): x += 9
            top, bot = int(G.ys_at(x)[0]), int(G.ys_at(x)[-1])
            g = runs(im[top:bot, x], 140)
            cs = [(top + np.mean(a), len(a)) for a in g if len(a) <= 9]
            cs = [(int(y), n, round(G.value(x, y)[1], 4)) for y, n in cs if not G.on_hline(x, y, 4)]
            print(p, mt, x, cs)
