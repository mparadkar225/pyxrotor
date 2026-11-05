"""Data-point symbols of the YP-47M efficiency cross-plots (efficiency against CP at constant
J), found with naca-archive's symbol detector (nacaarchive.symbols), then filtered: only
symbols inside the plot and outside the legend/label boxes, of the legend's shapes."""
import sys, json
import numpy as np
from PIL import Image, ImageDraw
sys.path.insert(0, '../../../../naca-archive')   # a naca-archive checkout beside pyxrotor
from nacaarchive.figures import find_grids
from nacaarchive import autodigitize as AD, symbols as SY


def detect(page):
    img = np.array(Image.open(page).convert('L'))
    g = find_grids(img)[0]
    raw, level = AD.ink_mask(img, g)
    x0, y0, x1, y1 = g.box
    with_grid = img[y0:y1 + 1, x0:x1 + 1] < level
    syms = SY.find_symbols(raw, (), raw=with_grid, vlines=[v - x0 for v in g.vlines],
                           hlines=[h - y0 for h in g.hlines])
    return img, g, [(s.kind, s.cx + x0, s.cy + y0, s.size) for s in syms]


def overlay(img, g, syms, out, tomap=None):
    o = Image.fromarray(img).convert('RGB'); d = ImageDraw.Draw(o)
    for i, (k, x, y, sz) in enumerate(syms):
        d.rectangle([x - sz / 2 - 3, y - sz / 2 - 3, x + sz / 2 + 3, y + sz / 2 + 3], outline=(255, 0, 0), width=3)
        d.text((x + sz / 2 + 4, y - 10), str(i), fill=(255, 0, 0))
    x0, y0, x1, y1 = g.box
    o.crop((x0 - 20, y0 - 20, x1 + 20, y1 + 20)).resize(((x1 - x0 + 40) // 2, (y1 - y0 + 40) // 2)).save(out)


PAGES = ["20030063991/pages/p23-23.png", "20030063991/pages/p33-33.png",
         "20030063996/pages/p18-18.png", "20030063996/pages/p27-27.png",
         "20030064146/pages/p13-13.png", "20030064146/pages/p22-22.png",
         "20030063262/pages/p21-21.png", "20030063262/pages/p31-31.png",
         "20030063262/pages/p41-41.png"]

if __name__ == "__main__":
    # syms_raw.json: every detected symbol (kind, x, y, size in page pixels) with the
    # plot's grid lines; symfull_*.png: numbered overlays used to assign each symbol a J
    res = {}
    for page in PAGES:
        img, g, syms = detect(page)
        key = page[:-4]
        res[key] = {"box": [int(v) for v in g.box], "vlines": [float(v) for v in g.vlines],
                    "hlines": [float(v) for v in g.hlines], "syms": syms}
        overlay(img, g, syms, "symfull_%s_%s.png" % (page[4:8], page.split("/")[-1][:-4]))
    json.dump(res, open("syms_raw.json", "w"))
