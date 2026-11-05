# Digitizing RM L51L28

Scripts as run on the 300 dpi page images (naca-archive's debug-figs workflow,
`19930086559:14,...,36`; the images are not kept here). Run from a folder holding
`19930086559/pages/p<N>-<N>.png`.

- **Figure 4 (blade form):** `track_fig4.py fig4_trace_config.json` traces h/b;
  `beta_pts.py` reads β column by column (runs off the grid lines, nearest to a hand-placed
  guide through the curve, outliers against a local quadratic dropped). b/D is a straight
  line at 0.0666 to x = 0.96.
- **Figures 6-8 (CT and CP against tip Mach number):** `run_mt2.py 16 17 18 19 20 21`.
  - `dig2.Grid` fits every grid line as a straight line, because the pages are skewed by up
    to 0.36° (12 px across a chart) and not uniformly (rows and columns lean differently).
    Pixels map to values by interpolating between the lines at the point.
  - `auto.py` finds curve pieces as connected components of the ink with the grid lines
    blanked, and follows each with `dig2.follow`, skipping grid crossings.
  - `merge.py` joins pieces of one curve (overlap or a consistent extrapolation across the
    gap). `run_mt.py` orders the curves by value at Mt 0.45, which is the blade-angle order.
  - `run_mt2.py` extends each end with `dig2.extend`, which may follow a curve along a grid
    line but drops a trailing stretch that is only grid line.
- `sample.py mt_full.json ..` writes `traces_fig6-8.json` and `static_fig6-8.csv`.
- Read by eye: figure 8(a), 12°, past Mt 0.89 (it runs along the CT = 0.105 grid line).
- Overlays of every trace were checked. Spot checks against the cross-plots of figures
  12-14 agree within about 0.002 in CT and 0.001 in CP.
