# Digitizing the propfan reports

The scripts are as they were run on the 300 dpi page images from naca-archive's debug-figs
workflow. The images themselves aren't kept here. Each script was run from a folder holding
`<NTRS id>/pages/p<N>-<NNN>.png`.

**CR-3505 figure 16** (SR-3 b/D, t/b, CLD, Δβ): `fig16.py 19820011304/pages/p58-058.png`
- The page is rotated upright. The grid is removed with a 4×4 binary opening, which keeps
  the 4–6 px curves and drops the 1–2 px grid lines.
- Calibration comes from the solid grid lines, with tilt and shear fitted.
- Δβ and t/b cross at 0.905R, so their outer parts are traced from seeds of their own.
- `steep_root()` traces the near-vertical t/b inboard of 0.305R row by row.

**CR-3505 figure 17** (aerodynamic sweep): `fig17.py`
- The grid is removed with the same opening.
- The curve is taken as the connected component through a seed, and points off a local
  quadratic are dropped.
- Calibration is bilinear between the frame corners.

**TM X-73612 figure 13(c)** (SR-2): `fig13c.py`
- The figure is gridless. Calibration comes from the t/b and r/R ticks.
- Steep parts are traced row by row.
- The b/D drop at the tip was read by eye.

**TM X-73612 figure 11** (nacelle local Mach number): `fig11.py`

**TM-79124 figures 12–13:** `run_panels.py 19790011898/pages`
- Each panel is calibrated from its own axis ticks, to within 0.025 points.
- The four curves are followed from seeds where they are apart.
- The traces stop where the curves cross.

Every trace was checked on an overlay.
