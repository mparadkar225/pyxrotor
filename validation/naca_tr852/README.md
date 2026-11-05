# NACA TR-852: drag of propeller-shank sections in flight

Barlow 1946, "Flight investigation at high speeds of the drag of three airfoils and a
circular cylinder representing full-scale propeller shanks" ([NTRS 19930091924](https://ntrs.nasa.gov/citations/19930091924)).

The models were 4.5 in thick (average shank size) and spanned 20.75 in. They were mounted
under an XP-51 wing at local Mach numbers 0.30–0.83 and Reynolds numbers 5e5–1.1e6 on the
diameter:
- a circular cylinder;
- NACA 16-025;
- NACA 16-040;
- NACA 16-040 with its rear quarter cut off.

**Cylinder** (figure 11, force tests). This is `pyxrotor.shank.TR852_CYLINDER`, the curve
the shank model uses. Section drag on the diameter:
- 0.58 at Mach 0.30, 0.51 at 0.45–0.48 (supercritical at these Reynolds numbers);
- rising to 0.95 at 0.60 and 1.09 at 0.70;
- 0.86 at 0.83.

The wake surveys gave up to 50% more; the report puts that on the oscillating wake. The force
and pressure-distribution tests agree.

**Airfoil fairings** (figure 16, `shank_sections_fig16.csv`). Profile drag on the thickness,
at zero angle of attack, traced on the steep parts; the low-speed levels are from the text.
`tests/measurement/test_shank_drag.py` compares them with `pyxrotor.sections.naca16`:
- **16-025.** Above its low-speed level, the drag rise follows the calibrated Mcrit
  (0.85 − 1.2 t + 0.15) within 0.015 (on the chord) from Mach 0.63 to 0.82.
- **16-040.** The same rule without the floor of 0.5 (0.85 − 1.2 t + 0.15 = 0.52) follows the
  rise within 0.03 above Mach 0.63.

  `pyxrotor.sections` keeps the floor (Mcrit ≥ 0.5 + 0.15), which starts the rise 0.1 in
  Mach too late: by Mach 0.74–0.79 the drag is 0.15–0.25 too low. Removing the floor costs
  0.06–0.11 in the P-47 comparisons, whose cuffed roots aren't modeled, so it stays.
- **Low-speed levels.** The estimate's drag is too high for these thick fairings: 0.047 on
  the chord for the 16-040, against 0.016 measured.

Digitizing: the page images (naca-archive debug-figs) were calibrated from their grid lines.
The 16-series curves were traced row by row on the steep parts, and the cylinder's faired
force-test curve was read by eye (±0.03).
