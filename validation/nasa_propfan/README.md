# XROTOR against the SR-2 and SR-3 propfan models

Both propfans are eight-blade models with a 62.23 cm reference diameter. They were tested in
the NASA Lewis 8x6 ft wind tunnel at Mach 0.45–0.85, behind a nacelle 0.35 D across.

| | SR-2 | SR-3 |
|---|---|---|
| tip sweep | 0° (straight) | 45° (aerodynamic sweep: −11° at the root, 46° at 0.95R) |
| activity factor, integrated design CL | 203, 0.08 | 235, 0.214 |
| sections | NACA 65/circular-arc to 0.37R, NACA 16 from 0.53R; t/b 0.02 at the tip | the same t/b |

Sources:
- NASA TM-79124, Jeracki, Mikkelson & Blaha 1979 ([NTRS 19790011898](https://ntrs.nasa.gov/citations/19790011898)): figure 7 (geometry), figures 12–13 (efficiency against J at 70–100% of design power loading, Mach 0.6 and 0.8), figure 14 (efficiency against Mach number at the design point).
- NASA CR-3505, Rohrbach et al. 1982 ([NTRS 19820011304](https://ntrs.nasa.gov/citations/19820011304)): SR-3 geometry (figures 14, 16, 17), Tables I and III, and the fixed-pitch power coefficients on p. 19.
- NASA TM X-73612, Mikkelson et al. 1977 ([NTRS 19770013165](https://ntrs.nasa.gov/citations/19770013165)): SR-2 geometry (figure 13(c)) and the nacelle flow (figure 11).

XROTOR is a straight lifting line with 2-D sections, so it has no sweep. Two questions:
- How far does it get on propfans?
- What would a sweep correction have to do?

## Data
- **`geometry.json`.** b/D, t/b, twist, design CL and sweep against r/R, traced
  (`digitizing/`). The fine-grid CR-3505 charts were traced after a binary opening removed
  the grid.
  - The traces reproduce the published integrals: activity factor 202 and 234 (published 203
    and 235), integrated design CL 0.075 and 0.212 (0.08 and 0.214).
  - The twist reads −0.11° and +0.02° at 0.75R, where it is 0 by definition.
  - The SR-3 tip diameter grows with blade angle: 1.019 Dref at 60° (CR-3505 fig. 14, read by
    eye).
- **`design_point_eta.csv`.**
  - SR-3: CR-3505 Table I, checked against the page image. Row 1 is at JREF 3.06, CPREF 1.695
    on the reference diameter; row 2 is the same on the true diameter.
  - SR-2 and SR-3: TM-79124 figure 14, symbols read by eye to about ±0.15 points. Its SR-3
    points agree with Table I within 0.7 points.
- **`power_loading_eta.csv`.** TM-79124 figures 12(a, e) and 13(a, e), traced
  (`traces_tm79124_fig12_13.json`), with CP = 0.05933·f·J³. Where the curves cross at low J
  they were not traced.
- **`sr3_fixed_angle.csv`.**
  - CR-3505 p. 19: at JREF 3.1 and β 60.5°, CPREF is 1.64 at Mach 0.70, 1.49 at 0.80 and 1.23
    at 0.85, with CT down 11% and 29%.
  - CR-3505 Table III: β 60.7°, on the true diameter. On Dref it is the same curve (CP 1.65 at
    J 3.06, Mach 0.70).
  - The design CP and J at Mach 0.80 need β = 61.7°.

## Model (`model.py`)
- **Sections.** The 16-series estimate of the RM L51L28 comparison (`validation/naca_l51l28`)
  is used at every station, including the 65/circular-arc root. Mcrit = 0.85 − 1.2 t, which
  is 0.826 at the 2% tip, plus a shift.
- **Nacelle.** XROTOR's body velocity is taken from the potential-flow axial Mach number
  about the nacelle (TM X-73612 fig. 11; UBODY −0.21 at the root to −0.01 at the tip).
- **Sweep options:**
  - `sweep="mcrit"` divides each section's Mcrit by cos(sweep).
  - `sweep="cosine"` also multiplies the lift slope by cos(sweep).
  - `sweep="core"` uses the compiled core's swept-blade model (`pyxrotor.sweep`). The
    quarter-chord line is built from the aerodynamic sweep, its planform adjusted until the
    core's effective sweep at the design point is the published one (within 0.2°). Each
    station's sweep then follows the local flow at every iteration, and the sections are
    evaluated in the normal plane: normal Mach number and angle, CL·cos², pressure drag ·cos³.
    The thickness is taken as the normal section's, t/b / cos(sweep).
- **Air.** Sea-level total conditions (the open-circuit tunnel draws atmospheric air).
- **Finding the blade angle.** XROTOR's power solve (POWE) often fails near the local-Mach
  cap, so the blade angle for a given CP is found outside it:
  - fixed-pitch solves on a 0.5° grid;
  - linear interpolation at the target CP;
  - near the cap the Newton iteration can settle into a small limit cycle. The mean of 8
    iterations is accepted when efficiency varies by less than 0.005 and CP by less than 2%.
  - A point is unusable if the bracketing solves are more than 2° apart, or differ in
    efficiency by more than 0.06.
  - `pyxrotor.robust` now has this method for any rotor (with stricter limit-cycle checks).
    This model keeps its own copy because the SR-3's diameter changes with blade angle, so
    each blade angle needs a new rotor.

## Results (`compare.py`, `results.json`)
**Design point, J 3.06, CP 1.7: model − test net efficiency** ("—" means no usable solution):

| | 0.60 | 0.70 | 0.75 | 0.80 | 0.85 |
|---|---|---|---|---|---|
| **SR-2** test | 0.781 | 0.780 | 0.773 | 0.759 | 0.734 |
| 7.55 | −0.021 | −0.294 | −0.493 | −0.641 | −0.744 |
| Mcrit + 0.15 | +0.011 | −0.011 | −0.056 | −0.079 | −0.080 |
| **Mcrit + 0.20** | +0.011 | +0.008 | −0.006 | −0.013 | +0.005 |
| no compressibility | +0.010 | +0.011 | +0.018 | +0.033 | +0.058 |
| **SR-3** test (Table I) | 0.802 | 0.805 | 0.796 | 0.788 | 0.770 |
| 7.55 | −0.016 | −0.209 | −0.387 | −0.534 | −0.696 |
| Mcrit + 0.15 | −0.007 | −0.018 | −0.042 | −0.037 | — |
| Mcrit + 0.20 | −0.007 | −0.016 | — | −0.005 | — |
| cosine sweep, Mcrit + 0.10 | −0.008 | −0.010 | — | +0.006 | +0.016 |
| core sweep (`pyxrotor.sweep`) | −0.015 | −0.023 | −0.027 | −0.066 | −0.182 |
| **core sweep, Mcrit + 0.10** | −0.014 | −0.017 | −0.009 | −0.003 | +0.004 |
| no compressibility | −0.009 | −0.011 | −0.002 | +0.006 | +0.024 |

- **SR-2, Mcrit + 0.20 follows the measured loss.** It is within 0.025 of its low-speed bias
  (+0.011) at every Mach number.
  - Figure 13(a), at Mach 0.8 and 70–100% loading, agrees: a mean error of +0.005.
  - This is the same shift the thin 16-series propellers of RM L51L28 wanted. Mcrit + 0.15,
    the P-47 correction for 5–7% Clark Y-like tips, loses 0.06–0.09 too much here.
  - 7.55 as configured is far off, with 0.3–0.75 too much loss. This is the P-47 result
    again, on 2% tips at helical tip Mach 1.14.
- **SR-3 loses little, and no setting can be told apart from no compressibility.** The test
  loses 3.2 points from Mach 0.6 to 0.85. The model without compressibility is within 0.024 at
  every Mach number.
  - Straight-blade settings up to Mcrit + 0.15 overpredict the loss.
  - Many solves at Mach ≥ 0.75 don't converge. Usable points at Mach 0.8 in figure 13(e):
    16–36 of 45 with compressibility, against 45 of 45 without.
  - The cosine sweep rule with Mcrit + 0.10 converges more often than the straight-blade
    settings and fits within 0.016.
- **SR-3 needs 5.8–8.2° more blade angle than XROTOR gives for the design power.** The test
  needs 61.7° at Mach 0.8; the model needs 53.5° to 55.9°. Hamilton Standard's own
  lifting-line programs were 4–5° low too (CR-3505). A straight lifting line overloads a
  swept blade.
  - On SR-2 the model's blade angle is 55.7–58.9° against about 58° in the earlier UTRC test
    (TM X-73612).
- **At fixed pitch, compressibility moves the power the wrong way.** From Mach 0.70 to 0.80
  and 0.85, SR-3's measured CP falls by 9% and 25%, and CT by 11% and 29%. CR-3505
  reads this as lower blade lift than designed, "sweep is not as effective as expected".
  - Every compressible XROTOR setting raises CP by 15–40%. This comes from its
    Prandtl-Glauert lift slope, with the local Mach number capped at 0.995 (the tip is at
    relative Mach 1.01–1.24).
  - So with a constant-speed governor the error appears as blade angle, not power. At a fixed
    blade angle the model's power and thrust trends with Mach number are wrong in sign.
- **Low-speed bias.** At Mach 0.6 the model is +0.015 on SR-2 and −0.012 on SR-3 (figs.
  12–13 means).

- **The core sweep model (stages 1–2 of a sweep model) moves both SR-3 discrepancies the
  right way, but only part of the way.**
  - **Efficiency.** With Mcrit + 0.10 the efficiency follows the test within 0.02 of the
    Mach 0.6 bias to Mach 0.85. Every solve converges, including all 45 points of figure
    13(e) at Mach 0.8 (16–36 with the straight-blade settings), at −0.011 on average.
  - **Blade angle.** The deficit at fixed pitch drops from 6.1° to 4.5°: the normal Mach
    number removes most of the Prandtl-Glauert lift gain, and cos L the lift slope.
  - **Power at fixed pitch.** From Mach 0.70 to 0.85 CP still rises: +10% against +22% straight,
    while the test falls 25%.
  - **What is left.** A sweep-aware strip theory, as Hamilton Standard's was, gets this far.
    The rest needs the induction of the swept lifting line (a curved bound vortex, the wake
    leaving from swept and raked positions) and a transonic lift loss.

`tests/measurement/test_propfan.py` pins these findings.

## Limits
- One section estimate for 65/CA and 16-series sections, with cascade effects at the root
  ignored. The nacelle velocity is from a potential-flow plot for the M0 = 0.8 case, used at
  every Mach number.
- The SR-3 diameter variation was read by eye.
- XROTOR has no sweep. `sweep="core"` adds it to the sections, from the blade's geometry and
  the local flow; the induction of a swept lifting line, tip relief and transonic lift loss
  are still missing. The cosine options touch only the section data.
- A sweep model would have to reproduce both the blade-angle offset and the fixed-pitch
  unloading with Mach number. Sweep (scimitar) blades at high tip Mach are not predicted by
  this model, beyond efficiency at a given power.
