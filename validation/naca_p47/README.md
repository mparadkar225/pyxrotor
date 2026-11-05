# XROTOR's compressibility model against P-47 flight tests

**Source:** NACA ACR L4L07 / Wartime Report L-177, A. W. Vogeley, *Climb and High-Speed
Tests of a Curtiss No. 714-1C2-12 Four-Blade Propeller on the Republic P-47C Airplane*,
Langley, December 1944 ([NTRS 19930093302](https://ntrs.nasa.gov/citations/19930093302)).

The propeller has four blades with Clark Y sections, a diameter of 12 ft 2 in and 2:1
gearing, driven by a P&W R-2800-21. Efficiency was measured in flight: thrust from a
slipstream total-pressure survey, torque from the engine torque meter.

## Data
- **`l177_table1.csv`, `l177_table2.csv`:** tables I (42 climb runs) and II (20 high-speed
  runs), transcribed by hand. The source images are pages 14 and 15, rendered at 300 dpi
  from naca-archive's cached PDF by its debug-figs workflow. The archive's OCR of these
  pages was unusable.
  - Each row passes η = C_T·J/C_P and M_t = M·√(1 + (π/J)²) to the tables' rounding
    (`test_tables_are_consistent`).
  - The runs cover tip Mach 0.78–1.10 and airplane Mach up to 0.71.
- **`geometry_l177.json`:** figure 1, the blade-form curves (b/D, h/b, θ − β), traced with
  `digitize_track.py`.
  - Calibration comes from the axis labels. As a check, θ − β reads −0.07° at 0.75R, where
    it is 0 by definition.
  - Root values (r/R ≤ 0.30) were read by eye. The cooling cuffs are not modeled.

## Model (`model.py`)
Each run is solved at the test's J, C_P, rpm, density and Mach number, with XROTOR POWE at
fixed rpm. The blade angle is found the same way the governor set it.

The sections use the same Clark Y estimate from thickness as the R&M 2357 and NACA 747
comparisons. That includes Mcrit = 0.85 − 1.2 t, which is 0.76 at the outer blade
(t = 0.072).

XROTOR's compressibility terms are the thing under test:
- the drag rise 10(M − Mcrit_eff)³, where Mcrit_eff falls by 0.25·|CL − CL at minimum drag|;
- the compressible CLmax limit;
- the Prandtl–Glauert lift slope;
- the local Mach number capped at 0.99.

## Results (`compare.py`, `results.json`)
Mean (model − test) efficiency by tip Mach, as [converged runs/runs]:

| setting | M_t < 0.80 | 0.80–0.85 | 0.85–0.90 | 0.90–0.96 | 0.96–1.02 | 1.02–1.10 |
|---|---|---|---|---|---|---|
| 7.55, Mcrit = 0.85 − 1.2t | +0.015 [9/9] | +0.003 [16/16] | −0.029 [14/14] | −0.086 [5/5] | −0.145 [5/5] | −0.472 [5/13] |
| Mcrit + 0.05 | +0.024 | +0.023 | +0.002 | −0.038 | −0.048 [4/5] | −0.909 [9/13] |
| Mcrit + 0.10 | +0.027 | +0.034 | +0.024 | −0.006 | +0.006 | −0.284 |
| **Mcrit + 0.15** | +0.027 | +0.038 | +0.038 | +0.017 | +0.043 | +0.023 [13/13] |
| Mcrit + 0.20 | +0.027 | +0.040 | +0.045 | +0.034 | +0.065 | +0.155 [10/13] |
| no compressibility | +0.027 | +0.040 | +0.050 | +0.054 | +0.088 | +0.242 |

- **Low-speed bias.** Without compressibility the model is 0.03–0.05 high at moderate tip
  Mach. That is the section estimate, as on R&M 2357, where efficiency was also high. Thrust
  coefficient is 0.006 high, and the blade angle is 0.6° low on average (range −2.0 to +0.4°).
  - The test's own compressibility loss is the growth of that difference above it:
    about 0.01 at M_t 0.85–0.90, 0.05 at 0.96–1.02 and 0.2 at 1.02–1.10.
  - That matches the report, which saw losses begin at tip Mach 0.83–0.88 and fall about
    7% per 0.1 of tip Mach.
- **7.55 as configured overpredicts the losses.** They start below tip Mach 0.85, a band the
  test shows nearly free of loss. By tip Mach 0.92–1.0 the model is 0.09–0.15 too low.
  - Above 1.02, efficiency comes out 0.3–0.4 too low on the converged runs (0.33 against a
    measured 0.71 at M_t 1.08).
  - Most of the remaining runs at that tip Mach don't converge; their efficiencies come out
    negative.
  - At a local Mach number of 1.0 a section with Mcrit 0.76 is past its compressible-stall
    limit: CLmax falls to about 0.26, and the cubic drag rise adds over 0.2 in CD
    (10 × (1.0 − 0.70)³, with Mcrit_eff about 0.70).
- **Mcrit + 0.15 reproduces the losses.** Raising every section's Mcrit by 0.15 (to 0.91 at
  the outer blade) puts the error back to the low-speed bias in every band, +0.02 to +0.04,
  up to tip Mach 1.10. Every run converges, including the low-power runs at efficiency
  0.4–0.6.
  - The shift is a lumped correction: it also absorbs 3-D tip relief, the drag law's shape,
    and the compressible-stall limit.
  - It holds for this propeller, about 7% thick near the tip, and for helical tip Mach up to
    1.1. It is not a general section Mcrit.

`tests/measurement/test_p47_compressibility.py` pins all three findings.

## Limits
- One propeller and one section family. Thrust comes from a wake survey.
- The airplane's spinner and cowl effects are part of the measurement and not of the model.
- The next checks are the other P-47 reports in naca-archive:
  - WR L-594: a Hamilton Standard 6507A-2 three-blade propeller, tip Mach up to 1.16.
  - RM L7D10: a Curtiss 838 three-blade propeller with trailing-edge extensions.

  Their table pages are already rendered on the archive's `debug-figs` branch.
