# XROTOR's compressibility model against P-47 flight tests

Three propellers were flight-tested on Republic P-47s at Langley, each with tables of every
run:

| key | report | propeller | blades, diameter | runs | tip Mach |
|---|---|---|---|---|---|
| `l177` | NACA ACR L4L07 / WR L-177, A. W. Vogeley, Dec 1944 ([NTRS 19930093302](https://ntrs.nasa.gov/citations/19930093302)) | Curtiss No. 714-1C2-12, Clark Y | 4, 12 ft 2 in | 62 | 0.78–1.10 |
| `l594` | NACA MR L5G09a / WR L-594, J. J. Gardner, 1945 ([NTRS 19930093653](https://ntrs.nasa.gov/citations/19930093653)) | Hamilton Standard No. 6507A-2, design CL 0.4 | 3, 12 ft 10 in | 112 | 0.81–1.17 |
| `l7d10` | NACA RM L7D10, J. J. Gardner, Apr 1947 ([NTRS 20050019331](https://ntrs.nasa.gov/citations/20050019331)) | Curtiss No. 838-1C2-18 with trailing-edge extensions, design CL 0.5 | 3, 12 ft 10 in | 112 | 0.83–1.14 |

All three have 2:1 gearing on a P&W R-2800. In each, thrust comes from a slipstream
total-pressure survey and torque from the engine torque meter. RM L7D10 notes that its
survey values are not corrected for slipstream rotation, so the reported thrust and
efficiency are 1.5–4% high. The earlier reports use the same method.

## Data
- **`<key>_table*.csv`:** every flight-data table, transcribed by hand. The source images
  are the table pages, rendered at 300 dpi from naca-archive's cached PDFs by its debug-figs
  workflow. The archive's OCR of these pages was unusable.
  - Each row passes η = C_T·J/C_P and M_t = M·√(1 + (π/J)²) to the tables' rounding
    (`test_tables_are_consistent`).
  - One printed tip Mach fails the check: RM L7D10 run 22-15 prints 1.090 where M and J give
    1.077. It is kept as printed; the model uses M.
  - RM L7D10's run labels "25 to 1" are written 25-1. L-594's run 5-5 has M printed .514
    (bold), which is .314 by its J and M_t.
- **`geometry_<key>.json`:** each figure 1, the blade-form curves (b/D, h/b, θ − β), traced
  with `digitize_track.py` (configs `fig1_<key>_trace_config.json`).
  - Calibration comes from axis labels or the frame. As a check, θ − β reads −0.07°, −0.01°
    and −0.001° at 0.75R, where it is 0 by definition.
  - Points read by eye, where the tracer loses a line on a heavy grid line or where curves
    crowd, are listed in each file's notes. Cuffs are not modeled.
  - L-594 and L7D10 don't name their blade sections.

## Model (`model.py`)
Each run is solved at the test's J, C_P, rpm, density and Mach number. The blade angle is
found the way the governor set it, at fixed rpm, with `pyxrotor.robust.solve_power`:
1. XROTOR POWE;
2. if that fails, POWE continued in Mach number from the solution without compressibility;
3. if that fails too, interpolation between fixed-pitch solves, accepting small limit
   cycles.

Before that solver was added, a failed run was only restarted once from the blade angle
without compressibility. That left 5 (L-594) and 9 (L7D10) of 112 runs unconverged with
Mcrit + 0.15.

All three use the same Clark Y estimate from thickness as the R&M 2357 and NACA 747
comparisons. That includes Mcrit = 0.85 − 1.2 t at the outer blade: 0.76 for the Curtiss
714 (t ≈ 0.072), 0.77 for the 6507A-2 (0.065) and 0.79 for the 838 (0.051).

XROTOR's compressibility terms are the thing under test:
- the drag rise 10(M − Mcrit_eff)³, where Mcrit_eff falls by 0.25·|CL − CL at minimum drag|;
- the compressible CLmax limit;
- the Prandtl–Glauert lift slope;
- the local Mach number capped at 0.99.

## Results (`compare.py`, `results.json` for l177, `results_<key>.json`)
Mean (model − test) efficiency by tip Mach, as [converged runs/runs]:

| | M_t < 0.85 | 0.85–0.90 | 0.90–0.96 | 0.96–1.02 | 1.02–1.08 | 1.08–1.17 |
|---|---|---|---|---|---|---|
| **Curtiss 714 (L-177)** | | | | | | |
| 7.55 as configured | +0.007 [25/25] | −0.029 [14/14] | −0.086 [5/5] | −0.145 [5/5] | −0.472 [5/11] | — [0/2] |
| Mcrit + 0.10 | +0.031 | +0.024 | −0.006 | +0.006 | −0.251 | −0.466 |
| **Mcrit + 0.15** | +0.034 | +0.038 | +0.017 | +0.043 | +0.027 [11/11] | +0.005 [2/2] |
| no compressibility | +0.035 | +0.050 | +0.054 | +0.088 | +0.227 | +0.324 |
| **Hamilton Std 6507A-2 (L-594)** | | | | | | |
| 7.55 as configured | +0.052 [7/7] | −0.005 [19/19] | −0.072 [31/31] | −0.153 [20/20] | −0.386 [15/19] | −0.812 [11/16] |
| Mcrit + 0.10 | +0.064 | +0.041 | +0.015 | −0.007 | −0.044 [18/19] | −0.087 [16/16] |
| **Mcrit + 0.15** | +0.064 | +0.045 | +0.033 | +0.026 | +0.036 [18/19] | +0.067 [16/16] |
| no compressibility | +0.065 | +0.049 | +0.045 | +0.053 | +0.107 | +0.185 |
| **Curtiss 838 (L7D10)** | | | | | | |
| 7.55 as configured | +0.021 [6/6] | −0.004 [25/25] | −0.057 [28/28] | −0.118 [21/21] | −0.296 [16/16] | −0.985 [11/16] |
| Mcrit + 0.10 | +0.027 | +0.023 | +0.010 | +0.005 [21/21] | −0.004 [16/16] | −0.126 [15/16] |
| **Mcrit + 0.15** | +0.027 | +0.024 | +0.021 | +0.029 [21/21] | +0.052 [16/16] | +0.058 [16/16] |
| no compressibility | +0.028 | +0.026 | +0.028 | +0.047 | +0.089 | +0.164 |

(For L-594 and L7D10 the first column is tip Mach 0.80–0.85; no run is lower.)
`compare.py` prints the full sweep, including Mcrit + 0.05 and + 0.20.

- **Low-speed bias.** Without compressibility the model is high by about +0.035 on the
  Curtiss 714, +0.05–0.065 on the Hamilton Standard and +0.027 on the Curtiss 838. This is
  the section estimate, as on R&M 2357, where efficiency was also high.
  - The Hamilton Standard's sections are probably not Clark Y; the report gives only its
    design CL of 0.4.
  - If the reports' thrust is 1.5–4% high from uncorrected slipstream rotation, the model's
    bias against the true efficiency is larger by about that.
- **7.55 as configured overpredicts the losses on all three, and by about as much.** Losses
  begin below tip Mach 0.90, where the tests show almost none.
  - At tip Mach 0.96–1.02 the model is 0.12–0.15 too low.
  - Above 1.02 it is 0.3–1.0 too low on the converged runs, and a third of those runs don't
    converge.
  - At a local Mach number of 1.0 a section with Mcrit 0.76 is past its compressible-stall
    limit: CLmax falls to about 0.26, and the cubic drag rise adds over 0.2 in CD
    (10 × (1.0 − 0.70)³, with Mcrit_eff about 0.70).
- **Mcrit + 0.15 reproduces the losses on all three.** In every band the error stays within
  about 0.03 of the propeller's own low-speed bias, up to tip Mach 1.10 on the Curtiss 714 and 1.17
  on the others.
  - Mcrit + 0.10 already fits up to tip Mach 1.02. Above that it leaves losses too large
    (by up to 0.47).
  - On the two three-blade propellers the error creeps up by about 0.02–0.03 above tip Mach
    1.02, so their best single shift is nearer +0.13.
  - The shift is a lumped correction: it also absorbs 3-D tip relief, the drag law's shape,
    and the compressible-stall limit. It holds for 5–7% thick tips and helical tip Mach up
    to about 1.15. It is not a general section Mcrit.
- **Convergence.** With Mcrit + 0.15, 1 of the 112 L-594 runs has no usable solution
  (tip Mach 1.0–1.02), and every L7D10 run has one. XROTOR's own POWE fails on more of them:
  more Newton iterations don't help, because the solution oscillates near the local-Mach cap.
  The Mach continuation and the fixed-pitch scan of `pyxrotor.robust` recover them.

`tests/measurement/test_p47_compressibility.py` pins the findings for each propeller.

## Limits
- One section estimate for three blade families.
- Thrust comes from wake surveys, uncorrected for slipstream rotation.
- The airplane's spinner and cowl effects are part of the measurement and not of the model.
- Static tests of thin 16-series blades (`validation/naca_l51l28`) want about +0.20.
