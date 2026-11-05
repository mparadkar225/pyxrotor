# XROTOR against the YP-39 high-speed propeller flight test

A. W. Vogeley, "Preliminary results of flight tests of a conventional three-blade propeller
at high speeds", NACA, 1942 ([NTRS 20150011182](https://ntrs.nasa.gov/citations/20150011182)).

- **Propeller.** A Hamilton Standard 3155-6 three-blade propeller on a Bell YP-39, 10 ft 3/4 in in
  diameter, with a spinner to 12 3/8 in radius (0.205R) and no cuffs. The blade is the
  Clark Y blade of NACA Report 747, so `model.py` reuses that comparison's traced geometry
  and section estimate (`validation/naca747`).
- **Runs.** The tests were level-flight runs at 273–468 mph and 15,800–21,800 ft density altitude,
  with the blade angle locked at the governed setting. That gives a flight Mach number of 0.39–0.66 and a
  helical tip Mach number of 0.95–1.06.
- **Efficiency.** It was measured with a rake behind the propeller. The report gives two values:
  - "entire propeller";
  - "working part": the thrust of the positive-thrust part of the disk over the full torque, so
    slightly low.

## Data (`flight_runs.csv`)
- Speed, Mach number and tip Mach number come from the figure 4 labels. Density altitude and
  brake horsepower come from Table I, checked on the page image.
- Both efficiencies were read by eye from figure 5, to about ±0.005.
- Engine speed 3000 rpm through a 1.8:1 gear gives 1666.7 propeller rpm. The reported tip Mach
  numbers agree, at 1656–1671 rpm.

## Model (`model.py`)
- **Air.** The speed of sound comes from the reported V and M, and the density from the density
  altitude (ISA).
- **Blade angle.** `pyxrotor.robust.solve_power` at fixed rpm sets it: XROTOR's power solve
  (POWE), continued in Mach number from the solution without compressibility if it fails,
  then a fixed-pitch scan. Every run has a solution with every setting.
- **Shanks.** The round blade shanks (spinner to about 0.35R) are not modeled. The section
  estimates there are Clark Y at the shank thickness. The test shows the shank drag growing
  with speed: entire minus working-part efficiency goes from −0.02 at 273 mph to −0.24 at
  468 mph. So the comparison is with the working-part efficiency.

## Results (`compare.py`, `results.json`)
Mean of (model − test working-part efficiency) by tip Mach number band. The brackets give
converged runs out of runs:

| setting | 0.94–0.99 | 0.99–1.035 | 1.035–1.07 |
|---|---|---|---|
| 7.55 (Mcrit = 0.85 − 1.2 t) | −0.115 [3/3] | −0.299 [3/3] | −0.885 [4/4] |
| Mcrit + 0.10 | +0.004 | −0.044 | −0.133 |
| **Mcrit + 0.15** | +0.026 [3/3] | +0.001 [3/3] | −0.029 [4/4] |
| Mcrit + 0.20 | +0.035 | +0.018 | +0.014 |
| no compressibility | +0.040 | +0.033 | +0.042 |

- **Mcrit + 0.15 and + 0.20** both follow the test to tip Mach 1.06, within 0.035. The
  run-to-run scatter of the test (±0.04) can't separate them. The P-47 propellers, also
  Clark Y-type sections, rule out + 0.20 (`validation/naca_p47`); `pyxrotor.sections`'
  calibrated rule gives + 0.15 here.
- **7.55** loses 0.1–0.9 too much.
- **Without compressibility** the model is 0.03–0.04 high, so the loss the test shows is real.
  It is small, however, compared with what 7.55 predicts.

`tests/measurement/test_yp39.py` pins these findings.

## Limits
- One flight test, with efficiencies read from a small figure.
- The shank drag is excluded by comparing with working-part efficiency, which is itself
  slightly low.
