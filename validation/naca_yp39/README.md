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
- **Shanks.** By default the round blade shanks are not modeled: the section estimates there
  are Clark Y at the shank thickness. The test shows the shank drag growing with speed:
  entire minus working-part efficiency goes from −0.02 at 273 mph to −0.24 at 468 mph. So the
  main comparison is with the working-part efficiency; the shank model is compared with the
  entire propeller below.

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

### The round shanks (`pyxrotor.shank`)
With `run(row, shank=(r_round, r_airfoil))` the shanks are round to r_round and faired
linearly into the Clark Y blade by r_airfoil. The round part carries TR-852's measured
cylinder drag (`validation/naca_tr852`). The model is then compared with the *entire*
propeller's efficiency. Mean (model − test), Mcrit + 0.15:

| shank | 0.94–0.99 | 0.99–1.035 | 1.035–1.07 |
|---|---|---|---|
| none | +0.111 on average | | |
| round to 0.25R, faired by 0.30R | +0.007 | −0.006 | −0.029 |
| round to 0.27R, faired by 0.35R | −0.022 | −0.079 | −0.145 |
| round to 0.30R, faired by 0.40R | −0.068 | −0.185 | −0.300 |

- **Shank extent.** The report doesn't give it. Report 747's h/b is off its chart (> 0.44)
  inboard of 0.30R, and the survey's total head is negative to about 0.33–0.4R. Round to
  0.25R and faired by 0.30R fits.
- **Growth with speed.** With that extent, the model's shank loss (entire minus blades) goes
  from 0.024 at 273 mph to 0.203 at 468 mph. The test's entire minus working-part efficiency
  goes from 0.024 to 0.243, and each run agrees within 0.045.
- **Sensitivity.** The result depends strongly on the shank's extent: moving the round end out by 0.05R
  and the fairing by 0.1R more than doubles the loss at the highest speeds.

`tests/measurement/test_yp39.py` and `tests/measurement/test_shank_drag.py` pin these
findings.

## Limits
- One flight test, with efficiencies read from a small figure.
- The shank's extent is fitted, not measured; its drag is TR-852's 2-D cylinder.
