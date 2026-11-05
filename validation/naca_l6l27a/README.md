# Shank drag: NACA RM L6L27a

Maynard 1947, "Aerodynamic characteristics at high speeds of full-scale propellers having
different shank designs" ([NTRS 20050019371](https://ntrs.nasa.gov/citations/20050019371)).

The report tests two 10 ft two-blade NACA propellers that differ only in their shanks:

| | NACA 10-(5)(08)-03 | NACA 10-(5)(08)-03R |
|---|---|---|
| shank | wide, thin 16-series sections to the spinner (h/b 0.15 at 0.217R) | conventional: an oval, h/b 0.73 at 0.2R, 0.3 at 0.3R |
| design CL inboard | falls to 0.16 at the spinner | 0.5 to the shank |

Outboard of about 0.45R the two blades are the same. Both use NACA 16-series sections, design
CL 0.5, 8% thick at 0.7R. They were run in the Langley 16-ft high-speed tunnel on a
2000 hp dynamometer, behind a spinner of 0.217 D.

## Data
- **`geometry.json`.** Figure 4 (b/D, h/b, blade angle, design CL), read by eye against its
  grid, with figure 7 for the 03R's shank sections. The 03R's h/b leaves the chart inboard of
  about 0.26R.
- **`fig28_difference.csv`.** Figure 28, traced: the efficiency difference
  (η₀₃ − η₀₃R)/η₀₃ at 1140 rpm and CP 0.2 (0.1 per blade), 225–450 mph, with its helical tip
  Mach numbers (0.61–0.79). These give the speed of sound, 341–344 m/s.
- **`fig27_efficiency.csv`.** Figure 27, read by eye to about ±0.01: efficiency against J at
  CP 0.10, 0.15 and 0.20, 1140 rpm, for both propellers.

## Model (`model.py`)
- **Sections.** `pyxrotor.sections.naca16` at each station's t/c and design CL, with the
  calibrated Mcrit.
- **The 03R's shank.** `pyxrotor.shank`: TR-852's cylinder drag, taken on the oval's
  thickness (h/b from figure 7), to 0.26R, faired linearly into the 16-series blade by 0.30R.
- **Operating point.** Constant power at 1140 rpm (`RoundShank.solve_power`, or
  `pyxrotor.robust.solve_power` for the 03).

## Results (`compare.py`, `results.json`)
Efficiency difference between the propellers, percent, at CP 0.2 and 1140 rpm:

| mph | tip Mach | test | model | model, no shank model |
|---|---|---|---|---|
| 225 | 0.61 | 2.9 | 2.8 | 0.6 |
| 275 | 0.64 | 3.8 | 3.4 | 0.9 |
| 325 | 0.68 | 5.1 | 4.6 | 1.4 |
| 375 | 0.72 | 7.1 | 7.4 | 2.2 |
| 400 | 0.75 | 8.4 | 10.6 | 2.8 |
| 425 | 0.77 | 9.9 | 14.6 | 3.7 |
| 450 | 0.79 | 11.7 | 19.1 | 5.2 |

- **The shank model reproduces the difference within 0.7 points up to 375 mph.** Without it,
  XROTOR's 16-series estimate at the shank's thickness gives a fifth to a half of the
  difference.
- **Above 400 mph the model loses too much.** At the shank the local Mach number reaches
  0.55–0.6 there, where TR-852's cylinder drag doubles. An oval of h/b 0.3–0.7 drawn as a
  cylinder of its thickness is an upper estimate; its real drag rise is later and smaller.
- **Sensitivity to the fairing.** Moving the round end and the fairing to (0.24R, 0.28R) or
  (0.28R, 0.33R) changes the difference at 375 mph from 7.4 to 5.7 or 9.2.
- **Absolute efficiency (figure 27).** The model is 0.027 below the test on the 03 (range
  0 to −0.043) and 0.043 below on the 03R. The test's 0.95 is a shaft-tension thrust
  behind a spinner designed for uniform flow.

`tests/measurement/test_shank_drag.py` pins these findings.
