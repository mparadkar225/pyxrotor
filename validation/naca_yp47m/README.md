# XROTOR against the YP-47M altitude-wind-tunnel propeller series

NACA tested propellers on a YP-47M fuselage with stub wings in the 20 ft Cleveland altitude
wind tunnel. The report series is "Altitude-wind-tunnel investigation of performance of
several propellers on YP-47M airplane at high blade loadings", 1946. Four of its four-blade
propellers are used here, all with NACA 16-series sections:

| propeller | report | NTRS | diameter | activity factor | design CL (outboard) |
|---|---|---|---|---|---|
| Aeroproducts H20C-162-X11M2 | I | [20030063991](https://ntrs.nasa.gov/citations/20030063991) | 12 ft 7 in | 107 | 0.5 |
| Curtiss 838-1C2-18R1 | II | [20030063996](https://ntrs.nasa.gov/citations/20030063996) | 13 ft | 120 | 0.6 |
| Curtiss 732-1C2-0 | IV | [20030064146](https://ntrs.nasa.gov/citations/20030064146) | 13 ft | 150 | 0.62 |
| Curtiss 836-14C2-18R1 | V | [20030063262](https://ntrs.nasa.gov/citations/20030063262) | 13 ft | 98 | 0.5 |

**Not used.**
- Report VI (Hamilton Standard 6507A-2, [20030065899](https://ntrs.nasa.gov/citations/20030065899)) plots efficiency against J at constant CP rather than the
  cross plots used here. Its blade is the one of WR L-594 (`validation/naca_p47`), and its
  figure 1 agrees with that geometry.
- Report III (Hamilton Standard 6543A-0) is not on NTRS.

**Test conditions.**
- Free-stream Mach number 0.3–0.5, density altitude 20,000–45,000 ft, up to 2500 bhp.
- Each report plots efficiency against power coefficient at constant J and Mach number
  ("cross plots"), with helical tip Mach numbers from 0.54 to 1.03.
- The reports warn that their efficiencies show trends only. There is no tunnel-wall
  correction, and the thrust uses one average installation drag coefficient.

## Data
- **`geometry.json`.** b/D, t/c, blade angle and design CL, traced from each report's figure 1
  (`digitizing/`). The Aeroproducts and Curtiss 732 traces reproduce their quoted activity
  factors, 107.5 and 150.2. The 838 and 836 traces give 102 and 81 against the quoted 120 and
  98; see `digitizing/README.md`.
- **`efficiency_vs_cp.csv`.** 175 points from the 9 cross plots (Mach 0.4 and 0.5 for each
  propeller; Mach 0.3 for the 836), each with its J.

## Model (`model.py`)
- **Sections.** The 16-series estimate of the RM L51L28 comparison (`validation/naca_l51l28`),
  at the traced t/c and design CL. Mcrit = 0.85 − 1.2 t plus a shift. The stall CL is the
  design CL + 0.85 (+ 1.1 above 12% thickness), unless `stall_margin` sets it.
- **Hub and body.** The hub is at 0.2R. There is no body: the cowl flow isn't known.
- **Air.** ISA at 30,000 ft. Only the Reynolds number depends on it.
- **Operating point.** Each point is solved at the test's M0, J and CP. XROTOR POWE runs at
  fixed rpm, continued in Mach number from the solution without compressibility, with a
  fixed-pitch scan when POWE fails.

## Installation offset
- Away from stall and at tip Mach < 0.8, the model's efficiency is above the test's, by
  0.0–0.5. The gap grows as CP falls and J rises.
- It behaves like a constant drag area. (η_model − η_test)·CP/J³ = ΔCT/J² is nearly constant
  for each propeller, across J from 2.1 to 4.4.
- So one offset per propeller is fitted: the median over the 104 such points of the model
  without compressibility. It is removed as η_model − f·J³/CP everywhere.

| | Aeroproducts | 838 | 732 | 836 |
|---|---|---|---|---|
| f = ΔCT/J² | 0.0030 | 0.0025 | 0.0029 | 0.0014 |

- On a 13 ft disk, f = 0.0025 is a drag area of about 0.85 ft² at the dynamic pressure.
  After removing it, most tip-Mach < 0.8 curves are within ±0.04. The exceptions are
  −0.04 to −0.09, at J 2.8 (Mach 0.4) and J ≥ 4.2.
- This offset absorbs any error that scales the same way, including profile drag.
  So absolute efficiency is not tested here, only the changes with CP and Mach number.

## Results (`compare.py`, `results.json`)
Mean residual (model − offset − test) on each curve, before the breakdown (before the test
falls more than 0.08 below the curve's peak):

| propeller | M0 | J | tip Mach | 7.55 | Mcrit + 0.15 | Mcrit + 0.20 | no compressibility |
|---|---|---|---|---|---|---|---|
| 836 | 0.5 | 1.75 | 1.03 | −0.127 | −0.017 [3/4] | +0.056 | +0.102 |
| 836 | 0.3 | 1.00 | 0.99 | −0.178 | −0.001 | +0.009 | +0.019 |
| 732 | 0.5 | 2.20 | 0.87 | −0.050 | −0.006 | −0.005 | −0.004 |
| 836 | 0.5 | 2.20 | 0.87 | −0.014 | +0.038 | +0.040 | +0.044 |
| 838 | 0.4 | 1.70 | 0.84 | −0.024 | −0.010 | −0.010 | −0.007 |
| 732 | 0.4 | 1.70 | 0.84 | −0.007 | +0.005 | +0.005 | +0.006 |
| 836 | 0.4 | 1.70 | 0.84 | −0.058 | −0.017 | −0.016 | −0.013 |
| Aeroproducts | 0.5 | 2.40 | 0.82 | +0.016 | +0.029 | +0.029 | +0.028 |
| Aeroproducts | 0.4 | 1.75 | 0.82 | +0.002 | +0.013 | +0.013 | +0.014 |

- **At sonic tips, Mcrit + 0.15 is right again.** On the Curtiss 836 at tip Mach 0.99–1.03,
  7.55 loses 0.13–0.18 too much.
  - Mcrit + 0.15 is within 0.02 on both curves.
  - At tip Mach 1.03, Mcrit + 0.20 is 0.056 high and no compressibility is 0.10 high.
  - This is the P-47 and YP-39 correction. Thinner 16-series blades (RM L51L28, SR-2) wanted
    + 0.20, but here the tips are 6% thick.
- **At tip Mach 0.82–0.87 the settings differ by less than 0.06, and 7.55 is a little low.**
  Over these seven curves 7.55 averages −0.019 and Mcrit + 0.15 +0.007.
- **Stall breakdown is missed.** At low J and high CP, the test's efficiency drops by 0.2–0.4
  within ΔCP ≈ 0.1. The reports call this blade stall, which their slipstream surveys show
  at the tips. The tip Mach number is only 0.72–0.84.
  - With its stall CL of design CL + 0.85 (≈ 1.35–1.45) the model stays near peak: it is
    0.04–0.25 high past the breakdown (Aeroproducts J 1.75 and 2.15: +0.22 to +0.25; Curtiss:
    +0.04 to +0.10).
  - A stall CL of design CL + 0.5 (`stall_margin=0.5`, "Mcrit + 0.15, stall cld+0.5" in
    `results.json`) puts the Aeroproducts breakdown right (−0.03 and +0.02).
  - The same setting drops the 838 and 836 curves 0.09–0.45 too far. The 732 barely
    changes.
  - No single stall CL fits all four propellers. At local Mach 0.6–0.7, 16-series sections
    stall well below XROTOR's default CLMAX. XROTOR reduces CLMAX with Mach number only
    near Mcrit + 0.215.
- **Curtiss 836 at Mach 0.3, J 1.0–1.3**, tip Mach 0.78–0.99: the efficiency falls steadily
  from 0.83 to 0.2 as CP rises from 0.1 to 1.0.
  - Before the breakdown, all settings except 7.55 are within 0.02.
  - Past it, Mcrit + 0.15 converges on 2 of 5 points at J 1.0, and is within 0.04 on J 1.3.

`tests/measurement/test_yp47m.py` pins these findings.

## Limits
- The test efficiencies have no wall correction. One installation offset per propeller is
  fitted, so only trends are compared.
- The cross-plot points were assigned to curves by eye, and symbols were read to about ±0.5
  points.
- The 838 and 836 geometry disagrees with the quoted activity factors by about 18%.
- No body or cowl flow, and the hub is at 0.2R.
- Report V's figure 5 legend misprints one J; it is taken as 1.30.
