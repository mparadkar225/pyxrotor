# Propeller comparisons with sections from the CR-176883 data banks

Do the transonic propeller comparisons still need the calibrated critical Mach number of
`pyxrotor.sections` when the sections come from measured-data banks instead of thickness
estimates? `compare.py` reruns every transonic dataset with `pyxrotor.databank`'s
`clark_y_section` / `naca16_section`: XROTOR sections fitted (`fit_section`, Mach 0.3,
Re 2e6) to NASA CR-176883's Clark Y and NACA 16-series banks. Each dataset's model is
unchanged apart from its section builder. Stall (CLMAX, CLMIN) and CMCON stay the estimates'
(the banks have none), thick shank sections (Clark Y > 0.35, 16-series > 0.40) stay the
estimates', and 16-series design CLs below 0.2 take the bank's drag terms at CLD 0.2 (its
camber-drag fit starts there).

Settings: `est+0` the estimate with Mcrit 0.85 - 1.2 t (as in 7.55 comparisons), `est+cal`
the calibrated estimate, `bank+S` the bank's fitted section with its Mcrit + S.
`results.json` holds every number; `tests/measurement/test_bank_sections.py` pins the findings.

## Results
P-47 (Clark Y): mean model - test efficiency minus the low-speed bias, tip Mach 1.02-1.08 / 1.08-1.2:

| | Curtiss 714 (L-177) | HS 6507A-2 (L-594) | Curtiss 838 (L7D10) |
|---|---|---|---|
| est+0 | -0.47 / - | -0.39 / -0.82 | -0.30 / -1.05 |
| est+cal | -0.01 / -0.03 | -0.01 / +0.02 | +0.02 / +0.04 |
| bank+0 | -0.60 / -0.82 | -0.21 / -0.31 | -0.20 / -0.41 |
| bank+0.10 | -0.10 / -0.12 | -0.02 / +0.01 | -0.01 / -0.04 |
| bank+0.15 | +0.03 / +0.07 | +0.03 / +0.09 | +0.03 / +0.05 |
| low-speed bias, est+cal / bank+0.15 | 0.036 / 0.038 | 0.050 / 0.049 | 0.025 / 0.024 |

YP-39 (Clark Y, working-part efficiency), tip Mach 0.94-0.99 / 0.99-1.035 / 1.035-1.07:
est+0 -0.11/-0.30/-0.88; est+cal +0.03/0.00/-0.03; bank+0 -0.07/-0.21/-0.45;
bank+0.15 +0.04/+0.01/-0.03; bank+0.20 +0.04/+0.03/+0.02.

16-series:

| | YP-47M Curtiss 836, sonic tips (M0 0.3 J 1.0 / M0 0.5 J 1.75) | SR-2 at Mach 0.85, minus its Mach 0.6 bias | L51L28 drift at Mt 1.0 (cld 0.3, 4/8/12 deg) |
|---|---|---|---|
| est+0 | -0.18 / -0.15 | -0.72 | 2.27 / 1.54 / 2.03 |
| est+cal | 0.00 / -0.02 | -0.01 | 1.03 / 0.92 / 1.00 |
| bank+0 | -0.32 / -0.19 | -1.14 | - / 1.25 / 1.54 |
| bank+0.15 | -0.01 / -0.04 | -0.30 | 1.04 / 0.93 / 0.99 |
| bank+0.20 | 0.00 / -0.01 | -0.14 | 0.96 / 0.93 / 0.97 |
| bank+0.25 | +0.01 / +0.02 | -0.03 | 0.94 / 0.92 / 0.97 |

## Findings
- **The banks don't remove the calibration.** With their own drag-rise Mach numbers the
  bank sections lose far too much at transonic tips: on the Clark Y propellers 40-125% of
  the estimate's excess loss, on the 16-series ones more than the estimate. They need a
  shift about as large as the estimates': Clark Y + 0.10-0.20, 16-series + 0.15-0.25.
- **So the shift is a propeller-level correction**, not an error in the thickness estimate:
  the 2-D section data put the drag rise where the estimate does (Clark Y -0.01/+0.07,
  16-series -0.08/+0.05). Candidates: three-dimensional tip relief, the shape of XROTOR's
  cubic drag law above Mcrit (the banks' rise is gradual, not cubic), and the compressible
  stall limit CLMAX = (Mcrit + 0.215 - M)/0.25 + CLDMIN that a higher Mcrit also moves.
- **The required shift tracks thinness**, as in the calibrated rule: 2-4% 16-series tips
  (SR-2, L51L28) want + 0.20-0.25 on the bank, 6% (YP-47M 836) + 0.15-0.20.
- **Low speed.** On the Clark Y propellers the efficiency level barely changes (the P-47
  biases move by 0.002 at most). The 16-series bank's drag due to lift (about 8x the
  estimate's DCDCL2) raises the static power of the L51L28 blades by 6-14%, narrowing the
  static under-prediction on the cambered blades (CP 0.94/0.73 -> 1.04/0.78 of the test at
  design CL 0.3, 4/8 deg); the symmetric blade stays 15-30% low.
