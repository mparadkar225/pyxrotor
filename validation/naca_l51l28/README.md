# XROTOR's compressibility model against static tests of thin 16-series propellers

NACA RM L51L28, J. H. Wood and J. M. Swihart, *The effect of blade-section camber on the static
characteristics of three NACA propellers*, 1952 ([NTRS 19930086559](https://ntrs.nasa.gov/citations/19930086559)).

| propeller | sections | design CL |
|---|---|---|
| NACA 10-(0)(066)-03 | NACA 16-series, constant from shank to tip | 0 |
| NACA 10-(3)(066)-03 | NACA 16-series, constant from shank to tip | 0.3 |
| NACA 10-(5)(066)-03 | NACA 16-series, constant from shank to tip | 0.5 |

- **The propellers.** Two blades, 10 ft diameter, h/b 0.066 at 0.70R (0.044 at 0.94R),
  constant chord (b/D 0.0666), solid duralumin.
- **The tests.** They were run on the Langley 6000-hp dynamometer's static stand, without a
  spinner. Blade angles were 0–16° at 0.75R, at 600–2200 rpm, which is tip Mach 0.25–1.02.
  At 16° the blades fluttered before full speed.
- **Why these tests.** They differ from the P-47 propellers in several ways: thin
  high-critical-Mach sections, two blades, zero advance and a test-stand measurement.
  Thrust and torque come from balances, not a wake survey.

## Data
- **`static_fig6-8.csv`.** The faired curves of figures 6–8 (C_T and C_P against tip Mach
  number at each blade angle), sampled every 0.05 in M_t.
  - `traces_fig6-8.json` holds the traces themselves.
  - Spot checks against the cross-plots of figures 12–14 agree within about 0.002 in C_T.
- **`geometry.json`.** Figure 4, the blade form, which is the same for all three propellers.
- **`digitizing/`.** The scripts and how they were run. Two things mattered:
  - Each page's grid is skewed by up to 0.36°, and not uniformly, so every grid line is
    fitted separately.
  - One stretch where a curve runs along a grid line was read by eye.
- **The report's text.** Its OCR, in naca-archive's `data` branch, gives the test
  conditions. Of note:
  - wind below 400 ft/min, uncorrected;
  - the spinner's effect was negligible in earlier tests;
  - the rise of C_T and C_P with M_t at low blade angles is "in part attributable to Mach
    number effect and a small part is probably caused by blade twist";
  - the 10-(0) twists 1.5 times as much as the 10-(5) (RM L51A19).

## Model (`model.py`)
- **Static operation.** XROTOR at fixed rpm with J = 0.002, within about 0.3% of the static
  limit.
- **Sections.** Estimated NACA 16-series characteristics:
  - zero-lift angle −c_ld/2π, the a = 1.0 mean line;
  - minimum drag at CL = c_ld;
  - CLmax = c_ld + 0.85.
- **Critical Mach number.** The generic Mcrit = 0.85 − 1.2 t of the P-47 comparisons, 0.80
  at the tip, so that their Mcrit + 0.15 is tested on different blades.
- **Not modeled:** blade twist under load, and the shank inside 0.2R.

## Results (`compare.py`; `results.json` at the test's blade angle, `results_matched.json` at its thrust)
At the test's blade angle the model can't follow the test.

**From M_t 0.4 to 1.0 the test's C_T rises 23–52% and its C_P 24–63%.** The model with
Mcrit + 0.20 gives +6–18% and +12–32%. The report puts part of this on blade twist under
load, and the test's rise is largest on the symmetric blade, which twists most. So the
comparison is made at the test's thrust instead: the model's blade angle is swept, and its
C_P is taken at the test's C_T.

The table gives model C_P / test C_P at the test's C_T. The columns at M_t 0.4 and 0.8 are
the ratio itself. The columns at 0.9 and 1.0 are its drift: the ratio divided by its value at
M_t 0.8, where 1.00 means the model's power rise above M_t 0.8 is the test's. A "—" means no
converged solution reaches the test's thrust.

| design CL, β | at 0.4 | at 0.8 | 7.55 drift, 0.9 / 1.0 | Mcrit + 0.10 | **Mcrit + 0.15** | **Mcrit + 0.20** | no compressibility |
|---|---|---|---|---|---|---|---|
| 0, 4° | 0.76 | 0.78 | 1.23 / 2.48 | 1.01 / 1.17 | 1.00 / 1.14 | 1.00 / 1.03 | 1.00 / 0.97 |
| 0, 8° | 0.62 | 0.69 | 1.36 / 2.75 | 1.02 / 1.42 | 1.00 / 1.14 | 1.00 / 1.03 | 1.00 / 1.02 |
| 0, 12° | 0.58 | 0.61 | 1.64 / — | 1.16 / 2.05 | 1.06 / 1.52 | 1.03 / 1.26 | 1.05 / — |
| 0.3, 4° | 0.94 | 0.98 | 1.44 / 2.27 | 1.05 / 1.45 | 0.99 / 1.15 | 0.98 / 1.01 | 0.97 / 0.94 |
| 0.3, 8° | 0.73 | 0.80 | 1.10 / 1.54 | 0.99 / 1.01 | 0.99 / 0.94 | 0.99 / 0.92 | 0.99 / 0.93 |
| 0.3, 12° | 0.69 | 0.68 | 1.30 / 2.03 | 1.08 / 1.20 | 1.06 / 1.04 | 1.06 / 0.97 | 1.05 / 0.98 |
| 0.5, 4° | 0.94 | 0.97 | 1.47 / 2.12 | 1.09 / — | 1.00 / 1.19 | 0.96 / 1.00 | 0.95 / 0.88 |
| 0.5, 8° | 0.86 | 0.87 | 1.13 / 1.40 | 1.00 / 1.02 | 0.98 / 0.94 | 0.98 / 0.91 | 0.98 / 0.90 |
| 0.5, 12° | 0.75 | 0.72 | 1.15 / 1.48 | 1.03 / 1.00 | 1.02 / 0.92 | 1.02 / 0.90 | 1.03 / 0.92 |

- **Static bias.** At low tip Mach number the model needs too little power for the test's
  thrust:
  - 38–42% too little on the symmetric blade at 8–12°;
  - 25–31% on the cambered blades at 12°, and 14–27% at 8°;
  - within 6% at 4° on the cambered blades.
  - At zero advance a propeller's inner blade is deeply stalled: at 12° the 0.2R section is
    at a blade angle of 42°.
  - XROTOR's post-stall drag, 2·(excess angle)², is about 0.2–0.5 there, well below a
    stalled section's.
  - This is a limit of the static case and the stall model, not of compressibility. The
    drift columns divide it out.
- **7.55 as configured overpredicts the losses here too.**
  - Power at the test's thrust rises 1.1–1.6 times too much by M_t 0.9 and 1.4–2.8 times by
    1.0.
  - The test's power at matched thrust hardly changes up to M_t 0.9.
  - This is the P-47 result again, on blades of a different family, count and operating
    condition.
- **Mcrit + 0.15, the P-47 correction.**
  - It fits the cambered blades at 8–12°, near their design lift coefficient, within 8%.
  - It overpredicts the rise by 14–19% at 4° on all three blades, and by 52% on the
    symmetric blade at 12°. In those cases the outer sections work far from c_ld (CL ≈ 0 on
    the cambered blades at 4°). XROTOR then lowers the effective Mcrit by
    0.25·|CL − c_ld|, up to 0.125 here, which these sections don't seem to show.
- **Mcrit + 0.20 is within 10% everywhere but the symmetric blade at 12° (+26%).** These
  thin 16-series blades want about 0.05 more shift than the P-47's Clark Y-like sections:
  - The 0.85 − 1.2 t estimate is low for 16-series sections, which were designed for a high
    critical Mach number.
  - Static operation also differs from flight: the tip sections see the full tip speed, but
    with large induced inflow and stronger three-dimensional relief at the tip.
- **Without compressibility the model misses a real loss.** At 4°, power falls 3–12%
  behind the test's by M_t 1.0, so the test does show a compressibility power rise there,
  and Mcrit + 0.20 matches it.

`tests/measurement/test_l51l28_static.py` pins these findings.

## Limits
- The section characteristics are estimates. NACA TN 1546 has measured 16-series data at
  Mach 0.3–0.8 and would replace them.
- Blade twist under load isn't modeled. RM L51A19 has the measured deflections.
- The matched-thrust comparison assumes twist moves the operating point along the
  propeller's own C_P–C_T curve. A twist that is not uniform along the blade would shift
  that curve a little.
- The pressure-distribution reports on the same three propellers (RM L50C03, L9L12, L50B21,
  in the Langley 16-ft tunnel at forward speed) would test the section loads directly.
