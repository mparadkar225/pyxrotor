# NACA Report 747: four-blade single- and dual-rotating propellers

D. Biermann and E. P. Hartman, *Wind-Tunnel Tests of Four- and Six-Blade Single- and
Dual-Rotating Tractor Propellers*, NACA Report 747, 1942 (Langley propeller-research
tunnel). Free from NTRS.

Two 10-ft Hamilton Standard propellers (3155-6 right-hand and 3156-6 left-hand blades,
Clark Y sections) on a 32-in streamline body, front and rear hubs 9.7 in (0.16 R) apart,
up to 550 rpm and 110 mph. The "four-blade single-rotating" propeller is the same pair of
2-blade hubs turning the same way; dual rotation turns the rear the other way, with its
blades set 1.1-1.6 deg lower (to absorb about the same power as the front at peak
efficiency). Coefficients use effective thrust (body drag added back) and both rotors'
total power.

## Files

* `geometry.json`: chord, thickness and twist from figure 1 (traced curves agree with the
  figure's blade-angle labels to 0.05 deg).
* `tests.json`: CT and CP against V/nD from figures 8, 9 (single) and 11, 13 (dual), no
  wing, digitized with `tools/digitize.py` at 0.02 steps on the attached-flow branches;
  the dual rear blade angles; peak efficiencies read from figures 7 and 10. Its `notes`
  record two corrections (crops that missed the CP = 0 or J = 0 grid line), found because
  CT J / CP came out above 1.
* `model.py`: the pyxrotor models. Both configurations are `RotorPair`s (co-rotating for
  "single"); `disk4` is one 4-blade disk for reference. Clark Y polars are estimated from
  thickness ratio, as for R&M 2357.
* `compare.py` runs the models at every digitized J (0.1 steps, ~10 min on 2 cores) into
  `results.json`; `report.py` prints the comparison and writes `summary.json`.

## Findings

Over the attached-flow range (J where both CT and CP were digitized), model / measured:

| | CT | CP | peak η, model − test |
|---|---|---|---|
| single 20–40° | −0.4 to −5.6% | −2.5 to −5.7% | −0.022 to +0.016 |
| single 45–50° | +5 to +8% | +1% | +0.025 to +0.033 |
| dual 20–30° | +6 to +7% | +3% | +0.02 to +0.04 |
| dual 35–45° | +11 to +18% | +7 to +13% | +0.02 to +0.03 |

* **The counter-rotation efficiency gain is right where it matters.** Peak efficiency,
  dual minus single (figures 10 and 7): test +0.005, −0.005, +0.001, +0.007, +0.011,
  +0.022, +0.029 at 20–50°; model +0.012, +0.013, +0.012, +0.017, +0.019, +0.023, +0.029.
  At high pitch, where the swirl is largest, the model matches to 0.001; at low pitch it
  is about 0.01 optimistic.
* **The dual propeller is loaded more in the model than in the test** at the same nominal
  blade angles, while the single one matches. Near zero thrust that is about 0.6° of
  effective pitch. Figure 12 (individual CPs, probed at a few J for 40°) shows front and
  rear absorbing nearly equal power near peak efficiency (rear/front 0.95 at J = 1.85),
  and the front below the model by about 15%. Since the circumferentially averaged
  upstream influence of the rear is the same for either rotation direction, this looks
  like the hardware's setting rather than the interaction; the report itself compares
  single and dual only after shifting the dual curves to the same zero-thrust J.
* Rear/front power in the model is 0.88 at J = 1.8 (test about 0.95), i.e. the model's
  rear is somewhat underloaded relative to the front. Slipstream contraction would move
  that further the wrong way at this spacing, so the test gives no evidence that it is
  needed here.
* Low-J CP (stalled blades) is underpredicted by up to 40%, as with R&M 2357's torque near
  stall: the section model's post-stall drag is not meant for it.

`tests/measurement/test_naca747.py` checks the data's internal consistency and pins
single-rotation CT/CP, the efficiency gain, and the known dual loading discrepancy at a
few points.
