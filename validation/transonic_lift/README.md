# The transonic lift loss on the P-47 flight tests

`pyxrotor.transonic` replaces XROTOR's Prandtl-Glauert lift past a lift-divergence Mach
number M_L = MLIFT − 0.116 |CL − CLDMIN|. The lift is the Prandtl-Glauert lift at M − x
less 0.352 x + 1.005 x², with x a smooth M − M_L (width 0.02). The law and
MLIFT = 1.090 − 3.415 t are fitted to the lift of CR-176883's Clark Y bank (the R&M 2036
tunnel data). `compare.py` runs the three P-47 propellers (`validation/naca_p47`) with it.

Mean (model − test) efficiency minus the low-speed bias, tip Mach 1.02–1.08 / 1.08–1.2,
with unconverged runs in brackets:

| setting | Curtiss 714 (L-177) | HS 6507A-2 (L-594) | Curtiss 838 (L7D10) |
|---|---|---|---|
| plain Mcrit | −0.47 / — [8] | −0.39 / −0.82 [9] | −0.30 / −1.05 [5] |
| plain, lift loss | −0.44 / — [8] | −0.41 / −0.78 [5] | −0.26 / −0.98 [6] |
| calibrated (+0.15) | −0.01 / −0.03 [0] | −0.01 / +0.02 [1] | +0.02 / +0.04 [0] |
| calibrated, lift loss +0 | −0.02 / −0.09 [1] | −0.05 / −0.12 [0] | −0.03 / −0.11 [0] |
| calibrated, lift loss +0.15 | +0.00 / −0.01 [3] | −0.03 / −0.07 [2] | +0.01 / +0.01 [3] |

- **The lift loss doesn't replace the calibration.** On the plain Mcrit it moves the losses by
  0.01–0.08: the drag rise dominates.
- **With the calibrated drag, the lift loss needs the same shift.** The unshifted lift
  divergence (from the 2-D data) loses 0.09–0.12 too much above tip Mach 1.08. Shifted by the
  same +0.15 it stays within 0.07 (better on the Curtiss 838, worse on the 6507A-2). So the
  propeller-level correction acts on the lift divergence as much as on the drag rise, as
  tip relief would.
- **Convergence.** The smooth onset is needed. With a kink at M_L, 11–41 runs per
  propeller didn't converge; now 0–3 don't. Far past the local-Mach cap the loss makes solves
  converge that 7.55's model can't (`tests/model/test_transonic.py`).

`tests/measurement/test_transonic_lift.py` pins these findings.
