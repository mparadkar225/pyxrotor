# ARC R&M 2357: Rotol wooden Spitfire IX propeller

A. B. Haines and P. B. Chater, *24-ft. Tunnel Tests on a Rotol Wooden Spitfire Propeller.
Test Results, and Data for Single Radius Calculations*, ARC R&M 2357, 1946. Free from
Cranfield's AERADE archive.

The propeller is Rotol hydulignum design RA.10046, tested as a 2-blade, 10.75 ft version in
the RAE 24-ft tunnel at blade angles of 10.1°, 15.1°, 22° and 27° at 0.7R. The section is
Clark Y over most of the blade.

* `tests.json`: K_T, K_Q and η against J from Tables 1–4, with the tip Mach number of each
  point. Tables 3–4 (22° and 27°) are complete for tip Mach ≈ 0.41 and 0.51; from Tables 1–2
  only the lowest-Mach row of each J group was transcribed.
* `geometry.json`: chord, thickness and blade-angle distributions digitized from Fig. 1.
  Checks against the text: t/c at 0.7R 9.4% (report 9.3%); solidity per blade 0.0335 (0.0338).

Notes:
* Two rows of Table 1 (θ = 10.1°, J = 0.302 and 0.400) have η inconsistent with
  J·K_T/(2π·K_Q) in the original; `test_rm2357.py` checks exactly these two are flagged.
* The report gives no section polars, so any XROTOR model needs estimated Clark Y sections.
  With the estimates in `tests/reference/test_ref755.py`, XROTOR 7.55 matches thrust to about
  1–3% over the normal operating range but underpredicts torque by 5–10%, so its efficiency
  is 0.04–0.06 high; near stall the torque error reaches about 35%.
