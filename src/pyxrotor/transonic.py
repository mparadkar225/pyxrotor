"""Transonic lift loss (not part of XROTOR).

XROTOR's lift grows with the Prandtl-Glauert factor 1/sqrt(1 - M^2) until the local Mach
number reaches the 0.99 cap (CLCDCM), so at fixed pitch a propeller's power keeps rising
with Mach number. Real sections lose lift past a lift-divergence Mach number. Here, in the
compiled core's CLCDCM:

    M_L = MLIFT - PEN |CL - CLDMIN|        (CL: the Prandtl-Glauert lift)
    x   = W ln(1 + exp((M - M_L)/W))        (a smooth M - M_L past M_L, ~0 below it)
    CL  = CL_PG(M - x) - sign (A x + B x^2), not past zero,

with each aero section's MLIFT (the lift-divergence Mach number at CLDMIN) set by
Rotor.set_transonic_lift, and PEN, A, B fitted to the lift of NASA CR-176883's Clark Y bank
(the R&M 2036 tunnel data, pyxrotor.databank.clcdcy) at t/c 0.05-0.12, angles 2-8 deg from
zero lift and Mach 0.3-0.95: CL within 0.027 rms (0.044 above Mach 0.75, where Prandtl-
Glauert is off by 0.58). The same fit gives MLIFT = 1.090 - 3.415 t (`mlift`). The smooth
onset (W = 0.02 in Mach number) is what lets Newton converge: a kink at M_L stalls it. The stall model,
drag and moment are XROTOR's, applied to the reduced lift. With a swept blade
(pyxrotor.sweep) M and CL are the normal-plane ones. The derivatives are exact (M_L's
dependence on CL included).

No 16-series data went into the fit; `mlift` applies the Clark Y law to any section.
"""
from __future__ import annotations

PEN = 0.116      # TL_PEN in _speedups.pyx
A = 0.352        # TL_A
B = 1.005        # TL_B
W = 0.02         # TL_W
M0, M1 = 1.0904, 3.415


FLOOR = 0.5      # as the Mcrit estimate's floor; the fit covers t/c 0.05-0.12 only


def mlift(t: float, shift: float = 0.0) -> float:
    """Lift-divergence Mach number at CLDMIN for thickness ratio t (the Clark Y bank fit,
    not below FLOOR: thick root sections), plus shift."""
    return max(M0 - M1 * float(t), FLOOR) + float(shift)
