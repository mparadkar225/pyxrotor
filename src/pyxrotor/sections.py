"""Section characteristics estimated from thickness ratio and design lift coefficient, as
AeroSection inputs for XROTOR's parametric section model (AERO menu).

Not a port: XROTOR takes every section characteristic from the user. These estimates are
the ones the measurement comparisons in validation/ use, collected here with a critical
Mach number calibrated against those tests.

    clark_y(t, xisect, mcrit_shift=None)         Clark Y family
    naca16(t, xisect, cld, mcrit_shift=None)     NACA 16-series of design CL cld
    mcrit(t, mcrit_shift=None)                   the critical Mach number either uses

Critical Mach number
--------------------
The thickness estimate Mcrit = 0.85 - 1.2 t (not below 0.5) with XROTOR's drag-rise law
(CDC = 10 (M - Mcrit_eff)^3, Mcrit_eff = Mcrit - 0.25 |CL - CLDMIN| - 0.0585) gives far too
much compressibility loss at transonic tips: 0.1-0.9 in efficiency at helical tip Mach
0.95-1.17 on every propeller tested. Adding a constant to each section's Mcrit fixes it:

    tips 5-7% thick (Clark Y: P-47 WR L-177, WR L-594, RM L7D10, YP-39;
                     16-series: YP-47M Curtiss 836)                       + 0.15
    tips 2-4.4% thick (16-series: RM L51L28 static, NASA SR-2 propfan)     + 0.20

So with mcrit_shift=None (the default) the shift depends on the section's own thickness:
CALIBRATED_SHIFT_THICK (0.15) at t >= CALIBRATED_T_THICK (0.05), CALIBRATED_SHIFT_THIN
(0.20) at t <= CALIBRATED_T_THIN (0.045), linear between. The data don't resolve the shape
of that step (no tested tip lies between 4.4% and 5%), and they cover helical tip Mach up
to about 1.17. tests/measurement/test_calibrated_sections.py checks the rule against every
dataset. It is a lumped correction (tip relief, the drag law's shape, the compressible
stall limit), not a section property: don't use it with measured section data.

A number for mcrit_shift gives max(0.85 - 1.2 t, 0.5) + mcrit_shift instead (0.0 is the
plain thickness estimate).

Other characteristics
---------------------
clark_y: zero-lift angle -0.37 deg per percent thickness (not below -6 deg), lift slope
6.3 - 4 (t - 0.066) per rad (not below 5.2), CLMAX 1.3 (1.1 at 20% and thicker), CLMIN
-0.5, minimum drag 0.0065 + 0.25 t^2 at CL 0.35. It reproduces R&M 2357's Clark Y sections.

naca16: zero-lift angle -cld / 2 pi, lift slope 6 per rad, CLMAX cld + 0.85 (cld + 1.1
above 12%), CLMIN cld - 0.85, minimum drag 0.0065 + 0.25 t^2 at cld.

Known gap: CLMAX doesn't fall with Mach number (XROTOR only limits it near Mcrit + 0.215),
so stall at local Mach 0.6-0.7 comes too late (YP-47M, validation/naca_yp47m).
"""
from __future__ import annotations

from .aero import AeroSection
from .constants import PI

CALIBRATED_SHIFT_THICK = 0.15
CALIBRATED_SHIFT_THIN = 0.20
CALIBRATED_T_THICK = 0.05
CALIBRATED_T_THIN = 0.045


def calibrated_shift(t: float) -> float:
    """The calibrated addition to the thickness estimate of Mcrit (see the module docstring)."""
    t = float(t)
    if t >= CALIBRATED_T_THICK:
        return CALIBRATED_SHIFT_THICK
    if t <= CALIBRATED_T_THIN:
        return CALIBRATED_SHIFT_THIN
    f = (t - CALIBRATED_T_THIN) / (CALIBRATED_T_THICK - CALIBRATED_T_THIN)
    return CALIBRATED_SHIFT_THIN + f * (CALIBRATED_SHIFT_THICK - CALIBRATED_SHIFT_THIN)


def mcrit(t: float, mcrit_shift: float | None = None) -> float:
    """Critical Mach number for thickness ratio t: max(0.85 - 1.2 t, 0.5) plus mcrit_shift,
    or plus calibrated_shift(t) when mcrit_shift is None."""
    t = float(t)
    shift = calibrated_shift(t) if mcrit_shift is None else float(mcrit_shift)
    return max(0.85 - 1.2 * t, 0.5) + shift


def clark_y(t: float, xisect: float, mcrit_shift: float | None = None) -> AeroSection:
    """Clark Y section of thickness ratio t at r/R = xisect."""
    t = float(t)
    a0 = max(-0.37 * 100.0 * t, -6.0)
    return AeroSection(
        xisect=float(xisect), a0=a0 * PI / 180.0, dclda=max(6.3 - 4.0 * (t - 0.066), 5.2),
        clmax=1.3 if t < 0.20 else 1.1, clmin=-0.5, dclda_stall=0.1, dcl_stall=0.1,
        cdmin=0.0065 + 0.25 * t * t, cldmin=0.35, dcdcl2=0.007 + 0.01 * max(t - 0.066, 0.0),
        cmcon=-0.1, mcrit=mcrit(t, mcrit_shift), reref=2.0e6, rexp=-0.15)


def naca16(t: float, xisect: float, cld: float, mcrit_shift: float | None = None) -> AeroSection:
    """NACA 16-series section of thickness ratio t and design CL cld at r/R = xisect."""
    t, cld = float(t), float(cld)
    return AeroSection(
        xisect=float(xisect), a0=-cld / (2.0 * PI), dclda=6.0,
        clmax=cld + (1.1 if t > 0.12 else 0.85), clmin=cld - 0.85, dclda_stall=0.1,
        dcl_stall=0.1, cdmin=0.0065 + 0.25 * t * t, cldmin=cld,
        dcdcl2=0.007 + 0.01 * max(t - 0.066, 0.0), cmcon=-0.25 * cld,
        mcrit=mcrit(t, mcrit_shift), reref=2.0e6, rexp=-0.15)
