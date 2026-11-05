"""pyxrotor models of the three two-blade NACA propellers of NACA RM L51L28 (Wood & Swihart
1952, static tests on the Langley 6000-hp propeller dynamometer):

    NACA 10-(0)(066)-03, 10-(3)(066)-03, 10-(5)(066)-03

10 ft diameter, NACA 16-series sections of design lift coefficient 0, 0.3 and 0.5 constant
along the blade, h/b 0.066 at 0.70R, constant chord (solidity 0.03 per blade at 0.70R). The
blade form (figure 4) is the same for all three (geometry.json). Tested at zero advance, at
beta(0.75R) = 0-16 deg and tip Mach numbers 0.25-1.02, without a spinner.

XROTOR needs a forward speed, so "static" is a small advance ratio, J_STATIC = 0.002.
CT falls about linearly with J (0.0634 at J -> 0, 0.0626 at J = 0.01 on the 10-(3) at 8 deg):
halving J to 0.001 changes CT and CP by 0.1%, so J_STATIC is within about 0.3% of static.

    static(cld, beta75, mt, mcrit_shift=0, mach_scale=1)
        CT, CP of one test point at its tip Mach number (rpm from Mt at sea level);
        mach_scale scales the Mach number by changing the speed of sound (J, Reynolds number
        kept), mach_scale=0.01 is the model without compressibility; mcrit_shift raises every
        section's critical Mach number.

Sections (estimated; RM L51L28 gives none, TN 1546 has the measured 16-series data):
  * zero-lift angle -cld/(2 pi) (the a = 1.0 uniform-load mean line), lift slope 6.0/rad;
  * CLmax = cld + 0.85 (RM L51L28: the 16-0xx and 16-5xx sections stall at about the same
    angle), 1.1 + cld inboard where t > 0.12;
  * minimum drag at CL = cld, CDmin 0.0065 + 0.25 t^2 and dCD/dCL^2 as the Clark Y estimate
    (validation/naca747/model.py);
  * Mcrit = 0.85 - 1.2 t, the same generic estimate as the P-47 comparisons
    (validation/naca_p47), so the Mcrit + 0.15 found there is tested on other blades.
XROTOR's compressibility model is the one under test: dCD = 10 (M - Mcrit_eff)^3,
compressible CLmax, Prandtl-Glauert lift slope, local Mach number capped at 0.99.

Not modeled: blade twist under load (the report says part of the rise of CT and CP with Mt at
low blade angles is probably blade twist), the shank inside 0.2R, ground and wind effects
(tests were run in winds below 400 ft/min, uncorrected).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from pyxrotor.aero import AeroSection
from pyxrotor.sections import naca16
from pyxrotor.constants import PI
from pyxrotor.rotor import Rotor

HERE = Path(__file__).resolve().parent
GEO = json.loads((HERE / "geometry.json").read_text())
FT = 0.3048
D = GEO["diameter_ft"] * FT
R = 0.5 * D
J_STATIC = 0.002
A_SL, RHO_SL, MU_SL = 340.294, 1.2250, 1.7894e-5


def section(t: float, xisect: float, cld: float, mcrit_shift: float | None = 0.0) -> AeroSection:
    """NACA 16-series section of design CL cld and thickness ratio t (pyxrotor.sections.naca16;
    mcrit_shift None is its calibrated Mcrit)."""
    return naca16(t, xisect, cld, mcrit_shift)


def sections(cld: float, mcrit_shift: float = 0.0):
    r = np.array(GEO["r_over_R"])
    t = np.array(GEO["h_over_b"])
    xs = [r[0], 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]
    return [section(float(np.interp(x, r, t)), float(x), cld, mcrit_shift) for x in xs]


def blade(beta75: float):
    """ARBI stations, 0.2R to the tip: r/R, c/R, blade angle (deg). Constant chord to 0.96R,
    then a rounded tip (c/R 0.10 at 0.98R, 0.03 at the tip); the twist beyond 0.94R is
    extrapolated along its last slope."""
    r = np.array(GEO["r_over_R"])
    c = 2.0 * np.array(GEO["b_over_D"])
    tw = np.array(GEO["twist_deg_from_0.75R"])
    rr = np.concatenate([np.linspace(r[0], 0.94, 16), [0.96, 0.98, 1.0]])
    cc = np.interp(rr, r, c)
    cc[-2], cc[-1] = 0.10, 0.03
    slope = (tw[-1] - tw[-2]) / (r[-1] - r[-2])
    tt = np.where(rr > r[-1], tw[-1] + slope * (rr - r[-1]), np.interp(rr, r, tw))
    return rr, cc, beta75 + tt


def static(cld: float, beta75: float, mt: float, mcrit_shift: float = 0.0,
           mach_scale: float = 1.0, j: float = J_STATIC) -> dict:
    """One static test point: CT, CP (and CT/CP) at tip Mach number mt (sea level)."""
    n = mt * A_SL / (PI * D)                     # test rpm: Mt = pi n D / a
    vel = j * n * D
    r = Rotor(fix_all=True)
    r.s.rho, r.s.vso, r.s.rmu = RHO_SL, A_SL / mach_scale, MU_SL
    r.set_sections(sections(cld, mcrit_shift))
    rr, cc, bb = blade(beta75)
    rr[0] = (rr[0] * R) / R                      # ARBI's XI0 = (hub R)/R, to the ulp
    r.arbi(nblds=GEO["nblds"], vel=vel, tip_radius=R, hub_radius=rr[0] * R,
           r_over_R=rr, c_over_R=cc, beta_deg=bb)
    ok = r.solve_rpm(60.0 * n)
    p = r.performance()
    ct = p["thrust"] / (RHO_SL * n * n * D ** 4)
    cp = p["power"] / (RHO_SL * n ** 3 * D ** 5)
    return {"converged": bool(ok), "CT": ct, "CP": cp, "FM": ct / cp,
            "Mt": mt * mach_scale, "rotor": r}


def tests() -> list[dict]:
    """The digitized test points (figures 6-8, faired curves sampled every 0.05 in Mt)."""
    with open(HERE / "static_fig6-8.csv") as f:
        return list(csv.DictReader(f))
