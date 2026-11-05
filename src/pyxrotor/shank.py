"""Round blade shanks: the drag of a propeller's cylindrical root, which XROTOR's airfoil
section model can't represent.

Not a port. XROTOR's parametric drag (CLCDCM: CD = CDMIN + DCDCL2 (CL - CLDMIN)^2 plus
10 (M - Mcrit_eff)^3) has no form for a circular cylinder. A cylinder has no lift, and in
the flight tests of NACA TR-852 (Barlow 1946, NTRS 19930091924: a 4.5 in cylinder under an
XP-51 wing, full-scale Reynolds number 5e5-1.1e6) its section drag coefficient (on the
diameter) is 0.51-0.58 at Mach 0.30-0.48. It then rises steeply to 1.09 at Mach 0.70 and
falls to 0.86 at 0.83. No setting of the cubic law follows that.

So the round part of the blade is given to XROTOR as a placeholder section with no lift
and no drag (round_section), and its drag is added from the TR-852 curve outside XROTOR:

    shank = RoundShank(r_round=0.26, r_airfoil=0.35)     # thickness_ratio for an oval shank
    rotor.set_sections([shank.round_section(hub), shank.round_section(0.26),
                        *airfoil sections from 0.35 outward])
    rotor.arbi(...)            # chord at the round stations = the shank diameter
    shank.drag(rotor)          # thrust, torque, power of the round shank (W, N, N m)
    shank.performance(rotor)   # rotor.performance() with the shank drag included
    shank.solve_power(rotor, power, rpm)   # constant speed, the shank's power included

Why this is exact for XROTOR's equations:
* A station's drag enters only TPQ's viscous thrust and power, DTV = -1/2 W c CD SI B dXI and
  DPV = 1/2 W c CD CI B dXI XI/ADV. APITER's residual (circulation against CL) does not
  depend on CD.
* So adding w * CD_cyl(M) at each station afterwards, with the station's own W, CI, SI from
  CSCALC, gives the thrust and power XROTOR would get if CLCDCM had returned that CD.
* w is the station's share of placeholder sections in GETCLCDCM's linear blend between
  sections: 1 inboard of r_round, falling linearly to 0 at r_airfoil. In between, the
  section is the airfoil's CL and CD times (1 - w) plus the cylinder's drag times w: a
  linear fairing.
* At fixed pitch nothing else changes. At fixed power, solve_power iterates XROTOR's power
  target (POWE at fixed rpm, through pyxrotor.robust) until the blades plus shank absorb the
  requested power.
tests/model/test_shank.py checks the formulas against XROTOR's own viscous totals (TVIS,
PVIS).

Limits:
* The curve is for M 0.30-0.83 at Reynolds number 5e5-1.1e6. Below Mach 0.30 it is held at
  0.58. At lower Reynolds numbers a real cylinder is subcritical (CD about 1.2), so low-speed
  shank drag at small scale is underestimated.
* TR-852's wake surveys gave up to 50% more drag than its force and pressure tests. The
  force tests are used here.
* The spanwise Mach gradient of a real shank and its interference with the spinner aren't
  in the test.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import robust
from .aero import AERODATA_ROW, AeroSection
from .wake import cscalc

# TR-852 figure 11, force tests (faired), section drag coefficient on the diameter
TR852_CYLINDER = ((0.30, 0.58), (0.35, 0.545), (0.40, 0.52), (0.45, 0.51), (0.48, 0.515),
                  (0.50, 0.555), (0.55, 0.765), (0.575, 0.87), (0.60, 0.95), (0.64, 1.04),
                  (0.70, 1.09), (0.75, 1.08), (0.80, 0.98), (0.83, 0.86))
PLACEHOLDER_MCRIT = 10.0


def cylinder_cd(mach: float, curve=TR852_CYLINDER) -> float:
    """Section drag coefficient of a round shank at local Mach number `mach` (TR-852),
    held at the end values outside the tested range."""
    m, c = zip(*curve)
    return float(np.interp(float(mach), m, c))


def round_section(xisect: float) -> AeroSection:
    """A section XROTOR treats as having no lift and no drag: lift slope 0.01 per radian
    (stall disabled, as the post-stall slope equals it), CD 0, Mcrit 10."""
    return AeroSection(xisect=float(xisect), a0=0.0, clmax=1.0, clmin=-1.0, dclda=0.01,
                       dclda_stall=0.01, dcl_stall=0.1, cdmin=0.0, cldmin=0.0, dcdcl2=0.0,
                       cmcon=0.0, mcrit=PLACEHOLDER_MCRIT, reref=1.0e6, rexp=0.0)


def _is_placeholder(s, n: int) -> bool:
    return (s.aerodata[AERODATA_ROW["cdmin"] - 1, n] == 0.0
            and s.aerodata[AERODATA_ROW["mcrit"] - 1, n] >= PLACEHOLDER_MCRIT)


@dataclass
class RoundShank:
    """A bluff shank to r_round (r/R), faired linearly into the airfoil by r_airfoil.

    thickness_ratio: the shank's thickness over its chord (h/b) against r/R, as
    (r_over_R, h_over_b) arrays; None is a round shank (h/b = 1). The cylinder's drag is
    taken on the frontal thickness h, so on the chord it is CD_cyl * h/b. That treats an
    oval shank as a cylinder of its thickness, an upper estimate."""
    r_round: float
    r_airfoil: float
    curve: tuple = field(default=TR852_CYLINDER)
    thickness_ratio: tuple | None = None

    def hb(self, xi: float) -> float:
        if self.thickness_ratio is None:
            return 1.0
        r, t = self.thickness_ratio
        return float(np.interp(xi, r, t))

    def round_section(self, xisect: float) -> AeroSection:
        return round_section(xisect)

    def weights(self, rotor) -> np.ndarray:
        """Each station's share of placeholder (round) sections in GETCLCDCM's blend."""
        s = rotor.s
        xa = np.array(s.xiaero[:s.naero], dtype=float)
        ph = np.array([_is_placeholder(s, n) for n in range(s.naero)], dtype=float)
        out = np.zeros(s.ii)
        for i in range(s.ii):
            xi = float(s.xi[i])
            if s.naero == 1 or xi <= xa[0]:
                out[i] = ph[0]
                continue
            n = int(np.searchsorted(xa, xi, side="right")) - 1      # XIAERO(N) <= XI
            if n >= s.naero - 1:
                out[i] = ph[-1]
                continue
            frac = (xi - xa[n]) / (xa[n + 1] - xa[n])
            out[i] = (1.0 - frac) * ph[n] + frac * ph[n + 1]
        return out

    def station_drag(self, rotor, cd=None):
        """Per-station viscous thrust and power (TPQ's DTV, DPV; nondimensional) for drag
        coefficients cd (default: w * cylinder CD at the local Mach number), and the
        local Mach numbers."""
        s = rotor.s
        blds = float(s.nblds)
        w = self.weights(rotor) if cd is None else None
        dtv, dpv, mach = np.zeros(s.ii), np.zeros(s.ii), np.zeros(s.ii)
        for i in range(s.ii):
            v = cscalc(s, i)
            mach[i] = v.w * s.vel / s.vso
            if cd is None:
                c = w[i] * cylinder_cd(mach[i], self.curve) * self.hb(float(s.xi[i]))
            else:
                c = float(cd[i])
            if c == 0.0:
                continue
            bdx = blds * float(s.dxi[i])
            hwc = 0.5 * v.w * float(s.ch[i])
            xx = float(s.xi[i]) / s.adv
            dtv[i] = -hwc * c * v.si * bdx
            dpv[i] = hwc * c * v.ci * bdx * xx
        return dtv, dpv, mach

    def drag(self, rotor) -> dict:
        """The shank's thrust (negative), power and torque, dimensional."""
        s = rotor.s
        dtv, dpv, mach = self.station_drag(rotor)
        r2 = s.rad * s.rad
        thrust = float(np.add.accumulate([0.0, *dtv])[-1]) * (s.rho * s.vel * s.vel * r2)
        power = float(np.add.accumulate([0.0, *dpv])[-1]) * (s.rho * s.vel ** 3 * r2)
        omega = s.vel / (s.rad * s.adv)
        return {"thrust": thrust, "power": power, "torque": power / omega,
                "max_mach": float(mach[self.weights(rotor) > 0].max(initial=0.0))}

    def performance(self, rotor) -> dict:
        """rotor.performance() with the shank's drag added to thrust, power and torque."""
        p = dict(rotor.performance())
        d = self.drag(rotor)
        p["thrust_blades"], p["power_blades"] = p["thrust"], p["power"]
        p["thrust"] += d["thrust"]
        p["power"] += d["power"]
        p["torque"] += d["torque"]
        p["efficiency"] = p["thrust"] * rotor.s.vel / p["power"]
        en, dia = p["rpm"] / 60.0, 2.0 * rotor.s.rad
        p["CT"] = p["thrust"] / (rotor.s.rho * en * en * dia ** 4)
        p["CP"] = p["power"] / (rotor.s.rho * en ** 3 * dia ** 5)
        p["shank"] = d
        return p

    def solve_power(self, rotor, power: float, rpm: float, iterations: int = 10,
                    rtol: float = 1e-5, **kw) -> robust.Solution:
        """Blade angle at which blades plus shank absorb `power` (W) at `rpm`: repeats
        robust.solve_power with the blades' target reduced by the shank's power. Returns
        that Solution with thrust, power, torque, efficiency, CT and CP including the shank.
        For a "scan" solution the shank's drag is taken at the lower bracketing solve."""
        power, rpm = float(power), float(rpm)
        start = rotor.s.copy()
        p_shank = 0.0
        sol = None
        for _ in range(iterations):
            rotor.s = start.copy()
            sol = robust.solve_power(rotor, power - p_shank, rpm, **kw)
            if not sol.usable:
                return sol
            new = self.drag(rotor)["power"]
            done = abs(new - p_shank) <= rtol * power
            p_shank = new
            if done:
                break
        d = self.drag(rotor)
        sol.thrust += d["thrust"]
        sol.power += d["power"]
        sol.torque += d["torque"]
        sol.efficiency = sol.thrust * rotor.s.vel / sol.power
        en, dia = rpm / 60.0, 2.0 * rotor.s.rad
        sol.CT = sol.thrust / (rotor.s.rho * en * en * dia ** 4)
        sol.CP = sol.power / (rotor.s.rho * en ** 3 * dia ** 5)
        sol.exact = sol.exact and math.isclose(sol.power, power, rel_tol=10 * rtol)
        return sol
