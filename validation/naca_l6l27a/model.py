"""pyxrotor model of the two full-scale NACA propellers of RM L6L27a (Maynard 1947, NTRS
20050019371). They differ only in their shanks:

    "NACA 10-(5)(08)-03"    airfoil (NACA 16-series) sections to the spinner, wide and thin
    "NACA 10-(5)(08)-03R"   conventional round shanks, faired into the blade by about 0.3R

Both are 10 ft, two blades, NACA 16-series, design CL 0.5, 8% thick at 0.7R, solidity 0.03
per blade at 0.7R. They were tested on a 2000 hp dynamometer in the Langley 16-ft high-speed
tunnel, with a spinner of 0.217 D (radius 0.217 R).

Geometry: geometry.json, read from figure 4. Sections: pyxrotor.sections.naca16 at each
station's t/c and design CL, with the calibrated Mcrit by default. The 03R's shank is a
thick oval (figure 7: h/b 0.73 at 0.2R, inside the spinner, and 0.295 at 0.3R). It is
modeled with pyxrotor.shank: TR-852's cylinder drag on the oval's thickness to
r_round = 0.26, faired linearly into a 16-series section by r_airfoil = 0.30.

The comparison is at the test's constant power coefficient and 1140 rpm (figures 27-28),
solved with pyxrotor.robust (and RoundShank.solve_power for the 03R). The speed of sound
is 342 m/s, from figure 28's helical tip Mach numbers at 1140 rpm (341-344 m/s from 225
to 450 mph). Density 1.2 kg/m^3 is assumed; the report doesn't give it, and only the
Reynolds number depends on it.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from pyxrotor import robust
from pyxrotor.constants import PI
from pyxrotor.rotor import Rotor
from pyxrotor.sections import naca16
from pyxrotor.shank import RoundShank

HERE = Path(__file__).resolve().parent
GEO = json.loads((HERE / "geometry.json").read_text())
P03, P03R = "NACA 10-(5)(08)-03", "NACA 10-(5)(08)-03R"
D = GEO["diameter_ft"] * 0.3048
R = 0.5 * D
HUB = GEO["spinner_r_over_R"]
A, RHO, MU = 342.0, 1.2, 1.8e-5
RPM = 1140.0
SHANK = (0.26, 0.30)


def sections(prop: str, mcrit_shift: float | None = None, shank: RoundShank | None = None):
    g = GEO[prop]
    r, t, c = (np.array(g[k]) for k in ("r_over_R", "h_over_b", "cl_design"))
    xs = [HUB, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]
    out = []
    if shank is not None:
        out = [shank.round_section(HUB), shank.round_section(shank.r_round)]
        xs = [shank.r_airfoil] + [x for x in xs if x > shank.r_airfoil + 1e-9]
    for x in xs:
        out.append(naca16(float(np.interp(x, r, t)), float(x), float(np.interp(x, r, c)), mcrit_shift))
    return out


def rotor(prop: str, J: float, mcrit_shift: float | None = None, shank: tuple | None = SHANK,
          mach_scale: float = 1.0, beta75: float | None = None):
    """A Rotor at advance ratio J and RPM. shank applies to the 03R only (None: no shank
    model, its root given 16-series sections at the plotted thickness, as XROTOR would)."""
    g = GEO[prop]
    n = RPM / 60.0
    vel = J * n * D
    sh = None
    if shank and prop == P03R:
        t = g["shank_h_over_b"]
        sh = RoundShank(*shank, thickness_ratio=(t["r_over_R"], t["h_over_b"]))
    r = Rotor(fix_all=True)
    r.s.rho, r.s.vso, r.s.rmu = RHO, A / mach_scale, MU
    r.set_sections(sections(prop, mcrit_shift, sh))
    rr = np.array(g["r_over_R"], dtype=float)
    rr[0] = (rr[0] * R) / R
    if beta75 is None:
        beta75 = math.degrees(math.atan(J / (0.75 * PI))) + 8.0
    r.arbi(nblds=GEO["nblds"], vel=vel, tip_radius=R, hub_radius=rr[0] * R, r_over_R=rr,
           c_over_R=2.0 * np.array(g["b_over_D"]),
           beta_deg=beta75 + np.array(g["twist_deg_from_0.75R"]))
    return r, sh


def run(prop: str, J: float, cp: float, mcrit_shift: float | None = None,
        shank: tuple | None = SHANK, mach_scale: float = 1.0) -> dict:
    """Efficiency at power coefficient cp, advance ratio J and 1140 rpm."""
    r, sh = rotor(prop, J, mcrit_shift, shank, mach_scale)
    n = RPM / 60.0
    power = cp * RHO * n ** 3 * D ** 5
    sol = (sh.solve_power(r, power, RPM) if sh is not None
           else robust.solve_power(r, power, RPM))
    vel = J * n * D
    return {"prop": prop, "J": J, "CP": cp, "usable": sol.usable, "method": sol.method,
            "eta": sol.efficiency, "beta75": sol.beta75_deg, "M": vel / A,
            "Mt": vel / A * math.sqrt(1.0 + (PI / J) ** 2) * mach_scale,
            "shank_power_fraction": (sh.drag(r)["power"] / sol.power) if (sh and sol.usable) else 0.0}


def table(name: str) -> list[dict]:
    with open(HERE / name) as f:
        return list(csv.DictReader(f))
