"""pyxrotor model of the four-blade propellers tested on a YP-47M in the NACA Cleveland
altitude wind tunnel (1946; "Altitude-wind-tunnel investigation of performance of several
propellers on YP-47M airplane at high blade loadings", parts I, II, IV, V):

    aeroproducts_H20C-162   NTRS 20030063991   12 ft 7 in
    curtiss_838-1C2-18R1    NTRS 20030063996   13 ft
    curtiss_732-1C2-0       NTRS 20030064146   13 ft
    curtiss_836-14C2-18R1   NTRS 20030063262   13 ft

All have NACA 16-series sections. Geometry (geometry.json) is traced from each report's
figure 1 (digitizing/). Sections are pyxrotor.sections.naca16 at the traced t/c and
design CL: Mcrit = 0.85 - 1.2 t plus a shift (or the calibrated rule, mcrit_shift=None).

The tests give efficiency against power coefficient at constant J and free-stream Mach
number (0.3-0.5). A point is reproduced at the test's M0, J and CP. The air is ISA at
30,000 ft (the tests ran at density altitudes of 20,000-45,000 ft; only the Reynolds number
depends on it). No body: the YP-47M cowl's velocity field isn't known.

    run(prop, M0, J, cp, mcrit_shift=0, mach_scale=1, stall_margin=None)
        The blade angle for the test's power at fixed rpm, from pyxrotor.robust.solve_power
        (XROTOR POWE; if that fails, continued in Mach number from the solution without
        compressibility; if that fails, interpolated between fixed-pitch solves).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from pyxrotor import robust
from pyxrotor.constants import PI
from pyxrotor.sections import naca16
from pyxrotor.rotor import Rotor

HERE = Path(__file__).resolve().parent
GEO = json.loads((HERE / "geometry.json").read_text())

INCH = 0.0254
HUB = 0.20
GAMMA, R_AIR = 1.4, 287.05
T_ALT, RHO_ALT = 228.71, 0.45831            # ISA, 30,000 ft


def sections(prop: str, mcrit_shift: float = 0.0, stall_margin: float | None = None):
    """16-series estimates at the traced t/c and design CL; stall_margin sets CLMAX to
    design CL + stall_margin (the estimate's own is + 0.85, + 1.1 above 12% thick)."""
    g = GEO[prop]
    r, t, c = (np.array(g[k]) for k in ("r_over_R", "h_over_b", "cl_design"))
    out = []
    for x in (HUB, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0):
        cld = float(np.interp(x, r, c))
        s = naca16(float(np.interp(x, r, t)), float(x), cld, mcrit_shift)
        if stall_margin is not None:
            s.clmax = cld + stall_margin
        out.append(s)
    return out


def blade(prop: str, beta75: float):
    g = GEO[prop]
    r = np.array(g["r_over_R"])
    m = r >= HUB - 1e-9
    return r[m].copy(), 2.0 * np.array(g["b_over_D"])[m], beta75 + np.array(g["twist_deg_from_0.75R"])[m]


def _rotor(prop, M0, J, beta, mcrit_shift, mach_scale, stall_margin=None):
    a = math.sqrt(GAMMA * R_AIR * T_ALT)
    mu = 1.458e-6 * T_ALT ** 1.5 / (T_ALT + 110.4)
    D = GEO[prop]["D_in"] * INCH
    R = 0.5 * D
    vel = M0 * a
    n = vel / (J * D)
    r = Rotor(fix_all=True)
    r.s.rho, r.s.vso, r.s.rmu = RHO_ALT, a / mach_scale, mu
    r.set_sections(sections(prop, mcrit_shift, stall_margin))
    rr, cc, bb = blade(prop, beta)
    rr[0] = (rr[0] * R) / R                       # ARBI's XI0 = (hub R)/R, to the ulp
    r.arbi(nblds=GEO[prop]["nblds"], vel=vel, tip_radius=R, hub_radius=rr[0] * R,
           r_over_R=rr, c_over_R=cc, beta_deg=bb)
    return r, n, D


def run(prop: str, M0: float, J: float, cp: float, mcrit_shift: float | None = 0.0,
        mach_scale: float = 1.0, stall_margin: float | None = None) -> dict:
    """One point at the test's M0, J and CP (pyxrotor.robust.solve_power at fixed rpm).
    mcrit_shift None is pyxrotor.sections' calibrated Mcrit."""
    # start from the advance angle at 0.75R plus 5 deg; POWE from far below it can settle
    # on a negative-thrust branch
    guess = math.degrees(math.atan(J / (0.75 * PI))) + 5.0
    r, n, D = _rotor(prop, M0, J, guess, mcrit_shift, mach_scale, stall_margin)
    sol = robust.solve_power(r, cp * RHO_ALT * n ** 3 * D ** 5, 60.0 * n)
    return {"prop": prop, "M0": M0, "J": J, "CP": cp, "method": sol.method, "usable": sol.usable,
            "eta": sol.efficiency, "CT": sol.thrust / (RHO_ALT * n * n * D ** 4),
            "beta75": sol.beta75_deg, "Mt": M0 * mach_scale * math.sqrt(1.0 + (PI / J) ** 2)}


def points() -> list[dict]:
    with open(HERE / "efficiency_vs_cp.csv") as f:
        return list(csv.DictReader(f))
