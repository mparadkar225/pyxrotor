"""pyxrotor models of three propellers flight-tested on the Republic P-47 at Langley, for
checking XROTOR's compressibility model against flight data:

    "l177"   Curtiss No. 714-1C2-12, four blades, 12 ft 2 in, Clark Y (NACA ACR L4L07 /
             WR L-177, Vogeley 1944; P-47C)
    "l594"   Hamilton Standard No. 6507A-2, three blades, 12 ft 10 in, design CL 0.4
             (NACA MR L5G09a / WR L-594, Gardner 1945; P-47D)
    "l7d10"  Curtiss No. 838-1C2-18, three blades, 12 ft 10 in, design CL 0.5, with
             trailing-edge extensions (NACA RM L7D10, Gardner 1947; P-47D-28)

Each test run (the reports' tables, transcribed in <prop>_table*.csv) is reproduced at its own
J, power coefficient, rotational speed, density and airplane Mach number: XROTOR POWE at fixed
rpm sets the blade angle to absorb the measured power, as the propeller's governor did, and
the efficiency, thrust coefficient and blade angle are compared with the measured ones.

    run(row, mach_scale=1)  one test run; mach_scale scales the airplane (and tip) Mach
                            number by changing the speed of sound, not V or n, so J, CP and
                            Reynolds number stay the test's: mach_scale=0.01 is the model
                            without compressibility, mcrit_shift raises every section's
                            critical Mach number
    runs(prop="l177")       all runs of a propeller's tables, as dicts

Sections: Clark Y with characteristics estimated from the thickness ratio, the same model
as the R&M 2357 and NACA Report 747 comparisons (validation/naca747/model.py), including its
critical Mach number Mcrit = 0.85 - 1.2 t (0.76 at the outer blade of the Curtiss 714,
t = 0.072). L-594 and L7D10 don't name their sections (only the design lift coefficients);
the same estimate is used for them, so the comparison tests the whole model. XROTOR's
compressibility model is the one under test: dCD = 10 (M - Mcrit_eff)^3, Prandtl-Glauert
lift slope, local Mach number capped at 0.99.

Not modeled: the cooling cuffs (root), the spinner and cowling, the measured thrust being
a slipstream total-pressure survey (the report's method, reference 1 there).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from pyxrotor.constants import PI
from pyxrotor.rotor import Rotor

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("naca747_model", HERE.parent / "naca747" / "model.py")
_m747 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m747)
clark_y = _m747.clark_y                              # the shared section model

PROPS = {
    "l177": ("geometry_l177.json", ("l177_table1.csv", "l177_table2.csv")),
    "l594": ("geometry_l594.json", ("l594_table1.csv", "l594_table2.csv", "l594_table3.csv")),
    "l7d10": ("geometry_l7d10.json", ("l7d10_table1.csv", "l7d10_table2.csv", "l7d10_table3.csv")),
}
FT = 0.3048
RHO_SL = 1.2250                   # sea-level density for sigma (kg/m^3)
GAMMA, R_AIR = 1.4, 287.05


def geometry(prop: str) -> dict:
    return json.loads((HERE / PROPS[prop][0]).read_text())


GEO = geometry("l177")            # the Curtiss 714 (kept for scripts that use it)


def sections(mcrit_shift: float = 0.0, geo: dict = GEO):
    """Clark Y sections at the blade-form thickness ratios; mcrit_shift is added to every
    section's critical Mach number (for the sensitivity study)."""
    r = np.array(geo["r_over_R"])
    t = np.array(geo["h_over_b"])
    xs = [r[0], 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 1.0]
    out = []
    for x in xs:
        s = clark_y(float(np.interp(x, r, t)), float(x))
        s.mcrit = s.mcrit + mcrit_shift
        out.append(s)
    return out


def blade(beta75: float, geo: dict = GEO):
    """ARBI stations, hub to tip: r/R, c/R, blade angle (deg). Beyond the last digitized
    point the chord and twist are extrapolated along their last slope (the chord not below
    0.02 R), and the chord is 0.02 R at the tip (a rounded tip)."""
    r = np.array(geo["r_over_R"])
    c = 2.0 * np.array(geo["b_over_D"])
    tw = np.array(geo["twist_deg_from_0.75R"])
    rr = np.concatenate([np.linspace(r[0], 0.90, 15), [0.95, 1.0]])
    cc = np.interp(rr, r, c)
    tt = np.interp(rr, r, tw)
    slope_c = (c[-1] - c[-2]) / (r[-1] - r[-2])
    slope_t = (tw[-1] - tw[-2]) / (r[-1] - r[-2])
    for k, x in enumerate(rr):
        if x > r[-1]:
            cc[k] = max(c[-1] + slope_c * (x - r[-1]), 0.02)
            tt[k] = tw[-1] + slope_t * (x - r[-1])
    cc[-1] = 0.02
    return rr, cc, beta75 + tt


def sutherland(t_k: float) -> float:
    return 1.458e-6 * t_k ** 1.5 / (t_k + 110.4)


def run(row: dict, mach_scale: float = 1.0, beta_guess: float | None = None,
        mcrit_shift: float = 0.0, prop: str | None = None, **fixes) -> dict:
    prop = prop or row.get("prop", "l177")
    geo = geometry(prop)
    D = geo["diameter_ft"] * FT
    R = 0.5 * D
    J, cp, n = float(row["J"]), float(row["CP"]), float(row["n_rps"])
    m = float(row["M"]) * mach_scale
    vel = J * n * D
    a_test = vel / float(row["M"])                     # the test's speed of sound
    a = vel / m
    rho = float(row["sigma"]) * RHO_SL
    rmu = sutherland(a_test * a_test / (GAMMA * R_AIR))  # viscosity of the test air
    power = cp * rho * n ** 3 * D ** 5
    beta0 = beta_guess
    if beta0 is None:
        beta0 = float(row["beta_deg"]) if row.get("beta_deg") else 30.0 + 7.0 * (J - 1.0)

    def solve(beta_start):
        r = Rotor(fix_all=True, **fixes)
        r.s.rho, r.s.vso, r.s.rmu = rho, a, rmu
        r.set_sections(sections(mcrit_shift, geo))
        rr, cc, bb = blade(beta_start, geo)
        rr[0] = (rr[0] * R) / R                        # ARBI's XI0 = (hub R)/R, to the ulp
        r.arbi(nblds=geo["nblds"], vel=vel, tip_radius=R, hub_radius=rr[0] * R,
               r_over_R=rr, c_over_R=cc, beta_deg=bb)
        return r.solve_power(power, fix="rpm", rpm=60.0 * n), r

    ok, r = solve(beta0)
    retried = False
    if not ok and beta_guess is None and mach_scale > 0.05:
        # start again from the blade angle of the same run without compressibility
        start = run(row, mach_scale=0.01, mcrit_shift=mcrit_shift, prop=prop, **fixes)
        if start["converged"]:
            ok, r = solve(start["beta75"])
            retried = True
    p = r.performance()
    beta75 = float(np.interp(0.75, r.s.xi[:r.s.ii], r.s.beta[:r.s.ii])) * 180.0 / PI
    ct = p["thrust"] / (rho * n * n * D ** 4)
    return {"prop": prop, "run": row["run"], "converged": bool(ok), "J": J, "CP": cp, "M": m,
            "Mt": m * math.sqrt(1.0 + (PI / J) ** 2), "eta": p["efficiency"], "CT": ct,
            "beta75": beta75, "retried": retried}


def table(name: str) -> list[dict]:
    with open(HERE / name) as f:
        return list(csv.DictReader(f))


def runs(prop: str = "l177") -> list[dict]:
    return [dict(row, table=t, prop=prop) for t in PROPS[prop][1] for row in table(t)]
