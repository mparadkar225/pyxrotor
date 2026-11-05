"""pyxrotor model of the Curtiss No. 714-1C2-12 four-blade propeller on the Republic P-47C
(NACA ACR L4L07 / WR L-177, Vogeley 1944), for checking XROTOR's compressibility model
against flight data.

Each test run (tables I and II, transcribed in l177_table*.csv) is reproduced at its own
J, power coefficient, rotational speed, density and airplane Mach number: XROTOR POWE at fixed
rpm sets the blade angle to absorb the measured power, as the propeller's governor did, and
the efficiency, thrust coefficient and blade angle are compared with the measured ones.

    run(row, mach_scale=1)  one test run; mach_scale scales the airplane (and tip) Mach
                            number by changing the speed of sound, not V or n, so J, CP and
                            Reynolds number stay the test's: mach_scale=0.01 is the model
                            without compressibility, mcrit_shift raises every section's
                            critical Mach number
    runs()                  all runs of tables I and II, as dicts

Sections: Clark Y with characteristics estimated from the thickness ratio, the same model
as the R&M 2357 and NACA Report 747 comparisons (validation/naca747/model.py), including its
critical Mach number Mcrit = 0.85 - 1.2 t (0.76 at the outer blade, t = 0.072). XROTOR's
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

GEO = json.loads((HERE / "geometry_l177.json").read_text())
FT = 0.3048
D = GEO["diameter_ft"] * FT
R = 0.5 * D
HUB = GEO["r_over_R"][0]
RHO_SL = 1.2250                   # sea-level density for sigma (kg/m^3)
GAMMA, R_AIR = 1.4, 287.05


def sections(mcrit_shift: float = 0.0):
    """Clark Y sections at the blade-form thickness ratios; mcrit_shift is added to every
    section's critical Mach number (for the sensitivity study)."""
    r = np.array(GEO["r_over_R"])
    t = np.array(GEO["h_over_b"])
    xs = [HUB, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 1.0]
    out = []
    for x in xs:
        s = clark_y(float(np.interp(x, r, t)), x)
        s.mcrit = s.mcrit + mcrit_shift
        out.append(s)
    return out


def blade(beta75: float):
    """ARBI stations, hub to tip: r/R, c/R, blade angle (deg). The chord falls from the
    last digitized point (0.915) to 0.02 R at the tip (a rounded tip); twist extrapolated."""
    r = np.array(GEO["r_over_R"])
    c = 2.0 * np.array(GEO["b_over_D"])
    tw = np.array(GEO["twist_deg_from_0.75R"])
    rr = np.concatenate([np.linspace(HUB, 0.90, 15), [0.95, 1.0]])
    cc = np.interp(rr, r, c)
    tt = np.interp(rr, r, tw)
    slope_c = (c[-1] - c[-2]) / (r[-1] - r[-2])
    slope_t = (tw[-1] - tw[-2]) / (r[-1] - r[-2])
    for k, x in enumerate(rr):
        if x > r[-1]:
            cc[k] = max(c[-1] + slope_c * (x - r[-1]), 0.02 * 2.0 / 2.0)
            tt[k] = tw[-1] + slope_t * (x - r[-1])
    cc[-1] = 0.02
    return rr, cc, beta75 + tt


def sutherland(t_k: float) -> float:
    return 1.458e-6 * t_k ** 1.5 / (t_k + 110.4)


def run(row: dict, mach_scale: float = 1.0, beta_guess: float | None = None,
        mcrit_shift: float = 0.0, **fixes) -> dict:
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
    r = Rotor(fix_all=True, **fixes)
    r.s.rho, r.s.vso, r.s.rmu = rho, a, rmu
    r.set_sections(sections(mcrit_shift))
    rr, cc, bb = blade(beta0)
    r.arbi(nblds=GEO["nblds"], vel=vel, tip_radius=R, hub_radius=HUB * R, r_over_R=rr,
           c_over_R=cc, beta_deg=bb)
    ok = r.solve_power(power, fix="rpm", rpm=60.0 * n)
    p = r.performance()
    beta75 = float(np.interp(0.75, r.s.xi[:r.s.ii], r.s.beta[:r.s.ii])) * 180.0 / PI
    ct = p["thrust"] / (rho * n * n * D ** 4)
    return {"run": row["run"], "converged": bool(ok), "J": J, "CP": cp, "M": m,
            "Mt": m * math.sqrt(1.0 + (PI / J) ** 2), "eta": p["efficiency"], "CT": ct,
            "beta75": beta75}


def table(name: str) -> list[dict]:
    with open(HERE / name) as f:
        return list(csv.DictReader(f))


def runs() -> list[dict]:
    return [dict(row, table=t) for t in ("l177_table1.csv", "l177_table2.csv") for row in table(t)]
