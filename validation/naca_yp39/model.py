"""pyxrotor model of the Hamilton Standard 3155-6 three-blade propeller flight-tested at high
speed on a Bell YP-39 (A. W. Vogeley, "Preliminary results of flight tests of a conventional
three-blade propeller at high speeds", NACA, 1942; NTRS 20150011182).

The blade is the one of NACA Report 747, whose traced geometry and Clark Y section estimate
are reused (validation/naca747). Here it is three blades on a 10 ft 3/4 in diameter, with
a spinner to 12 3/8 in radius (0.205R) and no cuffs.

Each run (flight_runs.csv) is reproduced at the test's true airspeed, density altitude
and power. Engine speed 3000 rpm with 1.8:1 gearing gives 1666.7 propeller rpm, and the
reported tip Mach numbers confirm it (1656-1671 rpm). The speed of sound comes from the
reported V and M; the density from the density altitude. XROTOR POWE at fixed rpm sets the
blade angle (the test's pitch was locked at the governed angle).

    run(row, mcrit_shift=0, mach_scale=1)   one run; mcrit_shift raises every section's
                                            critical Mach number, mach_scale=0.01 removes
                                            compressibility (speed of sound scaled up)
Not modeled: the round blade shanks between the spinner and about 0.35R. The test shows
their drag (negative total head inboard), which grows with speed. The model's
root sections are Clark Y estimates at the shank thickness (t/c 0.44), so its shank drag is far
smaller. Compare with the "working part" efficiency the report also gives (thrust of the
positive-thrust part of the disk over the full torque, so slightly low).
"""
from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

import numpy as np

from pyxrotor.constants import PI
from pyxrotor.rotor import Rotor

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("naca747_model", HERE.parent / "naca747" / "model.py")
m747 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m747)

INCH = 0.0254
D = (120.0 + 0.75) * INCH
R = 0.5 * D
HUB = 12.375 * INCH / R                 # spinner radius 12 3/8 in (figure 2)
RPM = 3000.0 / 1.8
GAMMA, R_AIR = 1.4, 287.05
HP = 745.7


def density_from_altitude(h_ft: float) -> float:
    """ISA troposphere density at density altitude h."""
    h = h_ft * 0.3048
    return 1.2250 * (1.0 - 2.25577e-5 * h) ** 4.2559


def sections(mcrit_shift: float = 0.0):
    out = []
    r = np.array(m747.GEO["r_over_R"])
    t = np.array(m747.GEO["h_over_b"])
    for x in [HUB, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 1.0]:
        s = m747.clark_y(float(np.interp(x, r, t)), float(x))
        s.mcrit = s.mcrit + mcrit_shift
        out.append(s)
    return out


def blade(beta75: float):
    g = m747.GEO
    r = np.array(g["r_over_R"])
    c = np.array(g["c_over_R"])
    tw = np.array(g["twist_deg_from_0.75R"])
    rr = np.concatenate([np.linspace(HUB, 0.95, 15), [1.0]])
    cc = np.interp(rr, r, c)
    tt = np.interp(rr, r, tw)
    cc[-1] = 0.055                       # rounded tip, as in the Report 747 model
    tt[-1] = tw[-1] + (tw[-1] - tw[-2])
    return rr, cc, beta75 + tt


def run(row: dict, mcrit_shift: float = 0.0, mach_scale: float = 1.0, beta_guess: float = 45.0) -> dict:
    V = float(row["V_mph"]) * 0.44704
    M = float(row["M"])
    a = V / M
    rho = density_from_altitude(float(row["density_altitude_ft"]))
    T = a * a / (GAMMA * R_AIR)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    n = RPM / 60.0
    power = float(row["bhp"]) * HP

    def solve(beta, ms):
        r = Rotor(fix_all=True)
        r.s.rho, r.s.vso, r.s.rmu = rho, a / ms, mu
        r.set_sections(sections(mcrit_shift))
        rr, cc, bb = blade(beta)
        rr[0] = (rr[0] * R) / R
        r.arbi(nblds=3, vel=V, tip_radius=R, hub_radius=rr[0] * R, r_over_R=rr, c_over_R=cc,
               beta_deg=bb)
        ok = r.solve_power(power, fix="rpm", rpm=RPM)
        b75 = float(np.interp(0.75, r.s.xi[:r.s.ii], r.s.beta[:r.s.ii])) * 180.0 / PI
        return ok, r, b75

    # continuation in Mach number from the incompressible blade angle
    ok, r, beta = solve(beta_guess, 0.01)
    for ms in [s for s in (0.6, 0.8, 0.9) if s < mach_scale] + [mach_scale]:
        if ms <= 0.01:
            break
        ok, r, b = solve(beta, ms)
        if not ok:
            break
        beta = b
    p = r.performance()
    J = V / (n * D)
    return {"V_mph": float(row["V_mph"]), "converged": bool(ok), "eta": p["efficiency"],
            "beta75": beta, "J": J, "CP": power / (rho * n ** 3 * D ** 5),
            "CT": p["thrust"] / (rho * n * n * D ** 4),
            "Mt": M * mach_scale * math.sqrt(1.0 + (PI / J) ** 2), "rotor": r}


def runs() -> list[dict]:
    with open(HERE / "flight_runs.csv") as f:
        return list(csv.DictReader(f))
