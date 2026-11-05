"""pyxrotor model of the four-blade propellers tested on a YP-47M in the NACA Cleveland
altitude wind tunnel (1946; "Altitude-wind-tunnel investigation of performance of several
propellers on YP-47M airplane at high blade loadings", parts I, II, IV, V):

    aeroproducts_H20C-162   NTRS 20030063991   12 ft 7 in
    curtiss_838-1C2-18R1    NTRS 20030063996   13 ft
    curtiss_732-1C2-0       NTRS 20030064146   13 ft
    curtiss_836-14C2-18R1   NTRS 20030063262   13 ft

All have NACA 16-series sections. Geometry (geometry.json) is traced from each report's
figure 1 (digitizing/). Sections are the 16-series estimate of the RM L51L28 comparison
(validation/naca_l51l28): Mcrit = 0.85 - 1.2 t plus a shift, at the traced design CL.

The tests give efficiency against power coefficient at constant J and free-stream Mach
number (0.3-0.5). A point is reproduced at the test's M0, J and CP. The air is ISA at
30,000 ft (the tests ran at density altitudes of 20,000-45,000 ft; only the Reynolds number
depends on it). No body: the YP-47M cowl's velocity field isn't known.

    run(prop, M0, J, cp, mcrit_shift=0, mach_scale=1, stall_margin=None)
        XROTOR POWE at fixed rpm, continued in Mach number from the blade angle without
        compressibility (mach_scale 0.01). If POWE fails at full Mach number, the blade angle
        is found from fixed-pitch solves on a 0.5 deg grid as in validation/nasa_propfan.
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
GEO = json.loads((HERE / "geometry.json").read_text())
_spec = importlib.util.spec_from_file_location("l51l28_model", HERE.parent / "naca_l51l28" / "model.py")
ml51 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ml51)

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
        s = ml51.section(float(np.interp(x, r, t)), float(x), cld, mcrit_shift)
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


def _state(r, n, D):
    p = r.performance()
    return {"CP": p["power"] / (RHO_ALT * n ** 3 * D ** 5), "CT": p["thrust"] / (RHO_ALT * n * n * D ** 4),
            "eta": p["efficiency"],
            "beta75": float(np.interp(0.75, r.s.xi[:r.s.ii], r.s.beta[:r.s.ii])) * 180.0 / PI}


def _powe(prop, M0, J, cp, beta, mcrit_shift, ms, sm=None):
    r, n, D = _rotor(prop, M0, J, beta, mcrit_shift, ms, sm)
    ok = r.solve_power(cp * RHO_ALT * n ** 3 * D ** 5, fix="rpm", rpm=60.0 * n)
    return bool(ok), _state(r, n, D)


def _fixed(prop, M0, J, beta, mcrit_shift, ms, sm=None):
    r, n, D = _rotor(prop, M0, J, beta, mcrit_shift, ms, sm)
    ok = r.solve_rpm(60.0 * n)
    return bool(ok), _state(r, n, D)


def run(prop: str, M0: float, J: float, cp: float, mcrit_shift: float = 0.0,
        mach_scale: float = 1.0, stall_margin: float | None = None) -> dict:
    out = {"prop": prop, "M0": M0, "J": J, "CP": cp, "method": None, "usable": False,
           "eta": float("nan"), "CT": float("nan"), "beta75": float("nan"),
           "Mt": M0 * mach_scale * math.sqrt(1.0 + (PI / J) ** 2)}
    # start from the advance angle at 0.75R plus 5 deg; POWE from far below it can settle
    # on a negative-thrust branch
    guess = math.degrees(math.atan(J / (0.75 * PI))) + 5.0
    ok, st = _powe(prop, M0, J, cp, guess, mcrit_shift, 0.01, stall_margin)
    if not ok:
        return out
    beta0 = beta = st["beta75"]
    if mach_scale > 0.01:
        for ms in [s for s in (0.6, 0.8, 0.9) if s < mach_scale] + [mach_scale]:
            ok, st = _powe(prop, M0, J, cp, beta, mcrit_shift, ms, stall_margin)
            if not ok:
                break
            beta = st["beta75"]
    if ok:
        out.update(st, CP=cp, method="POWE", usable=True)
        return out
    # fixed-pitch scan, linear interpolation at the target CP between converged neighbours
    pts = []
    for db in np.arange(-6.0, 4.01, 0.5):
        ok, st = _fixed(prop, M0, J, beta0 + db, mcrit_shift, mach_scale, stall_margin)
        if ok and math.isfinite(st["eta"]):
            pts.append(st)
    for lo, hi in zip(pts, pts[1:]):
        if lo["CP"] <= cp <= hi["CP"]:
            f = (cp - lo["CP"]) / (hi["CP"] - lo["CP"])
            for k in ("CT", "eta", "beta75"):
                out[k] = lo[k] + f * (hi[k] - lo[k])
            out["method"] = "scan"
            out["usable"] = hi["beta75"] - lo["beta75"] < 1.1 and abs(hi["eta"] - lo["eta"]) < 0.06
            break
    return out


def points() -> list[dict]:
    with open(HERE / "efficiency_vs_cp.csv") as f:
        return list(csv.DictReader(f))
