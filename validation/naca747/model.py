"""pyxrotor models of the NACA Report 747 propellers (Hamilton Standard 3155-6 / 3156-6).

    single(beta, J)            four-blade single rotation: two 2-blade rows in tandem,
                               co-rotating, 9.7 in apart (as tested: two 2-way hubs)
    dual(beta, J)              four-blade dual rotation: front at beta, rear at the
                               report's rear setting, counter-rotating, 9.7 in apart
    disk4(beta, J)             reference: one 4-blade disk (what "4-blade" would be ideally)
Each returns a dict with CT, CP, eta and the per-rotor CP (front, rear) where there are two.

Sections: Clark Y, with characteristics estimated from the thickness ratio (the same model
as the R&M 2357 comparison: it reproduces those sections at their thicknesses); this is the
largest uncertainty, which the dual-minus-single comparison largely cancels.
Operating points follow the test: sea level, n at most 550 rpm and V at most 110 mph.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pyxrotor.aero import AeroSection
from pyxrotor.constants import PI
from pyxrotor.pair import RotorPair
from pyxrotor.rotor import Rotor

HERE = Path(__file__).resolve().parent
GEO = json.loads((HERE / "geometry.json").read_text())
FT, INCH = 0.3048, 0.0254
D = GEO["diameter_ft"] * FT
R = 0.5 * D
HUB = 0.5 * GEO["spinner_diameter_in"] * INCH / R          # 0.267: blades start at the spinner
SPACING = GEO["dual_spacing_in"] * INCH
N_MAX = 550.0 / 60.0                                         # rps
V_MAX = 110.0 * 0.44704                                      # m/s
REAR = {20: 19.8, 25: 24.6, 30: 29.4, 35: 34.2, 40: 38.9, 45: 43.7, 50: 48.4, 55: 53.1, 60: 57.6}


def clark_y(t: float, xisect: float) -> AeroSection:
    """Clark Y characteristics for thickness ratio t (estimated; see the module docstring)."""
    a0 = max(-0.37 * 100.0 * t, -6.0)
    return AeroSection(
        xisect=xisect, a0=a0 * PI / 180.0, dclda=max(6.3 - 4.0 * (t - 0.066), 5.2),
        clmax=1.3 if t < 0.20 else 1.1, clmin=-0.5, dclda_stall=0.1, dcl_stall=0.1,
        cdmin=0.0065 + 0.25 * t * t, cldmin=0.35, dcdcl2=0.007 + 0.01 * max(t - 0.066, 0.0),
        cmcon=-0.1, mcrit=max(0.85 - 1.2 * t, 0.5), reref=2.0e6, rexp=-0.15)


def sections():
    r = np.array(GEO["r_over_R"])
    t = np.array(GEO["h_over_b"])
    xs = [HUB, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 1.0]
    return [clark_y(float(np.interp(x, r, t)), x if x < 1.0 else 1.0) for x in xs]


def blade(beta75: float):
    """ARBI stations from the spinner to the tip: r/R, c/R, blade angle (deg)."""
    r = np.array(GEO["r_over_R"])
    c = np.array(GEO["c_over_R"])
    tw = np.array(GEO["twist_deg_from_0.75R"])
    rr = np.concatenate([np.linspace(HUB, 0.95, 14), [1.0]])
    cc = np.interp(rr, r, c)
    tt = np.interp(rr, r, tw)
    cc[-1] = 0.055                       # rounded tip: chord falling as from 0.90-0.95
    tt[-1] = tw[-1] + (tw[-1] - tw[-2])  # twist extrapolated
    return rr, cc, beta75 + tt


def operating_point(J: float):
    """(V m/s, rpm) as in the test: n = 550 rpm unless V would exceed 110 mph."""
    n = N_MAX
    v = J * n * D
    if v > V_MAX:
        v = V_MAX
        n = v / (J * D)
    return v, 60.0 * n


def rotor(beta75: float, vel: float, nblds: int = 2, **fixes) -> Rotor:
    r = Rotor(fix_all=True, **fixes)
    r.atmosphere(0.0)
    r.set_sections(sections())
    rr, cc, bb = blade(beta75)
    r.arbi(nblds=nblds, vel=vel, tip_radius=R, hub_radius=HUB * R, r_over_R=rr, c_over_R=cc,
           beta_deg=bb)
    return r


def _coeffs(thrust, power, rpm, rho):
    n = rpm / 60.0
    ct = thrust / (rho * n * n * D ** 4)
    cp = power / (rho * n ** 3 * D ** 5)
    return ct, cp


def disk4(beta: float, J: float) -> dict:
    v, rpm = operating_point(J)
    r = rotor(beta, v, nblds=4)
    ok = r.solve_rpm(rpm)
    p = r.performance()
    ct, cp = _coeffs(p["thrust"], p["power"], rpm, r.s.rho)
    return {"converged": ok, "CT": ct, "CP": cp, "eta": ct * J / cp if cp else float("nan")}


def _pair(beta_f, beta_r, J, counter):
    v, rpm = operating_point(J)
    pair = RotorPair(rotor(beta_f, v), rotor(beta_r, v), SPACING, counter_rotating=counter)
    res = pair.solve_rpm(rpm)
    rho = pair.front.s.rho
    ct, cp = _coeffs(res.thrust, res.power, rpm, rho)
    _, cpf = _coeffs(res.front["thrust"], res.front["power"], rpm, rho)
    _, cpr = _coeffs(res.rear["thrust"], res.rear["power"], rpm, rho)
    return {"converged": res.converged, "CT": ct, "CP": cp, "eta": ct * J / cp if cp else float("nan"),
            "CP_front": cpf, "CP_rear": cpr, "iterations": res.iterations}


def single(beta: float, J: float) -> dict:
    """Two 2-blade rows co-rotating at the same blade angle (the tested 'single' 4-blade)."""
    return _pair(beta, beta, J, counter=False)


def dual(beta: float, J: float) -> dict:
    """Front at beta, rear at the report's setting (figure 5), counter-rotating."""
    return _pair(beta, REAR[int(round(beta))], J, counter=True)
