"""Propeller sizing study for a 500 mph fighter (not a validation; a design exploration).

Requirement (from the user): maximum speed about 500 mph at 25,000-35,000 ft, 2900 hp at
take-off / WEP, high-speed cruise about 350 mph. Layouts: 4- and 5-blade single rotation,
6- (3+3) and 8-blade (4+4) contra-rotation.

Assumptions (stated in README.md):
  * max-speed point 500 mph at 30,000 ft with the full 2900 hp (a turbo or two-stage engine
    holding WEP to its critical altitude); checked at 25,000 and 35,000 ft too;
  * cruise 350 mph at 25,000 ft, 1300 shp, engine and propeller rpm at 80% of maximum;
  * take-off: sea level, 2900 hp at maximum rpm, static and at 130 mph;
  * NACA 16-series sections, design CL as given, t/c 0.30 at the 0.2R root (faired, no round
    shank) to 0.055 at the tip; calibrated Mcrit (pyxrotor.sections), 7.55 fixes on;
  * spinner 0.22 R; contra rows 0.12 D apart, equal torques, same rpm.

Model accuracy (validation/): against wartime tests the calibrated model's efficiency is
0.03-0.06 high at low speed and follows the measured compressibility losses to helical tip
Mach ~1.15, so compare layouts with each other rather than trusting absolute values.
XROTOR's static (zero-speed) power is unreliable with stalled inner blades (RM L51L28), so
take-off is judged at 130 mph. Contra pair solves near the local-Mach cap sometimes fail
(recorded as not ok).

Each layout is designed with DESI (minimum induced loss) at the max-speed point for its
share of the power; contra rows are redesigned, each in the other's induced field
(pyxrotor.field), three times, then analysed together with RotorPair.

    PYTHONPATH=src python studies/fighter_500mph/study.py   -> results.json
"""
from __future__ import annotations

import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from pyxrotor import Rotor, robust, sections
from pyxrotor.pair import RotorPair

HERE = Path(__file__).resolve().parent
HP = 745.7
MPH = 0.44704
FT = 0.3048
P_MAX = 2900 * HP
V_MAX, ALT_MAX = 500 * MPH, 30000 * FT / 1000
V_CRU, ALT_CRU, P_CRU, RPM_CRU = 350 * MPH, 25000 * FT / 1000, 1300 * HP, 0.80
V_CLIMB = 130 * MPH
HUB = 0.22
SPACING = 0.12                                   # contra row spacing / D
THICKNESS = [(0.20, 0.30), (0.30, 0.18), (0.45, 0.11), (0.60, 0.085), (0.75, 0.07),
             (0.90, 0.06), (1.00, 0.055)]


def blade_sections(cl):
    return [sections.naca16(t, x, cl) for x, t in THICKNESS]


def rpm_for(d, mt):
    """rpm giving helical tip Mach mt at the max-speed point."""
    r = Rotor()
    r.atmosphere(ALT_MAX)
    a = r.s.vso
    u = math.sqrt((mt * a) ** 2 - V_MAX ** 2)
    return u / (math.pi * d) * 60.0


def new_rotor(cl):
    r = Rotor(fix_all=True)
    r.atmosphere(ALT_MAX)
    r.set_sections(blade_sections(cl))
    return r


def design(r, b, d, rpm, power, cl):
    """DESI at the max-speed point; designed without compressibility if DESI fails with it
    (it does near tip Mach 1), the sections' own Mcrit restored for analysis."""
    kw = dict(nblds=b, tip_radius=d / 2, hub_radius=HUB * d / 2, vel=V_MAX, rpm=rpm,
              power=power, cl=cl)
    try:
        if r.design(**kw) and np.all(np.isfinite(r.stations()["beta"])):
            return "compressible"
    except Exception:
        pass
    vso = r.s.vso
    r.s.vso = vso * 100.0
    try:
        ok = r.design(**kw)
    finally:
        r.s.vso = vso
    if not ok:
        raise RuntimeError("design failed")
    return "incompressible"


def geometry(r):
    st = r.stations()
    x, c = st["xi"], st["ch"]                         # c/R
    af = 1e5 / 16 * float(np.trapezoid(c / 2 * x ** 3, x))
    i75 = int(np.argmin(np.abs(x - 0.75)))
    return {"AF_per_blade": round(af, 1),
            "chord_0.75R_m": round(float(np.interp(0.75, x, c)) * r.s.rad, 3),
            "chord_max_m": round(float(c.max()) * r.s.rad, 3),
            "x": [round(float(v), 3) for v in x], "c_over_R": [round(float(v), 4) for v in c],
            "beta_deg": [round(math.degrees(float(v)), 2) for v in st["beta"]],
            "beta75_design_deg": round(math.degrees(float(st["beta"][i75])), 2)}


def _sol(sol):
    return {"ok": bool(sol.usable), "eta": _r(sol.efficiency), "thrust_N": _r(sol.thrust, 0),
            "beta75_deg": _r(sol.beta75_deg, 2), "method": sol.method}


def _r(v, n=4):
    return None if v is None or not math.isfinite(v) else round(float(v), n)


def single(b, d, mt, cl):
    rpm = rpm_for(d, mt)
    r = new_rotor(cl)
    how = design(r, b, d, rpm, P_MAX, cl)
    out = {"layout": f"{b}-blade", "D_m": d, "D_ft": round(d / FT, 2), "Mt_max": mt,
           "rpm": round(rpm, 1), "cl": cl, "design": how, "geometry": geometry(r)}
    for name, alt in (("max_25k", 25000), ("max", 30000), ("max_35k", 35000)):
        r.atmosphere(alt * FT / 1000)
        r.velocity(V_MAX)
        out[name] = _sol(robust.solve_power(r, P_MAX, rpm))
    r.atmosphere(ALT_CRU)
    r.velocity(V_CRU)
    out["cruise"] = _sol(robust.solve_power(r, P_CRU, RPM_CRU * rpm))
    r.atmosphere(0.0)
    for name, v in (("static", 0.002 * rpm / 60 * d), ("climb", V_CLIMB)):
        r.velocity(v)
        out[name] = _sol(robust.solve_power(r, P_MAX, rpm))
    out["static"]["lb_per_hp"] = _r(out["static"]["thrust_N"] / 4.448 / 2900, 3) \
        if out["static"]["thrust_N"] else None
    return out


def _pair_sol(pair, power, rpm):
    try:
        res = pair.solve_power(power, rpm)
    except Exception as e:
        return {"ok": False, "error": type(e).__name__}
    bf = math.degrees(float(np.interp(0.75, pair.front.s.xi[:pair.front.s.ii],
                                      pair.front.s.beta[:pair.front.s.ii])))
    br = math.degrees(float(np.interp(0.75, pair.rear.s.xi[:pair.rear.s.ii],
                                      pair.rear.s.beta[:pair.rear.s.ii])))
    return {"ok": bool(res.converged), "eta": _r(res.efficiency), "thrust_N": _r(res.thrust, 0),
            "front_thrust_share": _r(res.front["thrust"] / res.thrust, 3),
            "beta75_front_rear_deg": [_r(bf, 2), _r(br, 2)],
            "swirl_loss": _r(res.swirl_power_fraction)}


def contra(b, d, mt, cl, passes=3):
    rpm = rpm_for(d, mt)
    front, rear = new_rotor(cl), new_rotor(cl)
    hows = [design(front, b, d, rpm, P_MAX / 2, cl), design(rear, b, d, rpm, P_MAX / 2, cl)]
    pair = RotorPair(front, rear, spacing=SPACING * d)
    vso = front.s.vso
    for rotor in (front, rear):                    # the passes without compressibility:
        rotor.s.vso = vso * 100.0                  # the field is set by the loading
    for _ in range(passes):                        # each row redesigned in the other's field
        for which, rotor, other in (("rear", rear, front), ("front", front, rear)):
            if not other.solve_power(P_MAX / 2, fix="pitch"):
                raise RuntimeError(f"{which} pass: the other row did not converge")
            r, u, v = pair._imposed(which, other.induced_field())
            pair._impose(rotor, r, u, v)
            if not rotor.design(nblds=b, tip_radius=d / 2, hub_radius=HUB * d / 2, vel=V_MAX,
                                rpm=rpm, power=P_MAX / 2, cl=cl):
                raise RuntimeError(f"{which} redesign failed")
    for rotor in (front, rear):
        rotor.s.vso = vso
    out = {"layout": f"{2 * b}-blade contra ({b}+{b})", "D_m": d, "D_ft": round(d / FT, 2),
           "Mt_max": mt, "rpm": round(rpm, 1), "cl": cl,
           "design": "incompressible (rows redesigned in each other's field)",
           "geometry_front": geometry(front), "geometry_rear": geometry(rear)}
    for name, alt, v, p, n in (("max_25k", 25000, V_MAX, P_MAX, rpm),
                               ("max", 30000, V_MAX, P_MAX, rpm),
                               ("max_35k", 35000, V_MAX, P_MAX, rpm),
                               ("cruise", 25000, V_CRU, P_CRU, RPM_CRU * rpm),
                               ("climb", 0, V_CLIMB, P_MAX, rpm),
                               ("static", 0, 0.002 * rpm / 60 * d, P_MAX, rpm)):
        pair.atmosphere(alt * FT / 1000)
        pair.velocity(v)
        out[name] = _pair_sol(pair, p, n)
    if out["static"].get("thrust_N"):
        out["static"]["lb_per_hp"] = _r(out["static"]["thrust_N"] / 4.448 / 2900, 3)
    return out


def job(args):
    kind, b, d, mt, cl = args
    try:
        return (single if kind == "single" else contra)(b, d, mt, cl)
    except Exception as e:
        return {"layout": f"{kind} {b}", "D_m": d, "Mt_max": mt, "cl": cl, "error": repr(e)}


GRID = [(k, b, d, mt, cl)
        for k, b in (("single", 4), ("single", 5), ("contra", 3), ("contra", 4))
        for d in (3.40, 3.70, 4.00, 4.30)
        for mt in (0.95, 1.00, 1.05, 1.10, 1.15)
        for cl in (0.4, 0.5)]

if __name__ == "__main__":
    grid = GRID if len(sys.argv) < 2 else [g for g in GRID if g[0] + str(g[1]) in sys.argv[1:]]
    with ProcessPoolExecutor() as ex:
        res = list(ex.map(job, grid))
    (HERE / "results.json").write_text(json.dumps(res, indent=1))
    for r in res:
        if "error" in r:
            print(r["layout"], r["D_m"], r["Mt_max"], r["cl"], r["error"])
            continue
        e = {k: (r[k].get("eta") if r[k].get("ok") else None) for k in ("max", "cruise", "climb")}
        print(f"{r['layout']:22s} D {r['D_ft']:5.2f} ft Mt {r['Mt_max']:.2f} rpm {r['rpm']:6.0f} "
              f"cl {r['cl']}  eta max {e['max']} cruise {e['cruise']} climb {e['climb']}")
