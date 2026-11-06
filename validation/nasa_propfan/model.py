"""pyxrotor models of NASA's SR-2 (straight) and SR-3 (45 deg tip sweep) eight-blade propfan
models, 62.23 cm reference diameter, tested at Mach 0.45-0.85 in the NASA Lewis 8x6 ft wind
tunnel (NASA TM-79124, Jeracki, Mikkelson & Blaha 1979; NASA CR-3505, Rohrbach et al. 1982;
SR-2 also in NASA TM X-73612, Mikkelson et al. 1977).

XROTOR is a straight lifting line with 2-D sections: it has no sweep. The question is how far
that gets for propfans, and what a sweep correction would have to do.

    run(prop, M, J, cp=None, beta75=None, ...)
        one operating point at free-stream Mach number M and advance ratio J (reference
        diameter). With cp: the blade angle is found to absorb that power coefficient (XROTOR
        POWE at fixed rpm). With beta75: the blade angle is fixed (RPM). Returns CP, CT, net
        efficiency and the blade angle at 0.75R.

Model choices (each is an option to run with and without):
  * Geometry from the reports (geometry.json): b/D, t/b, twist, design CL against r/R, from
    0.24R (the spinner) to the tip. Activity factor 202 and 234 from the traces against 203
    and 235 published; integrated design CL 0.075 and 0.212 against 0.08 and 0.214.
  * Sections: the NACA 16-series estimate of the RM L51L28 comparison, applied to the
    65/circular-arc root sections as well: zero-lift angle -cld/2pi, minimum drag at
    CL = cld, CDmin 0.0065 + 0.25 t^2, CLmax = cld + 0.85 (1.1 + cld where t > 0.12),
    Mcrit = 0.85 - 1.2 t (0.826 at the 2% tip) plus mcrit_shift.
  * sweep="mcrit": the simple-sweep (cosine) rule on the drag rise only: each section's
    Mcrit is divided by cos(aerodynamic sweep). sweep="cosine" also multiplies the
    incompressible lift slope by cos(sweep), as simple sweep theory does (Hamilton
    Standard's strip method applied the cosine rule to the airfoil data and the local
    velocity). XROTOR's Prandtl-Glauert factor still uses the full Mach number, and the
    zero-lift angle and CLmax are left as they are.
  * sweep="core": the compiled core's swept-blade model (pyxrotor.sweep): the quarter-chord
    line from the aerodynamic sweep (s' = tan sweep), each station's effective sweep from the
    local flow at every iteration, sections evaluated in the normal plane (normal Mach number
    and angle, CL cos^2, pressure drag cos^3). The airfoil thickness is taken as the normal
    section's, t/b / cos(sweep).
  * nacelle=True: XROTOR's body velocity UBODY = M_local/M0 - 1 from the potential-flow
    axial Mach number about the 0.35 D nacelle at M0 = 0.8 (TM X-73612 fig. 11), the same
    fraction at every Mach number.
  * SR-3's tip radius grows with blade angle (swept blade pitched about its axis): the true
    diameter is 1.019 Dref at 60 deg (CR-3505 fig. 14, design rpm). The blade is solved at the
    true radius and the coefficients reported on Dref, as in the reports.

Air: the tunnel draws atmospheric air (open circuit), so total conditions are sea level;
static temperature and density follow from M. Only the Reynolds number depends on them.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from pyxrotor.constants import PI
from pyxrotor.sweep import Sweep
from pyxrotor.sections import naca16
from pyxrotor.rotor import Rotor

HERE = Path(__file__).resolve().parent
GEO = json.loads((HERE / "geometry.json").read_text())
section16 = naca16                            # the 16-series estimate (t, xisect, cld, shift)

DREF = 0.6223
GAMMA, R_AIR = 1.4, 287.05
T0, P0 = 288.15, 101325.0


def true_diameter(prop: str, beta75: float) -> float:
    if prop != "SR-3":
        return DREF
    t = GEO["SR-3_true_diameter"]
    return DREF * float(np.interp(beta75, t["beta75_deg"], t["D_over_Dref"]))


def sections(prop: str, mcrit_shift: float = 0.0, sweep: str | None = None):
    g = GEO[prop]
    r = np.array(g["r_over_R"])
    out = []
    for x in [r[0], 0.28, 0.32, 0.36, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
              0.85, 0.90, 0.95, 1.0]:
        t = float(np.interp(x, r, g["h_over_b"]))
        cld = float(np.interp(x, r, g["design_cl"]))
        if sweep == "core":                      # the airfoil normal to the swept line
            t = t / math.cos(math.radians(float(np.interp(x, r, g["sweep_deg"]))))
        s = section16(t, float(x), cld, mcrit_shift)
        if sweep in ("mcrit", "cosine"):
            lam = math.radians(float(np.interp(x, r, g["sweep_deg"])))
            s.mcrit = s.mcrit / math.cos(lam)            # drag rise on the normal Mach number
            if sweep == "cosine":
                s.dclda = s.dclda * math.cos(lam)        # simple-sweep lift slope
        out.append(s)
    return out


_SWEEP = {}


def blade_sweep(prop: str):
    """The swept quarter-chord line (pyxrotor.sweep) for sweep="core": offsets along the
    chord, crossing the pitch axis where the sweep changes sign (0.39R on SR-3). The
    planform angle is adjusted (3 passes) until the core's effective sweep at the design
    point (Mach 0.8, J 3.06, CP 1.695, solved without compressibility) is the published
    aerodynamic sweep: the arc of the tangential offset and the twist make the two differ
    by up to 5 deg."""
    if prop in _SWEEP:
        return _SWEEP[prop]
    g = GEO[prop]
    xi = np.array(g["r_over_R"], dtype=float)
    target = np.array(g["sweep_deg"], dtype=float)
    plan = target.copy()
    sw = Sweep.from_angle(xi, plan)
    for _ in range(3):
        _SWEEP[prop] = sw
        o = _solve(prop, 0.8, 3.06, 1.695, 58.0, 0.0, "core", True, 0.01)
        r = o["rotor"]
        eff = np.interp(xi, r.s.xi[:r.s.ii], r.effective_sweep_deg()
                        * np.sign(np.interp(r.s.xi[:r.s.ii], xi, target)))
        plan = plan + (target - eff)
        sw = Sweep.from_angle(xi, plan)
    _SWEEP[prop] = sw
    return sw


def blade(prop: str, beta75: float):
    g = GEO[prop]
    r = np.array(g["r_over_R"])
    rr = np.concatenate([np.linspace(r[0], 0.95, 18), [0.98, 1.0]])
    cc = 2.0 * np.interp(rr, r, g["b_over_D"])
    tt = np.interp(rr, r, g["twist_deg_from_0.75R"])
    return rr, cc, beta75 + tt


def air(M: float):
    ts = T0 / (1.0 + 0.2 * M * M)
    ps = P0 / (1.0 + 0.2 * M * M) ** 3.5
    a = math.sqrt(GAMMA * R_AIR * ts)
    rho = ps / (R_AIR * ts)
    mu = 1.458e-6 * ts ** 1.5 / (ts + 110.4)
    return a, rho, mu


def _solve(prop, M, J, cp, beta, mcrit_shift, sweep, nacelle, mach_scale, nitera=40):
    """One XROTOR solve at blade angle beta (fixed pitch) or from beta (POWE)."""
    a, rho, mu = air(M)
    vel = M * a
    n = vel / (J * DREF)
    g = GEO[prop]
    D = true_diameter(prop, beta)
    R = 0.5 * D
    r = Rotor(fix_all=True)
    r.s.rho, r.s.vso, r.s.rmu = rho, a / mach_scale, mu
    r.set_sections(sections(prop, mcrit_shift, sweep))
    rr, cc, bb = blade(prop, beta)
    rr[0] = (rr[0] * R) / R                         # ARBI's XI0 = (hub R)/R, to the ulp
    r.arbi(nblds=g["nblds"], vel=vel, tip_radius=R, hub_radius=rr[0] * R,
           r_over_R=rr, c_over_R=cc, beta_deg=bb)
    if sweep == "core":
        r.set_sweep(blade_sweep(prop))
    if nacelle:
        nac = GEO["nacelle_local_axial_mach_M0_0.8"]
        r.s.ubody[:r.s.ii] = np.interp(r.s.xi[:r.s.ii], nac["r_over_R"], nac["M_local"]) / 0.8 - 1.0
    r.s.nitera = nitera
    solve = (lambda: r.solve_rpm(60.0 * n)) if cp is None else \
        (lambda: r.solve_power(cp * rho * n ** 3 * DREF ** 5, fix="rpm", rpm=60.0 * n))

    def state():
        p = r.performance()
        return np.array([p["power"] / (rho * n ** 3 * DREF ** 5), p["thrust"] / (rho * n * n * DREF ** 4),
                         p["efficiency"],
                         float(np.interp(0.75, r.s.xi[:r.s.ii], r.s.beta[:r.s.ii])) * 180.0 / PI])
    ok = solve()
    cycle = None
    if not ok:
        # Near the local-Mach cap (CLCDCM limits M^2 to 0.99 with zero derivative) the Newton
        # iteration can settle into a small limit cycle instead of converging. Continue one
        # iteration at a time and accept the mean over 8 iterations if the cycle is small
        # (efficiency within 0.005, CP within 2%, blade angle within 0.1 deg).
        r.s.loprini = False
        r.s.nitera = 1
        hist = []
        for _ in range(8):
            solve()
            hist.append(state())
        hist = np.array(hist)
        spread = hist.max(axis=0) - hist.min(axis=0)
        cycle = {"dCP_rel": float(spread[0] / abs(hist[:, 0].mean())), "deta": float(spread[2]),
                 "dbeta": float(spread[3])}
        small = bool(np.all(np.isfinite(hist)) and cycle["dCP_rel"] < 0.02
                     and cycle["deta"] < 0.005 and cycle["dbeta"] < 0.1)
        vals = hist.mean(axis=0) if small else hist[-1]
    else:
        vals, small = state(), False
    return {"prop": prop, "M": M, "J": J, "converged": bool(ok), "limit_cycle": cycle,
            "usable": bool(ok) or small,
            "CP": float(vals[0]), "CT": float(vals[1]), "eta": float(vals[2]), "beta75": float(vals[3]),
            "D_over_Dref": D / DREF,
            "tip_rel_mach": M * math.sqrt(1.0 + (PI * D / (J * DREF)) ** 2), "rotor": r}


def run(prop: str, M: float, J: float, cp: float | None = None, beta75: float | None = None,
        mcrit_shift: float = 0.0, sweep: str | None = None, nacelle: bool = True,
        mach_scale: float = 1.0, scan: tuple = (-7.0, 5.0, 0.5)) -> dict:
    """One point; J and cp on the reference diameter. mach_scale scales the Mach number by
    changing the speed of sound (V, n, J kept): mach_scale = 0.01 is no compressibility.

    With beta75: one fixed-pitch solve. With cp: XROTOR's own power solve (POWE) often fails
    near the local-Mach cap, so the blade angle is found outside it instead: fixed-pitch
    solves on a blade-angle grid (scan, relative to the blade angle without compressibility),
    and CP, CT, efficiency and blade angle interpolated linearly at the target CP between the
    two usable solutions that bracket it. "bracket" is their spacing in degrees."""
    args = (mcrit_shift, sweep, nacelle)
    if cp is None:
        return _solve(prop, M, J, None, beta75, *args, mach_scale)
    # the blade angle without compressibility: XROTOR's POWE converges there
    ref = _solve(prop, M, J, cp, 58.0, mcrit_shift, sweep, nacelle, 0.01)
    if mach_scale <= 0.05:
        ref["bracket"] = 0.0
        return ref
    b0 = ref["beta75"]
    pts = []
    for db in np.arange(scan[0], scan[1] + 1e-9, scan[2]):
        o = _solve(prop, M, J, None, b0 + db, *args, mach_scale)
        if o["usable"]:
            pts.append(o)
    out = {"prop": prop, "M": M, "J": J, "usable": False, "converged": False, "bracket": None,
           "CP": float("nan"), "CT": float("nan"), "eta": float("nan"), "beta75": float("nan"),
           "points": len(pts)}
    for lo, hi in zip(pts, pts[1:]):
        if lo["CP"] <= cp <= hi["CP"]:
            f = (cp - lo["CP"]) / (hi["CP"] - lo["CP"])
            for k in ("CT", "eta", "beta75", "tip_rel_mach", "D_over_Dref"):
                out[k] = lo[k] + f * (hi[k] - lo[k])
            out["CP"] = cp
            out["bracket"] = hi["beta75"] - lo["beta75"]
            # at most one failed solve between, and no jump in efficiency (a limit-cycle mean
            # off the solution branch)
            out["usable"] = out["bracket"] < 2.1 and abs(hi["eta"] - lo["eta"]) < 0.06
            out["converged"] = lo["converged"] and hi["converged"]
            break
    return out


def table(name: str) -> list[dict]:
    with open(HERE / name) as f:
        return list(csv.DictReader(f))
