"""Solves that keep working near XROTOR's local-Mach cap.

Not a port: XROTOR's OPER commands are in rotor.Rotor. These functions wrap them for
transonic tips, where the plain commands often fail:

* CLCDCM caps each section's local M^2 at 0.99 with a zero derivative. Near the cap the
  Newton iteration (APITER) can settle into a small limit cycle instead of converging, and
  the power solve at fixed rpm (POWE, APER(3,1)) fails more often than the fixed-pitch solve
  (RPM, APER(4,2)).
* That was found on the P-47, YP-47M and propfan comparisons (validation/), each of which
  had its own work-around. These functions are those work-arounds in one place.

    solve_fixed_pitch(rotor, rpm)        RPM, accepting a small limit cycle
    solve_power(rotor, power, rpm)       blade angle for a power at fixed rpm (constant speed)
    scan_power(rotor, power, rpm, beta)  solve_power's last resort on its own

Both return a Solution. When Solution.exact is True the rotor holds a converged XROTOR
solution of the request. Otherwise the values are a limit-cycle mean or interpolated, and
the rotor holds the nearest XROTOR state used (see each function).

solve_power tries, in order:
1. "powe": XROTOR's POWE at fixed rpm, from the rotor's current state.
2. "continuation": POWE without compressibility (speed of sound scaled up 100x), then at a
   rising fraction of the Mach number (0.6, 0.8, 0.9 by default) and at the real one, each
   step from the previous solution.
3. "scan": fixed-pitch solves on a blade-angle grid around the solution without
   compressibility (6 deg below to 4 deg above, every 0.5 deg, by default). Every value is
   interpolated linearly at the target power between the two usable solves that bracket it.
   The point is usable only if those two are at most max_bracket_deg apart (one or two
   failed solves between them) and their efficiencies differ by at most max_eta_jump (a
   limit-cycle mean off the solution branch shows up as a jump).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .constants import PI

CONTINUATION = (0.01, 0.6, 0.8, 0.9)
SCAN = (-6.0, 4.0, 0.5)


@dataclass(frozen=True)
class CycleTolerance:
    """A limit cycle is accepted when, over `iterations` single Newton iterations, the
    efficiency varies by less than deta, the power and thrust by less than dpower_rel and
    dthrust_rel (relative to their means) and the blade angle at 0.75R by less than
    dbeta_deg. (The thrust test matters near zero thrust and at static conditions, where
    the efficiency says nothing.)

    The iterates must also still be moving like Newton steps near a solution: APITER's
    relaxation factor RLX at least min_rlx and the RMS Newton step in GAM (APITER's RMS,
    from DGAMOLD) at most max_step_rms. A stalled iteration (RLX 0, a large step that is
    never taken) has no spread at all but isn't a solution: on a static 16-series case it
    gave RMS 0.19, against at most 6e-3 in the limit cycles of SR-3 at Mach 0.8-0.85.
    (APITER converges at RMS 1e-7.)"""
    iterations: int = 8
    deta: float = 0.005
    dpower_rel: float = 0.02
    dthrust_rel: float = 0.02
    dbeta_deg: float = 0.1
    min_rlx: float = 0.05
    max_step_rms: float = 0.02


@dataclass
class Solution:
    usable: bool
    exact: bool                      # the rotor holds a converged solution of the request
    method: str | None               # "rpm", "cycle" (fixed pitch); "powe", "continuation", "scan"
    power: float = math.nan          # W
    thrust: float = math.nan         # N
    torque: float = math.nan         # N m
    efficiency: float = math.nan
    CT: float = math.nan
    CP: float = math.nan
    beta75_deg: float = math.nan     # blade angle at 0.75R
    rpm: float = math.nan
    limit_cycle: dict | None = None  # spread of the accepted cycle (or of the last tried)
    bracket_deg: float | None = None  # scan: spacing of the two bracketing solves
    solves: int = 0                  # XROTOR solve calls made
    points: list = field(default_factory=list, repr=False)   # scan: the usable solves

    def as_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if k != "points"}


_KEYS = ("power", "thrust", "torque", "efficiency", "CT", "CP", "beta75_deg", "rpm")


def beta75_deg(rotor) -> float:
    s = rotor.s
    return float(np.interp(0.75, s.xi[:s.ii], s.beta[:s.ii])) * 180.0 / PI


def _values(rotor) -> dict:
    p = rotor.performance()
    return {"power": p["power"], "thrust": p["thrust"], "torque": p["torque"],
            "efficiency": p["efficiency"], "CT": p["CT"], "CP": p["CP"],
            "beta75_deg": beta75_deg(rotor), "rpm": p["rpm"]}


def _finite(v: dict) -> bool:
    return all(math.isfinite(v[k]) for k in _KEYS)


def _cycle(rotor, step, tol: CycleTolerance):
    """Run `step` (one solve call) tol.iterations times at one Newton iteration each, without
    re-initialization; return (accepted, mean values, spread). Restores NITERA and LOPRINI."""
    s = rotor.s
    nitera, loprini = s.nitera, s.loprini
    s.loprini, s.nitera = False, 1
    try:
        hist, rlx, rms = [], [], []
        for _ in range(tol.iterations):
            step()
            hist.append(_values(rotor))
            d = s.dgamold[:s.ii]
            rlx.append(float(s.rlx))
            rms.append(math.sqrt(float(np.mean(d * d / (1.0 + 1.0 / (s.adv * s.adv))))))
    finally:
        s.nitera, s.loprini = nitera, loprini
    arr = {k: np.array([h[k] for h in hist]) for k in _KEYS}
    if not all(np.all(np.isfinite(a)) for a in arr.values()):
        return False, hist[-1], None
    spread = {"deta": float(np.ptp(arr["efficiency"])),
              "dpower_rel": float(np.ptp(arr["power"]) / abs(arr["power"].mean())),
              "dthrust_rel": float(np.ptp(arr["thrust"]) / abs(arr["thrust"].mean())),
              "dbeta_deg": float(np.ptp(arr["beta75_deg"])),
              "min_rlx": min(rlx), "max_step_rms": max(rms)}
    ok = (spread["deta"] < tol.deta and spread["dpower_rel"] < tol.dpower_rel
          and spread["dthrust_rel"] < tol.dthrust_rel and spread["dbeta_deg"] < tol.dbeta_deg
          and spread["min_rlx"] >= tol.min_rlx and spread["max_step_rms"] <= tol.max_step_rms)
    return ok, {k: float(a.mean()) for k, a in arr.items()}, spread


def solve_fixed_pitch(rotor, rpm: float, cycle: CycleTolerance | None = CycleTolerance()) -> Solution:
    """XROTOR RPM at the current blade angle. If it doesn't converge and cycle is not None,
    continue one Newton iteration at a time and accept the mean of cycle.iterations
    iterations if the limit cycle is small. The rotor is left at the last iterate."""
    ok = rotor.solve_rpm(rpm)
    v = _values(rotor)
    if ok:
        return Solution(usable=_finite(v), exact=True, method="rpm", solves=1, **v)
    if cycle is None:
        return Solution(usable=False, exact=False, method=None, solves=1, **v)
    small, v, spread = _cycle(rotor, lambda: rotor.solve_rpm(rpm), cycle)
    return Solution(usable=small, exact=False, method="cycle" if small else None,
                    limit_cycle=spread, solves=1 + cycle.iterations, **v)


def solve_power(rotor, power: float, rpm: float, *, continuation=CONTINUATION, scan=SCAN,
                max_bracket_deg: float = 2.1, max_eta_jump: float = 0.06,
                cycle: CycleTolerance | None = CycleTolerance()) -> Solution:
    """Blade angle giving `power` (W) at `rpm`, as a constant-speed propeller.

    continuation: Mach fractions for step 2 (the first should remove compressibility), or
    None to skip it. scan: (first, last, step) in degrees around the blade angle without
    compressibility (around the starting blade angle if that solve fails too), or None.

    The rotor ends at the solution for "powe" and "continuation". For "scan" it ends at the
    lower bracketing solve (a converged or limit-cycle state, not at the target); when no
    method works it is restored to its starting state."""
    power, rpm = float(power), float(rpm)
    start = rotor.s.copy()
    nsolves = 1
    if rotor.solve_power(power, fix="rpm", rpm=rpm):
        v = _values(rotor)
        return Solution(usable=_finite(v), exact=True, method="powe", solves=nsolves, **v)

    center = None
    if continuation:
        rotor.s = start.copy()
        vso = rotor.s.vso
        ok = False
        try:
            for k, ms in enumerate(list(continuation) + [1.0]):
                rotor.s.vso = vso / ms
                ok = rotor.solve_power(power, fix="rpm", rpm=rpm)
                nsolves += 1
                if k == 0 and ok:
                    center = beta75_deg(rotor)
                if not ok:
                    break
        finally:
            rotor.s.vso = vso
        if ok:
            v = _values(rotor)
            return Solution(usable=_finite(v), exact=True, method="continuation",
                            solves=nsolves, **v)

    out = Solution(usable=False, exact=False, method=None)
    if scan:
        rotor.s = start.copy()
        out = scan_power(rotor, power, rpm, beta75_deg(rotor) if center is None else center,
                         scan=scan, max_bracket_deg=max_bracket_deg, max_eta_jump=max_eta_jump,
                         cycle=cycle)
    out.solves += nsolves
    if out.method is None:
        rotor.s = start
    return out


def scan_power(rotor, power: float, rpm: float, center_deg: float, *, scan=SCAN,
               max_bracket_deg: float = 2.1, max_eta_jump: float = 0.06,
               cycle: CycleTolerance | None = CycleTolerance()) -> Solution:
    """Step 3 of solve_power on its own: fixed-pitch solves (solve_fixed_pitch) from the
    rotor's current state at blade angles (0.75R) center_deg + scan[0] ... + scan[1] every
    scan[2] degrees, and every value interpolated linearly at `power` between the first two
    consecutive usable solves that bracket it. The rotor ends at the lower bracketing solve,
    or in its starting state if nothing brackets the power."""
    power, rpm = float(power), float(rpm)
    start = rotor.s.copy()
    b_start = beta75_deg(rotor)
    out = Solution(usable=False, exact=False, method=None)
    pts = []
    for db in np.arange(scan[0], scan[1] + 1e-9, scan[2]):
        rotor.s = start.copy()
        rotor.change_pitch(float(center_deg) + float(db) - b_start)
        o = solve_fixed_pitch(rotor, rpm, cycle)
        out.solves += o.solves
        if o.usable:
            pts.append((o, rotor.s.copy()))
    out.points = [p[0] for p in pts]
    rotor.s = start
    for (lo, s_lo), (hi, _) in zip(pts, pts[1:]):
        if min(lo.power, hi.power) <= power <= max(lo.power, hi.power) and hi.power != lo.power:
            f = (power - lo.power) / (hi.power - lo.power)
            for k in _KEYS:
                setattr(out, k, getattr(lo, k) + f * (getattr(hi, k) - getattr(lo, k)))
            out.power = power
            out.bracket_deg = hi.beta75_deg - lo.beta75_deg
            out.limit_cycle = lo.limit_cycle or hi.limit_cycle
            out.method = "scan"
            out.usable = (out.bracket_deg <= max_bracket_deg
                          and abs(hi.efficiency - lo.efficiency) <= max_eta_jump)
            rotor.s = s_lo
            break
    return out
