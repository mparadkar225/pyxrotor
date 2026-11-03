"""OPER case sweeps and the case table (xoper.f, xcasepl.f).

Ported:
    SETCAS      setcas()       fill the case table for ASEQ (1), RSEQ (2), BSEQ (3), VSEQ (4)
    CASPLT      run_cases()    solve every stored case, then restore and re-solve the current
                               point (its computing part; the plots and scratch W1-W4 are not)
    OPER ADDC   add_case()     current operating point to the table
    OPER CLRC   clear_cases()
    OPER CASE   rerun_case()   rerun one stored case at the current blade
    SHOCAS      case_lines()   the case listing (OPER display and CPUT file format)
    GETCAS      read_case_lines()   CGET

The table is CASPAR(0:12, 100): CASPAR(0) encodes the case type, ITYP + 100*(engine
power-rpm line) + 1000*(RSEQ target: 1 thrust, 2 torque, 3 power), then ADV, VEL, tip BETA,
ALT, RHO, RMU, VSO and, once solved, power, thrust, torque, efficiency (999 = not solved).
Python index: s.caspar[k, n] is CASPAR(k, n+1).

7.55 behavior reproduced (see CLAUDE.md):
  * SETCAS's increments are named DADV and DBET, which are COMMON (the solver's ADV and
    blade-angle Newton steps), so ASEQ leaves its step in DADV and BSEQ its step, in degrees,
    in DBET.
  * CASPLT solves each non-engine case twice: APER(4,2) first, then the type's own solve.
  * ADDC and GETCAS do not set CASPAR(0), so such cases keep whatever type code was there.
  * GETCAS reads lines into 128 characters, but a SHOCAS line with results is 159 long, so
    CGET drops power, thrust, torque and efficiency (fix: read_case_lines(full_lines=True)).
"""
from __future__ import annotations

import math

import numpy as np

from . import fmath as fm
from .constants import ICASX, PI
from .fmath import div
from .fortio import fmt_g, fmt_i, getflt
from .oper import aper
from .state import State

PI_EXACT = 4.0 * math.atan(1.0)          # SHOCAS and GETCAS compute PI = 4*ATAN(1)
NOT_SOLVED = 999.0


def _fill(s: State, n: int, code: float | None, adv: float, vel: float, beta_tip: float) -> None:
    c = s.caspar
    if code is not None:
        c[0, n] = code
    c[1, n], c[2, n], c[3, n] = adv, vel, beta_tip
    c[4, n], c[5, n], c[6, n], c[7, n] = s.alt, s.rho, s.rmu, s.vso
    c[8:12, n] = NOT_SOLVED


def _count(first: float, last: float, step: float) -> tuple[float, int]:
    """SIGN the step toward last; NP = INT((last-first)/step + 0.5) + 1 (INTEGER*4)."""
    step = fm.sign(step, last - first)
    np_ = 1
    if step != 0.0:
        np_ = fm.iwrap(fm.ifix(div(last - first, step) + 0.5) + 1)
    return step, np_


def setcas(s: State, itype: int, first: float, last: float, step: float, *,
           engine_line: bool = False, fix: str = "A", constant_speed: bool = False,
           rpm: float | None = None, power: float | None = None) -> int:
    """Port of SUBROUTINE SETCAS(ITYPE,NINPUT,RINPUT) with all three values given (NINPUT=3).
    Appends cases to the table; returns how many were added.

    engine_line: answer to 'Use engine rpm/power line ?' (asked only if LPWRVAR).
    fix (RSEQ): 'T' thrust, 'Q' torque, 'P' power held at the current spec, else blade pitch.
    constant_speed (VSEQ): CS at rpm (default the current one) and power (default PSPEC,
    or the current power if PSPEC <= 0), instead of FP fixed pitch.
    """
    first, last, step = float(first), float(last), float(step)
    n0 = s.ncase
    s.kcase = 0
    if itype not in (1, 2, 3, 4):
        raise ValueError("itype must be 1 (ADV), 2 (RPM), 3 (blade angle) or 4 (velocity)")
    s.kcase = itype

    def common_step(v):                                # 7.55: the increment is COMMON
        if itype == 1:
            s.dadv = v                                 #   DADV (ASEQ)
        elif itype == 3:
            s.dbet = v                                 #   DBET (BSEQ, in degrees)

    common_step(step)                                  # set before the early RETURN
    if first == last:
        return 0
    step, np_ = _count(first, last, step)
    common_step(step)
    if np_ <= 0:
        return 0

    yes = bool(engine_line) if s.lpwrvar else False
    xans = 100.0 if yes else 0.0

    if itype == 2:
        f = (fix or " ").upper()[:1]
        xans = xans + {"T": 1000.0, "Q": 2000.0, "P": 3000.0}.get(f, 0.0)
    if itype == 4 and constant_speed:
        r = div(s.vel, ((s.rad * s.adv) * PI) / 30.0)
        if rpm is not None:
            r = float(rpm)
        s.adv = div(s.vel, ((s.rad * r) * PI) / 30.0)
        s.kcase = 5
        if xans != 100.0:
            if s.pspec <= 0.0 and s.ptot > 0.0:
                s.pspec = s.ptot * ((s.rho * (s.vel * s.vel * s.vel)) * (s.rad * s.rad))
            if power is not None:
                s.pspec = float(power)

    if s.ncase + np_ > ICASX:                          # 'Limiting number of cases...'
        np_ = ICASX - s.ncase
    code = xans + float(s.kcase)
    beta_tip = float(s.beta[s.ii - 1])
    for ip in range(np_):
        n = s.ncase
        s.ncase += 1
        x = first + step * float(ip)
        if itype == 1:
            _fill(s, n, code, x, s.vel, beta_tip)
        elif itype == 2:
            _fill(s, n, code, ((div(s.vel, x * s.rad)) * 30.0) / PI, s.vel, beta_tip)
        elif itype == 3:
            _fill(s, n, code, s.adv, s.vel, (x * PI) / 180.0)
        else:
            _fill(s, n, code, div(x * s.adv, s.vel), x, beta_tip)
    return s.ncase - n0


def _results(s: State, n: int) -> None:
    v, r = s.vel, s.rad
    c = s.caspar
    if s.conv:
        c[8, n] = ((s.ptot * s.rho) * (v * v * v)) * (r * r)
        c[9, n] = ((s.ttot * s.rho) * (v * v)) * (r * r)
        c[10, n] = ((s.qtot * s.rho) * (v * v)) * (r * r * r)
        c[11, n] = div(s.ttot, s.ptot)
    else:
        c[8:12, n] = NOT_SOLVED


def decode_type(code: float) -> tuple[int, int, int]:
    """CASPAR(0) -> (ITYP, IPWR, ITQP), with CASPLT's INTEGER*4 arithmetic."""
    itqp = fm.ifix(div(code, 1000.0))
    ipwr = fm.iwrap(fm.ifix(div(code, 100.0)) - fm.iwrap(10 * itqp))
    ityp = fm.iwrap(fm.iwrap(fm.ifix(code) - fm.iwrap(100 * ipwr)) - fm.iwrap(1000 * itqp))
    return ityp, ipwr, itqp


def run_cases(s: State) -> None:
    """The computing part of SUBROUTINE CASPLT: solve each stored case from the saved blade,
    store results in the table, then restore the operating point and re-solve it."""
    if s.ncase == 0:
        return
    ii = s.ii
    advsav, velsav, altsav, rhosav, vsosav, rmusav = s.adv, s.vel, s.alt, s.rho, s.vso, s.rmu
    betsav = s.beta[:ii].copy()
    bet0sav = s.beta0[:ii].copy()
    c = s.caspar
    for n in range(s.ncase):
        ityp, ipwr, itqp = decode_type(float(c[0, n]))
        s.adv, s.vel = float(c[1, n]), float(c[2, n])
        bet = float(c[3, n])
        s.alt, s.rho, s.rmu, s.vso = float(c[4, n]), float(c[5, n]), float(c[6, n]), float(c[7, n])
        delb = bet - float(betsav[ii - 1])
        s.beta[:ii] = betsav + delb
        s.beta0[:ii] = bet0sav + delb

        s.conv = False
        if ipwr == 0:
            if ityp != 5:                              # 7.55: an extra solve first
                aper(s, 4, 2, True)
            else:
                aper(s, 3, 1, True)
            if ityp == 1:
                aper(s, 4, 2, True)
            elif ityp == 2:
                if itqp == 1:
                    aper(s, 1, 1, True)
                elif itqp == 2:
                    aper(s, 2, 1, True)
                elif itqp == 3:
                    aper(s, 3, 1, True)
                else:
                    aper(s, 4, 2, True)
            elif ityp in (3, 4):
                aper(s, 4, 2, True)
            elif ityp == 5:
                aper(s, 3, 1, True)
        else:
            if ityp in (1, 2, 5):
                aper(s, 5, 1, True)
            elif ityp in (3, 4):
                aper(s, 5, 2, True)

        c[1, n] = s.adv
        c[3, n] = s.beta[ii - 1]
        _results(s, n)

    s.adv, s.vel, s.alt, s.vso, s.rho, s.rmu = advsav, velsav, altsav, vsosav, rhosav, rmusav
    s.beta[:ii] = betsav
    s.beta0[:ii] = bet0sav
    s.conv = False
    aper(s, 4, 2, True)


def add_case(s: State) -> bool:
    """OPER ADDC. Returns False if the table is full. (CASPAR(0) is left as it was.)"""
    if s.ncase >= ICASX:
        return False
    n = s.ncase
    s.ncase += 1
    _fill(s, n, None, s.adv, s.vel, float(s.beta[s.ii - 1]))
    return True


def clear_cases(s: State) -> None:
    """OPER CLRC."""
    s.ncase = 0
    s.kcase = 0


def rerun_case(s: State, icase: int) -> bool:
    """OPER CASE i (1-based): set the stored operating point and tip angle, solve at fixed
    advance ratio, and store the results if it converged. Returns CONV."""
    icase = int(icase)
    if s.ncase <= 0 or icase <= 0 or icase > s.ncase:
        return False
    n = icase - 1
    c = s.caspar
    s.adv, s.vel = float(c[1, n]), float(c[2, n])
    bet = float(c[3, n])
    s.alt, s.rho, s.rmu, s.vso = float(c[4, n]), float(c[5, n]), float(c[6, n]), float(c[7, n])
    ii = s.ii
    delb = bet - float(s.beta[ii - 1])
    s.beta[:ii] = s.beta[:ii] + delb
    s.beta0[:ii] = s.beta0[:ii] + delb
    s.conv = False
    aper(s, 4, 2, True)
    if s.conv:
        _results(s, n)
    return bool(s.conv)


# --------------------------------------------------------------------------- listing
HEADER = ("  n"
          "         V/wR         Btip            V          rpm"
          "          rho       mu*1e5       Vsound            h"
          "        P(kW)         T(N)       Q(N-m)          eff")


def case_lines(caspar: np.ndarray, ncase: int, rad: float, name: str) -> list[str]:
    """SUBROUTINE SHOCAS: the case listing, one string per record (CPUT writes these)."""
    out = [name.ljust(32)[:32], HEADER, "-" * 160]
    for n in range(ncase):
        p = [float(v) for v in caspar[:, n]]
        adv, vel = p[1], p[2]
        bet = (p[3] * 180.0) / PI_EXACT
        rmu = p[6] * 1.0e5
        rpm = (div(vel, rad * adv) * 30.0) / PI_EXACT
        vals = [adv, bet, vel, rpm, p[5], rmu, p[7], p[4]]
        if p[8] != NOT_SOLVED:
            vals += [p[8] * 0.001, p[9], p[10], p[11]]
        out.append(fmt_i(n + 1, 3) + "".join(" " + fmt_g(v, 12, 5) for v in vals))
    return out


def read_case_lines(lines: list[str], full_lines: bool = False) -> np.ndarray | None:
    """SUBROUTINE GETCAS: parse a case listing (after its 3 header records) into CASPAR
    columns 1..11 for each case; returns an (11, n) array, or None on a read error (7.55
    then sets NCASE = 0). full_lines=False keeps 7.55's 128-character records."""
    body = lines[3:]
    rows = []
    a = None
    for raw in body:
        line = raw if full_lines else raw[:128]
        if full_lines:
            vals, err = getflt(line.ljust(len(line) + 1), 13, limit=len(line) + 1)
        else:
            vals, err = getflt(line.ljust(128), 13)
        if err:
            return None
        n = len(vals)
        if n == 0:
            if a is None:
                raise ValueError("GETCAS: blank first case line (7.55 reads garbage)")
        else:
            a = (vals + [0.0] * 16)[:16] if a is None else vals + a[n:]
        adv, bet, vel = a[1], (a[2] * PI_EXACT) / 180.0, a[3]
        rho, rmu, vso, alt = a[5], div(a[6], 1.0e5), a[7], a[8]
        pow_ = thr = trq = eff = NOT_SOLVED
        if n == 13:
            pow_, thr, trq, eff = a[9] * 1000.0, a[10], a[11], a[12]
        rows.append([adv, vel, bet, alt, rho, rmu, vso, pow_, thr, trq, eff])
    if len(rows) > ICASX:
        raise ValueError(f"GETCAS: more than {ICASX} cases (7.55 writes past CASPAR)")
    return np.array(rows, dtype=float).T.reshape(11, len(rows))


def cases_from_lines(s: State, lines: list[str], full_lines: bool = False) -> None:
    """OPER CGET: replace the case table (CASPAR(0) not set, as in 7.55)."""
    par = read_case_lines(lines, full_lines)
    if par is None:
        s.ncase = 0
        return
    n = par.shape[1]
    s.caspar[1:12, :n] = par
    s.ncase = n
