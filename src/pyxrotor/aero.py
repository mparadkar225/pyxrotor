"""Blade section aerodynamics: Drela's parametric CL, CD, CM model (xaero.f).

Ported routines:
    CLCDCM     clcdcm()        section CL, CD, CM and derivatives for one parameter set
    GETCLCDCM  getclcdcm()     the same at a blade station, interpolated between sections
    GETALF     getalf()        inverse: angle of attack for a given CL (Newton)
    SETIAERO   setiaero()      which section each station takes its data from
    SORTAR     sort_sections() order sections by radius

Every expression keeps the Fortran's operation order, so results match XROTOR 7.55 bit for
bit; see CLAUDE.md before "simplifying" anything. Operations that can meet a zero divisor,
overflow, or NaN use pyxrotor.fmath, which returns the IEEE inf/NaN the Fortran produces.

Indexing: stations and sections are 0-based here. An IAERO value N in Fortran is N-1 in
`iaero`; Fortran's "unset" value 0 (from INIT) is -1.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, fields
from typing import NamedTuple, Sequence

import numpy as np

from . import _accel
from . import fmath as fm
from .fmath import div as _div, pow as _pow   # names used throughout this module

# Compressibility drag model constants (CLCDCM, "HHY 10/23/00")
CDMFACTOR = 10.0
CLMFACTOR = 0.25
MEXP = 3.0
CDMDD = 0.0020
CDMSTALL = 0.1000

# GETALF Newton iteration (DATA statements)
GETALF_NITER = 10
GETALF_EPS = 1.0e-5




@dataclass
class AeroSection:
    """One aero section: XIAERO(N) and AERODATA(1:13,N). Defaults are SETDEF's."""
    xisect: float = 0.0           # XIAERO       r/R of the section
    a0: float = 0.0               # AERODATA( 1) zero-lift angle of attack (rad)
    clmax: float = 1.5            # AERODATA( 2) stall CL
    clmin: float = -0.5           # AERODATA( 3) negative stall CL
    dclda: float = 6.28           # AERODATA( 4) lift-curve slope (/rad)
    dclda_stall: float = 0.1      # AERODATA( 5) stalled lift-curve slope (/rad)
    dcl_stall: float = 0.1        # AERODATA( 6) CL increment from incipient to total stall
    cdmin: float = 0.013          # AERODATA( 7) minimum CD
    cldmin: float = 0.5           # AERODATA( 8) CL at minimum CD
    dcdcl2: float = 0.004         # AERODATA( 9) d(CD)/d(CL**2)
    cmcon: float = -0.1           # AERODATA(10) section CM
    reref: float = 200000.0       # AERODATA(11) Reynolds number at which CD applies
    rexp: float = -0.4            # AERODATA(12) exponent for Re scaling of CD
    mcrit: float = 0.8            # AERODATA(13) critical Mach number


# Fortran AERODATA row (1-based) for each AeroSection field: a0 -> 1, ..., mcrit -> 13
AERODATA_ROW = {f.name: i for i, f in enumerate(fields(AeroSection)) if f.name != "xisect"}


class SectionAero(NamedTuple):
    """CLCDCM outputs. clmax/clmin are the Mach-reduced stall limits that CLCDCM writes back
    into its CLMAX/CLMIN arguments."""
    clift: float
    cl_alf: float
    cl_w: float
    stallf: bool
    cdrag: float
    cd_alf: float
    cd_w: float
    cd_rey: float
    cmom: float
    cm_al: float
    cm_w: float
    clmax: float
    clmin: float
    mach_limited: bool        # local Mach was >= 1 and was limited (CLCDCM prints a warning)


class StationAero(NamedTuple):
    """GETCLCDCM outputs at a blade station."""
    clift: float
    cl_alf: float
    cl_w: float
    clmax: float
    clmin: float
    dcl_stall: float
    stallf: bool
    cdrag: float
    cd_alf: float
    cd_w: float
    cd_rey: float
    cmom: float
    cm_al: float
    cm_w: float


def clcdcm(alf: float, w: float, rey: float, vel: float, vso: float, sec: AeroSection) -> SectionAero:
    """Port of SUBROUTINE CLCDCM: CL(alpha), CD(alpha), CM for one section parameter set.

    alf  angle of attack (rad)       w    local speed / VEL       rey  local Reynolds number
    vel  flight speed (COMMON VEL)   vso  speed of sound (COMMON VSO)
    rey <= 0 disables Reynolds-number scaling of CD.
    """
    alf, w, rey, vel, vso = float(alf), float(w), float(rey), float(vel), float(vso)
    if _accel.enabled:
        return SectionAero(*_accel.ext.clcdcm(alf, w, rey, vel, vso, sec))
    a0, clmax, clmin = sec.a0, sec.clmax, sec.clmin
    dclda, dclda_stall, dcl_stall = sec.dclda, sec.dclda_stall, sec.dcl_stall
    cdmin, cldmin, dcdcl2, cmcon = sec.cdmin, sec.cldmin, sec.dcdcl2, sec.cmcon
    mcrit, reref, rexp = sec.mcrit, sec.reref, sec.rexp

    # ---- Prandtl-Glauert compressibility factor
    msq = w * w * (vel * vel) / (vso * vso)
    msq_w = 2.0 * w * (vel * vel) / (vso * vso)
    mach_limited = False
    if msq >= 1.0:
        mach_limited = True                      # 'CLFUNC: Local Mach number limited to 0.99'
        msq = 0.99
        msq_w = 0.0
    pg = 1.0 / fm.sqrt(1.0 - msq)
    pg_w = 0.5 * msq_w * (pg * pg * pg)

    # ---- Mach number and dependence on velocity
    mach = fm.sqrt(msq)
    mach_w = 0.0
    if mach != 0.0:
        mach_w = 0.5 * msq_w / mach

    # --- CL from dCL/dAlpha and Prandtl-Glauert scaling
    cla = dclda * pg * (alf - a0)
    cla_alf = dclda * pg
    cla_w = dclda * pg_w * (alf - a0)

    # --- Effective CLmax is limited by Mach effects
    dmstall = _pow(CDMSTALL / CDMFACTOR, 1.0 / MEXP)
    clmaxm = fm.fmax(0.0, (mcrit + dmstall - mach) / CLMFACTOR) + cldmin
    clmax = fm.fmin(clmax, clmaxm)
    clminm = fm.fmin(0.0, -(mcrit + dmstall - mach) / CLMFACTOR) + cldmin
    clmin = fm.fmax(clmin, clminm)

    # --- CL limiter function (turns on after +-stall)
    ecmax = fm.exp(fm.fmin(200.0, _div(cla - clmax, dcl_stall)))
    ecmin = fm.exp(fm.fmin(200.0, _div(clmin - cla, dcl_stall)))
    cllim = dcl_stall * fm.log((1.0 + ecmax) / (1.0 + ecmin))
    cllim_cla = ecmax / (1.0 + ecmax) + ecmin / (1.0 + ecmin)

    # --- Subtract a (nearly unity) fraction of the limited CL function
    fstall = _div(dclda_stall, dclda)
    clift = cla - (1.0 - fstall) * cllim
    cl_alf = cla_alf - (1.0 - fstall) * cllim_cla * cla_alf
    cl_w = cla_w - (1.0 - fstall) * cllim_cla * cla_w

    stallf = clift > clmax or clift < clmin

    # --- CM from CMCON and Prandtl-Glauert scaling
    cmom = pg * cmcon
    cm_al = 0.0
    cm_w = pg_w * cmcon

    # ---- Reynolds number scaling factor
    if rey <= 0:
        rcorr = 1.0
        rcorr_rey = 0.0
    else:
        rcorr = _pow(rey / reref, rexp)
        rcorr_rey = rexp / rey

    # --- Linear-range drag: CD0 + quadratic in (CL - CLDMIN)
    dcl = clift - cldmin
    cdrag = (cdmin + dcdcl2 * (dcl * dcl)) * rcorr
    cd_alf = (2.0 * dcdcl2 * dcl * cl_alf) * rcorr
    cd_w = (2.0 * dcdcl2 * dcl * cl_w) * rcorr
    cd_rey = cdrag * rcorr_rey

    # --- Post-stall drag
    fstall = _div(dclda_stall, dclda)
    dcdx = _div((1.0 - fstall) * cllim, pg * dclda)
    dcd = 2.0 * (dcdx * dcdx)
    dcd_alf = _div(4.0 * dcdx * (1.0 - fstall) * cllim_cla * cla_alf, pg * dclda)
    dcd_w = 4.0 * dcdx * (_div((1.0 - fstall) * cllim_cla * cla_w, pg * dclda) - dcd / pg * pg_w)

    # --- Compressibility drag above an effective Mcrit that falls with |CL - CLDMIN|
    dmdd = _pow(CDMDD / CDMFACTOR, 1.0 / MEXP)
    critmach = mcrit - CLMFACTOR * abs(clift - cldmin) - dmdd
    critmach_alf = -CLMFACTOR * abs(cl_alf)
    critmach_w = -CLMFACTOR * abs(cl_w)
    if mach < critmach:
        cdc = 0.0
        cdc_alf = 0.0
        cdc_w = 0.0
    else:
        cdc = CDMFACTOR * _pow(mach - critmach, MEXP)
        cdc_w = _div(MEXP * mach_w * cdc, mach) - _div(MEXP * critmach_w * cdc, critmach)
        cdc_alf = -_div(MEXP * critmach_alf * cdc, critmach)

    # --- Total drag (FAC is XROTOR's hook for Mach-scaled profile drag; kept at 1)
    fac = 1.0
    fac_w = 0.0
    cdrag = fac * cdrag + dcd + cdc
    cd_alf = fac * cd_alf + dcd_alf + cdc_alf
    cd_w = fac * cd_w + fac_w * cdrag + dcd_w + cdc_w
    cd_rey = fac * cd_rey

    return SectionAero(clift, cl_alf, cl_w, stallf, cdrag, cd_alf, cd_w, cd_rey,
                       cmom, cm_al, cm_w, clmax, clmin, mach_limited)


def setiaero(xi: Sequence[float], sections: Sequence[AeroSection]) -> np.ndarray:
    """Port of SUBROUTINE SETIAERO: index of the inner bounding section for each station."""
    iaero = np.zeros(len(xi), dtype=np.int64)
    for i, x in enumerate(xi):
        iaero[i] = 0
        for n, sec in enumerate(sections):
            if sec.xisect <= x:
                iaero[i] = n
    return iaero


def sort_sections(sections: Sequence[AeroSection]) -> list[AeroSection]:
    """Port of SUBROUTINE SORTAR applied to XIAERO/AERODATA: stable sort by r/R.
    (SORTAR is a bubble sort that swaps only strictly out-of-order pairs, so it is stable.)"""
    return sorted(sections, key=lambda s: s.xisect)


def getclcdcm(is_: int, alf: float, w: float, rey: float, *, xi: Sequence[float],
              iaero: np.ndarray, sections: Sequence[AeroSection], vel: float, vso: float,
              fix_fallback: bool = False) -> StationAero:
    """Port of SUBROUTINE GETCLCDCM: section properties at station is_ (0-based).

    Data comes from section iaero[is_] and, if there is one, the next section outward,
    linearly interpolated in r/R. Like the Fortran, the interpolation fraction is not
    clipped, so a stale iaero extrapolates (see "Known issues" in CLAUDE.md).

    If iaero[is_] is not a valid section index, it is recomputed and written back to
    iaero, faithfully including 7.55's fallback behavior: the search leaves N at the first
    section *outboard* of the station (the outer, not the inner, bounding section), and
    stores the inner index in IAERO only if the station lies outboard of section 1. A station
    outboard of the last section finds no bounding pair and gets section 1 (the root one).

    fix_fallback=True instead takes the inner bounding section, as SETIAERO would (section 1
    for a station inboard of it), and stores it in IAERO.
    """
    naero = len(sections)
    xis = float(xi[is_])

    # --- Check for installed aero data section index
    n = int(iaero[is_])
    if (n < 0 or n >= naero) and fix_fallback:
        n = 0
        for k in range(naero):
            if sections[k].xisect <= xis:
                n = k
        iaero[is_] = n
    if n < 0 or n >= naero:                           # 7.55's search
        jumped = False
        if naero > 1:
            # --- Find lower index of aero data sections XIAERO(N) bounding XI(IS)
            for n in range(naero):
                if sections[n].xisect <= xis:
                    iaero[is_] = n
                else:
                    jumped = True                     # GO TO 10 with N at this section
                    break
            # (otherwise XROTOR prints 'Aero section not found for station')
        if not jumped:
            n = 0
            iaero[is_] = n

    # --- Data for the inner bounding section
    sec1 = sections[n]
    xisect1 = sec1.xisect
    r1 = clcdcm(alf, w, rey, vel, vso, sec1)
    clift, cl_alf, cl_w, stallf = r1.clift, r1.cl_alf, r1.cl_w, r1.stallf
    cdrag, cd_alf, cd_w, cd_rey = r1.cdrag, r1.cd_alf, r1.cd_w, r1.cd_rey
    cmom, cm_al, cm_w = r1.cmom, r1.cm_al, r1.cm_w
    clmax, clmin, dcl_stall = r1.clmax, r1.clmin, sec1.dcl_stall

    # --- If there is an outer bounding section, interpolate linearly to the station
    if n < naero - 1:
        sec2 = sections[n + 1]
        frac = _div(xis - xisect1, sec2.xisect - xisect1)
        r2 = clcdcm(alf, w, rey, vel, vso, sec2)
        stallf = stallf or r2.stallf
        clift = (1.0 - frac) * clift + frac * r2.clift
        cl_alf = (1.0 - frac) * cl_alf + frac * r2.cl_alf
        cl_w = (1.0 - frac) * cl_w + frac * r2.cl_w
        clmax = (1.0 - frac) * clmax + frac * r2.clmax
        clmin = (1.0 - frac) * clmin + frac * r2.clmin
        dcl_stall = (1.0 - frac) * dcl_stall + frac * sec2.dcl_stall

        cmom = (1.0 - frac) * cmom + frac * r2.cmom
        cm_al = (1.0 - frac) * cm_al + frac * r2.cm_al
        cm_w = (1.0 - frac) * cm_w + frac * r2.cm_w

        cdrag = (1.0 - frac) * cdrag + frac * r2.cdrag
        cd_alf = (1.0 - frac) * cd_alf + frac * r2.cd_alf
        cd_w = (1.0 - frac) * cd_w + frac * r2.cd_w
        cd_rey = (1.0 - frac) * cd_rey + frac * r2.cd_rey

    return StationAero(clift, cl_alf, cl_w, clmax, clmin, dcl_stall, stallf,
                       cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w)


class AlphaSolution(NamedTuple):
    alf: float
    alf_cl: float
    alf_w: float
    stallf: bool
    converged: bool           # False where XROTOR prints 'GETALF: ... inversion failed'


def getalf(is_: int, clift: float, w: float, *, xi: Sequence[float], iaero: np.ndarray,
           sections: Sequence[AeroSection], vel: float, vso: float,
           fix_fallback: bool = False) -> AlphaSolution:
    """Port of SUBROUTINE GETALF: angle of attack giving CL = clift at station is_.

    Newton iteration on GETCLCDCM with Reynolds scaling off, started from section 1's
    zero-lift angle at every station (as in 7.55). Up to 10 iterations, |dalpha| < 1e-5.
    """
    clift = float(clift)
    stallf = False
    a0 = sections[0].a0
    rey = 0.0

    alf = a0
    alf_cl = alf_w = 0.0
    for _ in range(GETALF_NITER):
        r = getclcdcm(is_, alf, w, rey, xi=xi, iaero=iaero, sections=sections, vel=vel, vso=vso,
                      fix_fallback=fix_fallback)
        stallf = r.stallf
        dalf = -_div(r.clift - clift, r.cl_alf)
        alf = alf + dalf
        alf_cl = _div(1.0, r.cl_alf)
        alf_w = -_div(r.cl_w, r.cl_alf)
        if abs(dalf) < GETALF_EPS:
            return AlphaSolution(alf, alf_cl, alf_w, stallf, True)
    return AlphaSolution(alf, alf_cl, alf_w, stallf, False)
