"""Local velocity triangle and wake ("equivalent prop") geometry (xoper.f, xrotor.f).

Ported routines:
    UVADD   uvadd()    externally added axial/tangential inflow at a radius
    CSCALC  cscalc()   velocity components, W and phi at a station, with sensitivities
    XWINIT  xwinit()   initial equivalent-prop radii XW and wake swirl VWAK
    SETXW   setxw()    equivalent-prop radii and swirl from angular momentum, with
                       sensitivities to ADV, ADW and every GAM(J)

CSCALC and SETXW are compiled (_speedups.pyx); the functions here call it.

Not ported: ADWCORR, which computes a local array and is never called in XROTOR 7.55.

All take the State and a 0-based station index; arrays follow State's conventions
(Fortran VIND(3,I) is s.vind[2, i]).
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from . import _speedups
from . import fmath as fm
from .fmath import div
from .numerics import seval
from .state import State


def uvadd(s: State, xiw: float) -> tuple[float, float]:
    """Port of SUBROUTINE UVADD(XIW,WA,WT): added inflow (normalized by VEL) at r/R = xiw,
    zero outside the tabulated radius range or with fewer than two points."""
    if s.nadd <= 1:
        return 0.0, 0.0
    rdim = float(xiw) * s.rad
    n = s.nadd
    if s.radd[0] <= rdim <= s.radd[n - 1]:
        wa = div(seval(rdim, s.uadd, s.uaddr, s.radd, n), s.vel)
        wt = div(seval(rdim, s.vadd, s.vaddr, s.radd, n), s.vel)
        return wa, wt
    return 0.0, 0.0


class Velocities(NamedTuple):
    """CSCALC outputs (Fortran argument order)."""
    utot: float
    wa: float
    wt: float
    vt: float
    vt_adw: float
    va: float
    va_adw: float
    vd: float
    vd_adw: float
    ci: float
    ci_adv: float
    ci_vt: float
    si: float
    si_va: float
    w: float
    w_adv: float
    w_vt: float
    w_va: float
    phi: float
    p_adv: float
    p_vt: float
    p_va: float


def cscalc(s: State, i: int) -> Velocities:
    """Port of SUBROUTINE CSCALC: velocity components at station i on the real prop."""
    return Velocities(*_speedups.cscalc(s, i))


def xwinit(s: State) -> None:
    """Port of SUBROUTINE XWINIT: initial estimate of the equivalent-prop radii XW and the
    wake swirl VWAK, from mass continuity with the undisturbed axial velocity."""
    uduct = 0.0
    if s.duct:
        uduct = s.urduct - 1.0
    ii = s.ii
    xm = s.xw0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for i in range(ii):
            xi = float(s.xi[i])
            urat = 1.0 + uduct + float(s.ubody[i])
            s.dxw[i] = fm.sqrt(xm * xm + 2.0 * urat * xi * float(s.dxi[i])) - xm
            xp = xm + float(s.dxw[i])
            s.xw[i] = 0.5 * (xp + xm)

            s.xw_adv[i] = 0.0
            s.xw_adw[i] = 0.0
            s.xw_gam[i, :ii] = 0.0

            xw = float(s.xw[i])
            s.vwak[i] = div(float(s.vind[2, i]) * xi, xw)
            s.vw_adv[i] = 0.0
            s.vw_adw[i] = div(float(s.vind_adw[2, i]) * xi, xw)
            s.vw_gam[i, :ii] = s.vind_gam[2, i, :ii] * xi / xw
            xm = xp
    s.xwtip = xm


def setxw(s: State) -> None:
    """Port of SUBROUTINE SETXW: equivalent-prop radius XW and swirl VWAK at each station.

    The radial stream function (axial velocity * r dr) sets the equivalent prop's radial
    coordinate; angular momentum is conserved to set its swirl. Also returns, in the State,
    the sensitivities of XW, DXW and VWAK to ADV, ADW and each circulation GAM(J).
    """
    _speedups.setxw(s)
