"""Local velocity triangle and wake ("equivalent prop") geometry (xoper.f, xrotor.f).

Ported routines:
    UVADD   uvadd()    externally added axial/tangential inflow at a radius
    CSCALC  cscalc()   velocity components, W and phi at a station, with sensitivities
    XWINIT  xwinit()   initial equivalent-prop radii XW and wake swirl VWAK
    SETXW   setxw()    equivalent-prop radii and swirl from angular momentum, with
                       sensitivities to ADV, ADW and every GAM(J)

Not ported: ADWCORR, which computes a local array and is never called in XROTOR 7.55.

All take the State and a 0-based station index; arrays follow State's conventions
(Fortran VIND(3,I) is s.vind[2, i]).
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

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
    vt = float(s.vind[2, i])
    vt_adw = float(s.vind_adw[2, i])
    va = float(s.vind[0, i])
    va_adw = float(s.vind_adw[0, i])

    # ---- duct effect on freestream and induced axial velocity
    uduct = 0.0
    vaduct_va = 1.0
    if s.duct:
        uduct = s.urduct - 1.0
        vaduct_va = 2.0 * s.urduct
    vd = va * (vaduct_va - 1.0)
    vd_va = vaduct_va - 1.0
    vd_adw = vd_va * va_adw

    # ---- freestream, body-induced and added inflow velocities
    xi = float(s.xi[i])
    utot = 1.0 + uduct + float(s.ubody[i])
    wa, wt = uvadd(s, xi)

    adv = s.adv
    ci = div(xi, adv) - wt - vt
    ci_adv = -div(xi, adv * adv)
    ci_vt = -1.0

    si = utot + wa + va + vd
    si_va = 1.0 + vd_va

    wsq = ci * ci + si * si
    w = fm.sqrt(wsq)
    w_adv = div(ci * ci_adv, w)
    w_vt = div(ci * ci_vt, w)
    w_va = div(si * si_va, w)

    phi = fm.atan2(si, ci)
    p_adv = div(-(si * ci_adv), wsq)
    p_vt = div(-(si * ci_vt), wsq)
    p_va = div(ci * si_va, wsq)

    return Velocities(utot, wa, wt, vt, vt_adw, va, va_adw, vd, vd_adw, ci, ci_adv, ci_vt,
                      si, si_va, w, w_adv, w_vt, w_va, phi, p_adv, p_vt, p_va)


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
    ii = s.ii
    adv, adw, urduct = s.adv, s.adw, s.urduct
    xwm = s.xw0
    xwm_gam = np.zeros(ii)
    xwm_adv = 0.0
    xwm_adw = 0.0

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for i in range(ii):
            xi = float(s.xi[i])
            xdx = xi * float(s.dxi[i])
            v = cscalc(s, i)
            vt, wa, si = v.vt, v.wa, v.si

            # ---- first guess for XWO
            dxwo = fm.sqrt(xwm * xwm + 2.0 * v.utot * xdx) - xwm
            xwo = xwm + 0.5 * dxwo

            # ---- Newton loop for XWO (CW, the equivalent-prop swirl, is computed but
            #      never used in the Fortran and is omitted here)
            for _ in range(30):
                vw = div(vt * xi, xwo)
                vw_xwo = -div(vw, xwo)
                utotw = urduct
                vaw = div(vw * xwo, adw)
                vaw_vw = div(xwo, adw)
                vaw_xwo = div(vw, adw) + vaw_vw * vw_xwo
                sw = utotw + wa + vaw
                sw_xwo = vaw_xwo
                rez = sw * xwo * 2.0 * (xwo - xwm) - si * xdx
                rez_xwo = sw * 2.0 * (2.0 * xwo - xwm) + sw_xwo * xwo * 2.0 * (xwo - xwm)
                delxwo = -div(rez, rez_xwo)
                # RLX here is the COMMON variable (OPER_X), so SETXW leaves its last value
                # in the solver state; reproduced deliberately.
                rlx = 1.0
                if abs(delxwo) > 0.2 * (xwo - xwm):
                    rlx = div(0.2 * (xwo - xwm), abs(delxwo))
                s.rlx = rlx
                xwo = xwo + rlx * delxwo
                if abs(delxwo) < 1.0e-6:
                    break
            # (else XROTOR prints 'SETXW: Xw convergence failed.')

            # ---- Vw( xwo , Vt(Adw Gj) )
            vw = div(vt * xi, xwo)
            vw_xwo = -div(vw, xwo)
            vw_vt = div(xi, xwo)

            # ---- axial velocity on equivalent prop (derived from swirl)
            utotw = urduct
            vaw = div(vw * xwo, adw)
            vaw_vw = div(xwo, adw)
            vaw_xwo = div(vw, adw)
            vaw_adw = -div(vaw, adw)

            sw = utotw + wa + vaw
            sw_xwo = vaw_xwo + vaw_vw * vw_xwo
            sw_vt = vaw_vw * vw_vt
            sw_adw = vaw_adw

            # ---- Res( xwo, xwm, Sw(Adw Vt xw), S(Adw Vt) )
            z_xwo = sw * 2.0 * (2.0 * xwo - xwm)
            z_xwm = -(sw * xwo * 2.0)
            z_sw = xwo * 2.0 * (xwo - xwm)
            z_si = -xdx

            z_xwo = z_sw * sw_xwo + z_xwo
            z_vt = z_sw * sw_vt
            z_va = z_si * v.si_va
            z_adw = z_sw * sw_adw

            # ---- Res( xwo, Adv, Adw, Gj )
            z_adv = z_xwm * xwm_adv
            z_adw = z_xwm * xwm_adw + z_vt * v.vt_adw + z_va * v.va_adw + z_adw
            z_gam = z_xwm * xwm_gam + z_vt * s.vind_gam[2, i, :ii] + z_va * s.vind_gam[0, i, :ii]

            # ---- xwo( Adv, Adw, Gj )
            s.xw_adv[i] = -div(z_adv, z_xwo)
            s.xw_adw[i] = -div(z_adw, z_xwo)
            s.xw_gam[i, :ii] = -(z_gam / z_xwo)

            # ---- Vw( Adv, Adw, Gj )
            s.vwak[i] = vw
            s.vw_adv[i] = vw_xwo * s.xw_adv[i]
            s.vw_adw[i] = vw_xwo * s.xw_adw[i] + vw_vt * v.vt_adw
            s.vw_gam[i, :ii] = vw_xwo * s.xw_gam[i, :ii] + vw_vt * s.vind_gam[2, i, :ii]

            s.xw[i] = xwo

            # ---- dxw( xwo(Adv Adw Gj), xwm(Adv Adw Gj) )
            s.dxw[i] = 2.0 * (xwo - xwm)
            s.dxw_adv[i] = 2.0 * (s.xw_adv[i] - xwm_adv)
            s.dxw_adw[i] = 2.0 * (s.xw_adw[i] - xwm_adw)
            s.dxw_gam[i, :ii] = 2.0 * (s.xw_gam[i, :ii] - xwm_gam)

            # ---- new xwm( Adv Adw Gj ) for the next station
            xwm = 2.0 * xwo - xwm
            xwm_adv = 2.0 * float(s.xw_adv[i]) - xwm_adv
            xwm_adw = 2.0 * float(s.xw_adw[i]) - xwm_adw
            xwm_gam = 2.0 * s.xw_gam[i, :ii] - xwm_gam

    s.xwtip = xwm
