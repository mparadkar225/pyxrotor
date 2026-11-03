"""Minimum-induced-loss rotor design (xdesi.f).

Ported routines:
    SETCLD  setcld()   linear design-CL distribution, root to tip
    DESINI  desini()   initial circulation and chord (Larrabee, then Goldstein if not FAST)
    DEITER  deiter()   Newton iteration for the MIL design
    DESMIL  desmil()   = DEITER (the Fortran only differs in its failure message)
    DESGEN  desgen()   set up the design point from the DESI parameters and design the rotor

ICASE 2: fix the CL distribution, vary twist and chord (what DESGEN uses).
ICASE 3: fix the chord distribution, vary twist and CL.

Not ported: SETCLD's file option (CL1 = 999), SAVCLD, MODCL (cursor editing), and OUTPUT
(DESGEN's printout; it changes only scratch arrays).

Ducted designs: DEITER's DUCT branch reads GAMDES, which is not in COMMON and is never
initialized in 7.55, so its result is undefined; pyxrotor raises instead of guessing.
"""
from __future__ import annotations

import numpy as np

from . import fmath as fm
from .aero import getalf, getclcdcm
from .constants import IQ, PI
from .fmath import div
from .induced import helico, vcalc, vrtxco
from .numerics import filter_, gauss
from .oper import induced_model, sections, tpq
from .state import State
from .geometry import setx
from .wake import setxw, uvadd, xwinit


class BetzLimitError(ValueError):
    """DESGEN's 'Disk loading above ultimate Betz limits' check failed."""


def setcld(s: State, cl1: float, cl2: float) -> None:
    """Port of SUBROUTINE SETCLD: CLDES(I) = (1-XI)*CL1 + XI*CL2 on the *current* XI grid.
    (7.55: before a design that is the previous grid; on a first design all XI are 0, so
    every station gets CL1. s.cldes_linear records (CL1, CL2) for Fixes.linear_cl.)"""
    if cl1 == 999.0:
        raise NotImplementedError("SETCLD: reading the CL distribution from a file is not ported")
    cl1, cl2 = float(cl1), float(cl2)
    ii = s.ii
    s.cldes[:ii] = (1.0 - s.xi[:ii]) * cl1 + s.xi[:ii] * cl2
    s.cldes_linear = (cl1, cl2)


def _past_stall(s: State) -> bool:
    sec1 = sections(s)[0]
    return any(c > sec1.clmax or c < sec1.clmin for c in s.cldes[:s.ii])


def _set_grid(s: State) -> None:
    """CALL SETX: radial stations into the State."""
    g = setx(s.ii, s.xi0, s.xitip, s.ixspac)
    s.dt = g.dt
    s.t[:s.ii], s.xi[:s.ii], s.dxi[:s.ii], s.xv[:s.ii + 1] = g.t, g.xi, g.dxi, g.xv


def desini(s: State) -> None:
    """Port of SUBROUTINE DESINI: elliptic-like chord to get CD, Larrabee's circulation with
    Prandtl's tip factor, then (unless FAST) a Goldstein circulation from the potential
    or vortex formulation; chord from circulation and the design CL."""
    ii = s.ii
    secs = sections(s)
    adv = s.adv
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # ---- elliptical-like planform for the CD calculation
        for i in range(ii):
            xi = float(s.xi[i])
            xx = div(xi, adv)
            s.ch[i] = div(0.1 * xx * xx * fm.sqrt(s.xitip * s.xitip - xi * xi), 1.0 + xx * xx)

        # ---- CD array
        for i in range(ii):
            t = div(float(s.xi[i]), adv)
            w = fm.sqrt(1.0 + t * t)
            s.re[i] = div(float(s.ch[i]) * w * s.rho * s.vel * s.rad, s.rmu)
            a = getalf(i, float(s.cl[i]), w, xi=s.xi, iaero=s.iaero, sections=secs, vel=s.vel,
                       vso=s.vso, fix_fallback=s.fixes.section_fallback)
            s.stall[i] = a.stallf
            r = getclcdcm(i, a.alf, w, float(s.re[i]), xi=s.xi, iaero=s.iaero, sections=secs,
                          vel=s.vel, vso=s.vso, fix_fallback=s.fixes.section_fallback)
            s.cl[i], s.cd[i], s.stall[i] = r.clift, r.cdrag, r.stallf

        # ---- Larrabee: Prandtl's factor F for the tip (none for ducted designs)
        blds = float(s.nblds)
        ri1 = ri2 = rj1 = rj2 = 0.0
        sfac = fm.sqrt(1.0 + div(1.0, s.adw * s.adw))
        for i in range(ii):
            dx = float(s.dxw[i])
            xw = float(s.xw[i])
            yy = div(xw, adv)
            yyinv = div(1.0, 1.0 + yy * yy)
            dol = div(float(s.cd[i]), float(s.cl[i]))
            if s.duct:
                f = 1.0
            else:
                arg = fm.fmin(20.0, 0.5 * blds * (1.0 - div(xw, s.xwtip)) * sfac)
                ek = fm.exp(-arg)
                fk = fm.sqrt(1.0 - ek * ek)
                f = fm.atan2(fk, ek) * 2.0 / PI           # XROTOR.INC PI (truncated)
            s.gam[i] = f * yy * yy * yyinv
            gam = float(s.gam[i])
            ri1 = ri1 + 4.0 * xw * gam * (1.0 - div(dol, yy)) * dx
            ri2 = ri2 + 2.0 * xw * gam * (1.0 - div(dol, yy)) * dx * yyinv
            rj1 = rj1 + 4.0 * xw * gam * (1.0 + dol * yy) * dx
            rj2 = rj2 + 2.0 * xw * gam * (1.0 + dol * yy) * dx * yy * yy * yyinv

        # ---- displacement velocity
        if s.dest:
            disc = fm.fmax(0.0, 1.0 - div(8.0 / PI * s.tdes * ri2, ri1 * ri1))
            zeta = div(0.5 * ri1, ri2) * (1.0 - fm.sqrt(disc))
        if s.desp:
            disc = fm.fmax(0.0, 1.0 + div(8.0 / PI * s.pdes * rj2, rj1 * rj1))
            zeta = div(0.5 * rj1, rj2) * (fm.sqrt(disc) - 1.0)
        if not (s.dest or s.desp):
            raise ValueError("DESINI: neither thrust (DEST) nor power (DESP) is specified")

        s.gam[:ii] = s.gam[:ii] * div(2.0 * PI * adv, blds) * zeta

        # ---- better guess for the wake advance ratio
        s.effinv = 1.0 + 0.5 * zeta
        if s.free:
            s.adw = adv * s.effinv

        if not s.fast:
            # ---- Goldstein circulation: radially constant displacement velocity v'
            if not s.vrtx:
                helico(s)
            else:
                vrtxco(s)
            q = np.zeros((IQ, IQ))
            for i in range(ii):
                yy = div(float(s.xw[i]), adv)
                s.gam[i] = div(0.5 * zeta * yy, 1.0 + yy * yy)
                q[i, :ii] = s.vind_gam[2, i, :ii]
            gauss(q, s.gam, ii)

        # ---- circulation and chord
        for i in range(ii):
            t = div(float(s.xi[i]), adv)
            w = fm.sqrt(1.0 + t * t)
            s.ch[i] = div(2.0 * float(s.gam[i]), w * float(s.cl[i]))
    s.lrotor = True


def deiter(s: State, icase: int) -> None:
    """Port of SUBROUTINE DEITER(ICASE): Newton iteration on circulation, induced efficiency
    EFFINV and wake advance ratio, enforcing a radially constant induced efficiency on the
    equivalent prop and the design thrust (DEST) or power (DESP). Up to NITERD iterations."""
    if s.duct:
        raise NotImplementedError("DEITER: ducted design reads an uninitialized GAMDES in 7.55")
    ii = s.ii
    k1, k2 = ii, ii + 1                                 # Fortran K1 = II+1, K2 = II+2
    q = np.zeros((IQ, IQ))
    eps = 5.0e-07

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for it in range(1, max(s.niterd, 1) + 1):
            if s.free or it == 1:
                induced_model(s)
            vcalc(s)
            setxw(s)
            tpq(s, icase)

            s.gresmx = 0.0
            s.fresmx = 0.0
            s.aresmx = 0.0
            adv, adw = s.adv, s.adw

            # ---- constant induced efficiency at each station of the equivalent prop
            for i in range(ii):
                xi = float(s.xi[i])
                xw = float(s.xw[i])
                vw = float(s.vwak[i])
                vaw = div(vw * xw, adw)
                vaw_vw = div(xw, adw)
                vaw_xw = div(vw, adw)
                vaw_adw = -div(vaw, adw)
                utotw = s.urduct
                wa, wt = uvadd(s, xi)

                cw = div(xi, adv) - wt - vw
                cw_vw = -1.0
                cw_xw = 0.0
                sw = utotw + wa + vaw
                sw_adw = vaw_adw
                sw_vw = vaw_vw
                sw_xw = vaw_xw

                eiloc = div(sw, cw) * div(xi, adv)
                ei_sw = div(1.0, cw) * div(xi, adv)
                ei_cw = -div(eiloc, cw)
                ei_xw = 0.0
                ei_vw = ei_sw * sw_vw + ei_cw * cw_vw
                ei_xw = ei_sw * sw_xw + ei_cw * cw_xw + ei_xw
                ei_adw = ei_sw * sw_adw

                if ei_vw == 0.0:
                    q[i, :ii] = 0.0
                    s.dq[i] = 0.0
                    q[i, i] = 1.0
                    q[i, k1] = 0.0
                    q[i, k2] = 0.0
                else:
                    s.dq[i] = eiloc - s.effinv
                    q[i, :ii] = ei_xw * s.xw_gam[i, :ii] + ei_vw * s.vw_gam[i, :ii]
                    q[i, k1] = -1.0
                    q[i, k2] = ei_xw * float(s.xw_adw[i]) + ei_vw * float(s.vw_adw[i]) + ei_adw
                s.gresmx = fm.fmax(s.gresmx, abs(float(s.dq[i])))

            q[k1, :k2 + 1] = 0.0
            q[k2, :k2 + 1] = 0.0

            # ---- design thrust or power, defined with the equivalent prop
            if s.dest:
                s.dq[k1] = s.twak + s.tvis - s.tdes
                q[k1, :ii] = s.tw_gam[:ii] + s.tv_gam[:ii]
                q[k1, k2] = s.tw_adw + s.tv_adw
                s.fresmx = fm.fmax(s.fresmx, abs(div(float(s.dq[k1]), s.tdes)))
            elif s.desp:
                s.dq[k1] = s.pwak + s.pvis - s.pdes
                q[k1, :ii] = s.pw_gam[:ii] + s.pv_gam[:ii]
                q[k1, k2] = s.pw_adw + s.pv_adw
                s.fresmx = fm.fmax(s.fresmx, abs(div(float(s.dq[k1]), s.pdes)))

            if s.free:
                # ---- wake advance ratio from the equivalent prop's thrust and power
                s.adwfctr = 1.0
                advfact = 1.0
                s.dq[k2] = div(s.adwfctr * adw * s.twak, s.pwak) - adv * advfact
                z_tw = div(s.adwfctr * adw, s.pwak)
                z_pw = -div(s.adwfctr * adw * s.twak, s.pwak * s.pwak)
                q[k2, :ii] = z_tw * s.tw_gam[:ii] + z_pw * s.pw_gam[:ii]
                q[k2, k2] = z_tw * s.tw_adw + z_pw * s.pw_adw + div(s.adwfctr * s.twak, s.pwak)
                s.aresmx = fm.fmax(s.aresmx, abs(div(float(s.dq[k2]), adv)))
            else:
                s.dq[k2] = 0.0
                q[k2, k2] = 1.0

            # ---- solve the linearized Newton system
            gauss(q, s.dq, k2 + 1)

            # ---- under-relaxation against excessive circulation changes
            rlx = 1.0
            if it <= 2:
                rlx = 0.2
            if s.niterd == 0:
                rlx = 0.0
            for i in range(ii):
                s.dgam[i] = -s.dq[i]
                dgam, gam = float(s.dgam[i]), float(s.gam[i])
                if div(rlx * dgam, gam) < -0.3:
                    rlx = div(-0.3 * gam, dgam)
                if div(rlx * dgam, gam) > 0.3:
                    rlx = div(0.3 * gam, dgam)
            s.deff = -float(s.dq[k1])
            s.dadw = -float(s.dq[k2])
            if rlx * s.dadw > 0.5 * adw:
                rlx = fm.fmin(rlx, div(0.5 * adw, s.dadw))
            if rlx * s.dadw < -0.3 * adw:
                rlx = fm.fmin(rlx, div(-0.3 * adw, s.dadw))
            s.rlx = rlx

            # ---- update circulation, induced efficiency, wake advance ratio
            rms = 0.0
            for i in range(ii):
                dgam = float(s.dgam[i])
                s.gam[i] = s.gam[i] + rlx * dgam
                rms = rms + div(dgam * dgam, 1.0 + div(1.0, adv * adv))
            s.rms = fm.sqrt(div(rms, float(ii)))
            s.effinv = s.effinv + rlx * s.deff
            s.adw = s.adw + rlx * s.dadw

            if rlx < 0.2:
                filter_(s.gam, 0.2 * ii, ii)

            # ---- convergence
            if fm.fmax(s.rms, abs(s.dadw)) <= eps and abs(s.deff) <= eps * s.effinv:
                if s.free:
                    induced_model(s)
                vcalc(s)
                tpq(s, icase)
                s.conv = True
                return


def desmil(s: State, icase: int) -> None:
    """Port of SUBROUTINE DESMIL(ICASE). (All three formulations call DEITER.)"""
    deiter(s, icase)


def desgen(s: State, betz_ok: bool = False, refresh_sections: bool | None = None,
           fix_linear_cl: bool | None = None) -> bool:
    """Port of SUBROUTINE DESGEN: design an MIL rotor from the DESI parameters in the State
    (NBLDS, RADDES, R0DES, RWDES, VELDES, ADVDES or RPMDES, TDDES/DEST or PDDES/DESP, CLDES).

    Returns False where DESGEN returns early without designing: design CL past section 1's
    stall limits, or (unless betz_ok, the answer to XROTOR's prompt) thrust or power beyond
    the Betz windmill limits.

    refresh_sections (default s.fixes.design_sections): False reproduces 7.55, which keeps
    IAERO from the *previous* radial grid (stale section indices, see CLAUDE.md). True calls
    SETIAERO on the new grid.

    fix_linear_cl (default s.fixes.linear_cl): True re-evaluates the last SETCLD's linear CL
    (s.cldes_linear) on the new grid, and moves the stall-limit check after SETX so that it
    tests those values; a design CL past the limits then returns False after the operating
    point and grid have been set up.
    """
    if refresh_sections is None:
        refresh_sections = s.fixes.design_sections
    if fix_linear_cl is None:
        fix_linear_cl = s.fixes.linear_cl
    fix_cl = fix_linear_cl and s.cldes_linear is not None
    ii = s.ii
    if not fix_cl and _past_stall(s):
        return False                                   # 'Design CL is past stall limit'

    # ---- operating parameters
    s.conv = False
    s.rad = s.raddes
    s.vel = s.veldes
    if s.advdes == 0.0:
        s.adv = div(s.veldes, s.raddes * s.rpmdes * PI / 30.0)
    else:
        s.adv = s.advdes

    # ---- reasonable thrust or power inputs (Betz limits)
    v2, r2 = s.vel * s.vel, s.rad * s.rad
    tmin = -0.5 * PI * s.rho * v2 * r2
    pmin = -0.5 * PI * s.rho * (v2 * s.vel) * r2 * 16.0 / 27.0
    if (s.dest and s.tddes < tmin) or (s.desp and s.pddes < pmin):
        if not betz_ok:
            return False

    # ---- non-dimensional design thrust or power
    s.tdes = div(s.tddes, s.rho * v2 * r2)
    s.pdes = div(s.pddes, s.rho * (v2 * s.vel) * r2)

    # ---- discrete coordinates
    s.xi0 = div(s.r0des, s.rad)
    s.xw0 = div(s.rwdes, s.rad)
    _set_grid(s)
    xwinit(s)
    if refresh_sections:
        from .aero import setiaero
        s.iaero[:ii] = setiaero(s.xi[:ii], sections(s))
    if fix_cl:
        cl1, cl2 = s.cldes_linear
        s.cldes[:ii] = (1.0 - s.xi[:ii]) * cl1 + s.xi[:ii] * cl2
        if _past_stall(s):
            return False

    # ---- CL array from the design CL
    s.cl[:ii] = s.cldes[:ii]
    s.stall[:ii] = False

    # ---- wake advance ratios
    s.effinv = 1.0
    if s.free and s.dest:
        s.effinv = 0.50 * fm.sqrt(abs(1.0 + s.tdes * 2.0 / PI)) + 0.50
    if s.free and s.desp:
        s.effinv = 0.25 * fm.sqrt(abs(1.0 + s.pdes * 4.0 / PI)) + 0.75
    s.adw = s.effinv * s.adv

    # ---- initialize the rotor, then the MIL design (fixed CL, varying chord and twist)
    if s.ldesini or not s.lrotor:
        desini(s)
    desmil(s, 2)
    return True
