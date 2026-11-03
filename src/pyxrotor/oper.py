"""Operating-point solver (xoper.f).

Ported routines:
    TPQ     tpq()      thrust, torque, power and their sensitivities (real and equivalent prop)
    APINIT  apinit()   initial circulation and wake advance ratio (graded momentum)
    APITER  apiter()   Newton iteration for an operating point
    APER    aper()     APINIT (optional) + APITER

ISPEC (target): 1 thrust TSPEC, 2 torque QSPEC, 3 power PSPEC, 4 fixed advance ratio,
                5 power from the engine power-rpm line.
ICON (constraint): 1 advance ratio (rpm) fixed, 2 blade pitch fixed.

All state lives in the State; outputs are the COMMON variables XROTOR sets (TTOT, PTOT,
QTOT, GAM, CL, CD, ..., CONV). Faithfully reproduced 7.55 quirks are marked "7.55:".
"""
from __future__ import annotations

import math

import numpy as np

from . import fmath as fm
from .aero import AERODATA_ROW, AeroSection, getalf, getclcdcm
from .constants import IQ, PI
from .fmath import div
from .induced import gradmo, helico, vcalc, vrtxco
from .numerics import deval, filter_, gauss, seval
from .state import State
from .wake import cscalc, setxw, uvadd


# --------------------------------------------------------------------------- helpers
def sections(s: State) -> list[AeroSection]:
    """The aero sections stored in XIAERO/AERODATA/NAERO."""
    out = []
    for n in range(s.naero):
        kw = {name: float(s.aerodata[row - 1, n]) for name, row in AERODATA_ROW.items()}
        out.append(AeroSection(xisect=float(s.xiaero[n]), **kw))
    return out


def _station_aero(s: State, secs, i: int, alfa: float, w: float, rey: float):
    """CALL GETCLCDCM(I,ALFA,W,REY, CL(I),..,STALL(I), CD(I),.., CM(I),..): stores CL, CD,
    CM, STALL at station i and returns the full result."""
    r = getclcdcm(i, alfa, w, rey, xi=s.xi, iaero=s.iaero, sections=secs, vel=s.vel, vso=s.vso,
                  fix_fallback=s.fixes.section_fallback)
    s.cl[i], s.cd[i], s.cm[i], s.stall[i] = r.clift, r.cdrag, r.cmom, r.stallf
    return r


def _rey(s: State, ch: float, w: float) -> float:
    """REY = CH(I)*ABS(W) * RHO*VEL*RAD/RMU, evaluated left to right."""
    return div(ch * abs(w) * s.rho * s.vel * s.rad, s.rmu)


def induced_model(s: State) -> None:
    """The APITER/DESI selection: GRADMO if FAST, else HELICO, or VRTXCO if VRTX. Sets IWTYP."""
    if s.fast:
        gradmo(s)
        s.iwtyp = 1
    elif not s.vrtx:
        helico(s)
        s.iwtyp = 2
    else:
        vrtxco(s)
        s.iwtyp = 3


# ------------------------------------------------------------------------------- TPQ
def tpq(s: State, itype: int) -> None:
    """Port of SUBROUTINE TPQ(ITYPE): thrust, power and torque, and their sensitivities to
    GAM(J), ADV, ADW and the blade-angle change.

    ITYPE 1: analysis (fixed blade angles). 2: design, fixed CL, chord from circulation.
    3: design, fixed chord, angles from CL.

    Not ported: momentum-theory checks (TMOM, PMOM, VAAavg, VATavg, TDIM, PDIM) that TPQ
    computes into unused locals, and the unused equivalent-prop axial velocity SW.
    """
    ii = s.ii
    secs = sections(s)
    s.tinv = s.pinv = 0.0
    s.twak = s.pwak = 0.0
    s.tvis = s.pvis = 0.0
    s.ti_adv = s.pi_adv = s.ti_adw = s.pi_adw = 0.0
    s.tw_adv = s.pw_adv = s.tw_adw = s.pw_adw = 0.0
    s.tv_adv = s.pv_adv = s.tv_adw = s.pv_adw = 0.0
    s.tv_dbe = s.pv_dbe = 0.0
    for arr in (s.ti_gam, s.pi_gam, s.tw_gam, s.pw_gam, s.tv_gam, s.pv_gam):
        arr[:ii] = 0.0

    cosr = math.cos(s.rake)
    blds = float(s.nblds)
    adv, adw = s.adv, s.adw
    rfac = (s.rho, s.vel, s.rad, s.rmu)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for i in range(ii):
            xi = float(s.xi[i])
            bdx = blds * float(s.dxi[i])
            xx = div(xi, adv)
            xx_adv = -div(xx, adv)

            v = cscalc(s, i)
            w, ci, si, phi = v.w, v.ci, v.si, v.phi
            alfa = float(s.beta[i]) - phi

            if itype == 1:
                # ---- analysis: fixed local beta (except for pitch change)
                alfa = float(s.beta[i]) - phi
                al_gi = 0.0
                al_dbe = 1.0
                al_adv = -v.p_adv
                al_vt = -v.p_vt
                al_va = -v.p_va
                rey = _rey(s, float(s.ch[i]), w)
                r = _station_aero(s, secs, i, alfa, w, rey)
                cl_al, cl_w = r.cl_alf, r.cl_w
                cl_gi = cl_al * al_gi
                cl_dbe = cl_al * al_dbe
                cl_adv = cl_al * al_adv + cl_w * v.w_adv
                cl_vt = cl_al * al_vt + cl_w * v.w_vt
                cl_va = cl_al * al_va + cl_w * v.w_va
                ch_gi = ch_adv = ch_vt = ch_va = 0.0
            elif itype == 2:
                # ---- design: fixed local CL, chord from circulation
                a = getalf(i, float(s.cl[i]), w, xi=s.xi, iaero=s.iaero, sections=secs,
                           vel=s.vel, vso=s.vso, fix_fallback=s.fixes.section_fallback)
                alfa, s.stall[i] = a.alf, a.stallf
                al_gi = 0.0
                al_dbe = 0.0
                al_adv = a.alf_w * v.w_adv
                al_vt = a.alf_w * v.w_vt
                al_va = a.alf_w * v.w_va
                cl_gi = cl_dbe = cl_adv = cl_vt = cl_va = 0.0
                gam = float(s.gam[i])
                cli = float(s.cl[i])
                chnew = div(2.0 * gam, w * cli)
                if chnew <= 0.0:
                    # ---- chord going zero or negative: use nearby station data this pass
                    if i == 0:
                        s.ch[i] = s.ch[i + 1]
                    elif i == ii - 1:
                        s.ch[i] = s.ch[i - 1]
                    else:
                        s.ch[i] = 0.5 * (s.ch[i - 1] + s.ch[i + 1])
                    ch_gi = ch_adv = ch_vt = ch_va = 0.0
                else:
                    s.ch[i] = div(2.0 * gam, w * cli)
                    ch_gi = div(2.0, w * cli)
                    chw = -div(float(s.ch[i]), w)
                    ch_adv = chw * v.w_adv
                    ch_vt = chw * v.w_vt
                    ch_va = chw * v.w_va
                s.beta[i] = alfa + phi
                s.beta0[i] = s.beta[i]
            elif itype == 3:
                # ---- design: fixed local chord, angles from CL
                ch = float(s.ch[i])
                s.cl[i] = div(2.0 * float(s.gam[i]), w * ch)
                cl_gi = div(2.0, w * ch)
                cl_dbe = 0.0
                clw = -div(float(s.cl[i]), w)
                cl_adv = clw * v.w_adv
                cl_vt = clw * v.w_vt
                cl_va = clw * v.w_va
                a = getalf(i, float(s.cl[i]), w, xi=s.xi, iaero=s.iaero, sections=secs,
                           vel=s.vel, vso=s.vso, fix_fallback=s.fixes.section_fallback)
                alfa, s.stall[i] = a.alf, a.stallf
                al_gi = a.alf_cl * cl_gi
                al_dbe = a.alf_cl * cl_dbe
                al_adv = a.alf_cl * cl_adv + a.alf_w * v.w_adv
                al_vt = a.alf_cl * cl_vt + a.alf_w * v.w_vt
                al_va = a.alf_cl * cl_va + a.alf_w * v.w_va
                ch_gi = ch_adv = ch_vt = ch_va = 0.0
                s.beta[i] = alfa + phi
                s.beta0[i] = s.beta[i]
            else:
                raise ValueError(f"TPQ: ITYPE must be 1, 2 or 3, not {itype}")

            ch = float(s.ch[i])
            s.re[i] = div(ch * abs(w) * s.rho * s.vel * s.rad, s.rmu)
            re_w = div(ch * s.rho * s.vel * s.rad, s.rmu)
            re_ch = div(abs(w) * s.rho * s.vel * s.rad, s.rmu)
            re_gi = re_ch * ch_gi
            re_adv = re_ch * ch_adv + re_w * v.w_adv
            re_vt = re_ch * ch_vt + re_w * v.w_vt
            re_va = re_ch * ch_va + re_w * v.w_va

            # ---- CD (and CM) with sensitivities
            r = _station_aero(s, secs, i, alfa, w, float(s.re[i]))
            cd_al, cd_w, cd_re = r.cd_alf, r.cd_w, r.cd_rey
            cd_gi = cd_al * al_gi + cd_re * re_gi
            cd_adv = cd_al * al_adv + cd_re * re_adv + cd_w * v.w_adv
            cd_vt = cd_al * al_vt + cd_re * re_vt + cd_w * v.w_vt
            cd_va = cd_al * al_va + cd_re * re_va + cd_w * v.w_va
            cd_dbe = cd_al * al_dbe

            cli, cdi = float(s.cl[i]), float(s.cd[i])
            # ---- local efficiencies (corrected for rake)
            eff = div(div(cli * ci - cdi * si, cdi * ci + cli * si), xx)
            eff = eff * cosr
            effi = div(ci, si * xx)
            effi = effi * cosr
            s.effp[i] = div(eff, effi)

            hwc = 0.5 * w * ch
            hwc_w = 0.5 * ch
            hwc_ch = 0.5 * w

            # ---- viscous thrust and power on the real prop
            cosrv = 1.0
            dtv = -hwc * cdi * si * bdx * cosrv
            dtv_cd = -hwc * si * bdx * cosrv
            dtv_si = -hwc * cdi * bdx * cosrv
            dtv_w = -hwc_w * cdi * si * bdx * cosrv
            dtv_ch = -hwc_ch * cdi * si * bdx * cosrv
            dtv_gi = dtv_cd * cd_gi + dtv_ch * ch_gi
            dtv_dbe = dtv_cd * cd_dbe
            dtv_adv = dtv_cd * cd_adv + dtv_ch * ch_adv + dtv_w * v.w_adv
            dtv_vt = dtv_cd * cd_vt + dtv_ch * ch_vt + dtv_w * v.w_vt
            dtv_va = dtv_cd * cd_va + dtv_ch * ch_va + dtv_si * v.si_va + dtv_w * v.w_va

            s.tvis = s.tvis + dtv
            s.tv_adv = s.tv_adv + dtv_adv
            s.tv_adw = dtv_vt * v.vt_adw + dtv_va * v.va_adw     # 7.55: not accumulated
            s.tv_dbe = s.tv_dbe + dtv_dbe
            s.tv_gam[i] = s.tv_gam[i] + dtv_gi
            s.tv_gam[:ii] = s.tv_gam[:ii] + dtv_vt * s.vind_gam[2, i, :ii] + dtv_va * s.vind_gam[0, i, :ii]

            dpv = hwc * cdi * ci * bdx * xx
            dpv_cd = hwc * ci * bdx * xx
            dpv_ci = hwc * cdi * bdx * xx
            dpv_w = hwc_w * cdi * ci * bdx * xx
            dpv_ch = hwc_ch * cdi * ci * bdx * xx
            dpv_gi = dpv_cd * cd_gi + dpv_ch * ch_gi
            dpv_dbe = dpv_cd * cd_dbe
            dpv_adv = (dpv_cd * cd_adv + dpv_ch * ch_adv + dpv_ci * v.ci_adv + dpv_w * v.w_adv
                       + hwc * cdi * ci * bdx * xx_adv)
            dpv_vt = dpv_cd * cd_vt + dpv_ch * ch_vt + dpv_ci * v.ci_vt + dpv_w * v.w_vt
            dpv_va = dpv_cd * cd_va + dpv_ch * ch_va + dpv_w * v.w_va

            s.pvis = s.pvis + dpv
            s.pv_adv = s.pv_adv + dpv_adv
            s.pv_adw = dpv_vt * v.vt_adw + dpv_va * v.va_adw     # 7.55: not accumulated
            s.pv_dbe = s.pv_dbe + dpv_dbe
            s.pv_gam[i] = s.pv_gam[i] + dpv_gi
            s.pv_gam[:ii] = s.pv_gam[:ii] + dpv_vt * s.vind_gam[2, i, :ii] + dpv_va * s.vind_gam[0, i, :ii]

            # ---- inviscid thrust and power on the real prop
            gam = float(s.gam[i])
            cosri = 1.0
            dti = gam * ci * bdx * cosri
            dti_ci = gam * bdx * cosri
            dti_gi = ci * bdx * cosri
            dti_vt = dti_ci * v.ci_vt
            dti_adv = dti_ci * v.ci_adv
            dti_adw = dti_vt * v.vt_adw
            s.tinv = s.tinv + dti
            s.ti_adv = s.ti_adv + dti_adv
            s.ti_adw = s.ti_adw + dti_adw
            s.ti_gam[i] = s.ti_gam[i] + dti_gi
            s.ti_gam[:ii] = s.ti_gam[:ii] + dti_vt * s.vind_gam[2, i, :ii]

            dpi = gam * si * bdx * xx
            dpi_si = gam * bdx * xx
            dpi_gi = si * bdx * xx
            dpi_xx = gam * si * bdx
            dpi_va = dpi_si * v.si_va
            dpi_adv = dpi_xx * xx_adv
            dpi_adw = dpi_va * v.va_adw
            s.pinv = s.pinv + dpi
            s.pi_adv = s.pi_adv + dpi_adv
            s.pi_adw = s.pi_adw + dpi_adw
            s.pi_gam[i] = s.pi_gam[i] + dpi_gi
            s.pi_gam[:ii] = s.pi_gam[:ii] + dpi_va * s.vind_gam[0, i, :ii]

            # ---- inviscid thrust and power on the equivalent prop (same omega and GAM)
            vw = float(s.vwak[i])
            wa, wt = uvadd(s, xi)
            xw = float(s.xw[i])
            cw = div(xw, adv) - wt - vw
            cw_adv = -div(xw, adv * adv)
            cw_vw = -1.0

            dxw = float(s.dxw[i])
            dtw = gam * cw * blds * dxw
            dtw_gi = cw * blds * dxw
            dtw_cw = gam * blds * dxw
            dtw_dxw = gam * cw * blds
            dtw_vw = dtw_cw * cw_vw
            dtw_adv = dtw_cw * cw_adv + dtw_vw * float(s.vw_adv[i]) + dtw_dxw * float(s.dxw_adv[i])
            dtw_adw = dtw_vw * float(s.vw_adw[i]) + dtw_dxw * float(s.dxw_adw[i])
            s.twak = s.twak + dtw
            s.tw_adv = s.tw_adv + dtw_adv
            s.tw_adw = s.tw_adw + dtw_adw
            s.tw_gam[i] = s.tw_gam[i] + dtw_gi
            s.tw_gam[:ii] = s.tw_gam[:ii] + dtw_vw * s.vw_gam[i, :ii] + dtw_dxw * s.dxw_gam[i, :ii]

            dpw = div(gam * si * bdx * xi, adv)
            dpw_si = div(gam * bdx * xi, adv)
            dpw_gi = div(si * bdx * xi, adv)
            dpw_adv = -div(dpw, adv)
            dpw_va = dpw_si * v.si_va
            dpw_adw = dpw_va * v.va_adw
            s.pwak = s.pwak + dpw
            s.pw_adv = s.pw_adv + dpw_adv
            s.pw_adw = s.pw_adw + dpw_adw
            s.pw_gam[i] = s.pw_gam[i] + dpw_gi
            s.pw_gam[:ii] = s.pw_gam[:ii] + dpw_va * s.vind_gam[0, i, :ii]

            # ---- blade thrust and power contributions (per blade, per span)
            s.dtii[i] = div(dti, bdx)
            s.dpii[i] = div(dpi, bdx)
            s.dtwi[i] = div(dtw, bdx)
            s.dpwi[i] = div(dpw, bdx)
            s.dtvi[i] = div(dtv, bdx)
            s.dpvi[i] = div(dpv, bdx)

    s.ttot = s.twak + s.tvis
    s.ptot = s.pwak + s.pvis
    s.qtot = s.ptot * s.adv


# ---------------------------------------------------------------------------- APINIT
def apinit(s: State) -> None:
    """Port of SUBROUTINE APINIT: circulation from the blade-element lift with no induced
    velocity, then up to 10 graded-momentum passes converging the wake advance ratio."""
    ii = s.ii
    secs = sections(s)
    blds = float(s.nblds)
    s.dbeta = 0.0
    uduct = 0.0
    vaduct_va = 1.0
    if s.duct:
        uduct = s.urduct - 1.0
        vaduct_va = 2.0 * s.urduct
    s.adw = s.adv * (1.0 + uduct)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # ==== section circulation neglecting induced velocity
        tsum = 0.0
        for i in range(ii):
            xi = float(s.xi[i])
            utot = s.urduct + float(s.ubody[i])
            wa, wt = uvadd(s, xi)
            si = utot + wa
            ci = div(xi, s.adv) - wt
            w = fm.sqrt(ci * ci + si * si)
            phi = fm.atan2(si, ci)
            alfa = float(s.beta[i]) - phi
            rey = _rey(s, float(s.ch[i]), w)
            _station_aero(s, secs, i, alfa, w, rey)
            s.gam[i] = 0.5 * float(s.cl[i]) * w * float(s.ch[i])
            tsum = tsum + blds * float(s.gam[i]) * ci * float(s.dxi[i])

        # ---- momentum-theory estimate of axial velocity sets the wake advance ratio
        vhsq = 0.5 * tsum / PI
        vhsq = fm.fmax2(vhsq, -0.25)                      # NaN -> -0.25 (minmax_sites)
        s.adw = s.adv * 0.5 * (1.0 + fm.sqrt(1.0 + 4.0 * vhsq))

        vcalc(s)

        # ==== refine with graded momentum
        for _ in range(10):                               # NITERG
            gradmo(s)
            tsum = 0.0
            t_adw = 0.0
            dclmax = 0.0
            rlxmin = 1.0
            for i in range(ii):
                xi = float(s.xi[i])
                gam = float(s.gam[i])
                vt_gam = float(s.vind_gam[2, i, i])
                vt = vt_gam * gam
                vt_adw = float(s.vind_adw[2, i])
                va_gam = float(s.vind_gam[0, i, i])
                va = va_gam * gam
                va_adw = float(s.vind_adw[0, i])

                uduct = 0.0
                vaduct_va = 1.0
                if s.duct:
                    uduct = s.urduct - 1.0
                    vaduct_va = 2.0 * s.urduct
                utot = 1.0 + uduct + float(s.ubody[i])
                wa, wt = uvadd(s, xi)

                ci = div(xi, s.adv) - wt - vt
                ci_vt = -1.0
                si = utot + wa + va * vaduct_va
                si_va = vaduct_va
                wsq = ci * ci + si * si
                w = fm.sqrt(wsq)
                w_vt = div(ci * ci_vt, w)
                w_va = div(si * si_va, w)
                phi = fm.atan2(si, ci)
                p_vt = div(-(si * ci_vt), wsq)
                p_va = div(ci * si_va, wsq)

                alfa = float(s.beta[i]) - phi
                al_vt = -p_vt
                al_va = -p_va

                ch = float(s.ch[i])
                rey = _rey(s, ch, w)
                r = _station_aero(s, secs, i, alfa, w, rey)

                # ---- Res( CL(AL W), W, GAM )
                rez = ch * float(s.cl[i]) * w - 2.0 * gam
                z_cl = ch * w
                z_w = ch * float(s.cl[i])
                z_g = -2.0
                z_al = z_cl * r.cl_alf
                z_w = z_cl * r.cl_w + z_w
                z_vt = z_w * w_vt + z_al * al_vt
                z_va = z_w * w_va + z_al * al_va
                z_adw = z_vt * vt_adw + z_va * va_adw
                z_g = z_vt * vt_gam + z_va * va_gam + z_g

                delg = -div(rez, z_g)
                dcl = div(2.0 * delg, ch * w)

                # ---- limit the GAM update by the CL change (RLX is COMMON)
                rlx = 1.0
                if rlx * abs(dcl) > 0.2:
                    if dcl != 0.0:
                        rlx = fm.fmin(rlx, div(0.2, abs(dcl)))
                s.rlx = rlx
                if abs(dcl) > abs(dclmax):
                    dclmax = dcl
                if abs(rlx) < rlxmin:
                    rlxmin = rlx

                s.gam[i] = gam + rlx * delg
                g_adw = -div(z_adw, z_g)

                cosr = 1.0
                gam = float(s.gam[i])
                tsum = tsum + blds * gam * ci * float(s.dxi[i]) * cosr
                t_g = blds * ci * float(s.dxi[i]) * cosr
                t_vt = blds * gam * ci_vt * float(s.dxi[i]) * cosr
                t_adw = t_adw + (t_g + t_vt * vt_gam) * g_adw + t_vt * vt_adw

            # ---- momentum-theory estimate of induced axial velocity
            vhsq = 0.5 * tsum / PI
            vhsq = fm.fmax2(vhsq, -0.2499)                # NaN -> -0.2499
            vhsq_t = 0.5 / PI
            root = fm.sqrt(1.0 + 4.0 * vhsq)
            rez = s.adw - s.adv * 0.5 * (1.0 + root)
            z_adw = 1.0 - div(s.adv, root) * vhsq_t * t_adw
            s.dadw = -div(rez, z_adw)
            s.dadw = fm.fmin2(s.dadw, 10.0 * s.adw)       # NaN -> the ADW bound
            s.dadw = fm.fmax2(s.dadw, -0.9 * s.adw)
            s.adw = s.adw + s.dadw

            if rlxmin < 0.2:
                filter_(s.gam, 0.2 * ii, ii)
            if abs(dclmax) < 0.001:
                break


# ---------------------------------------------------------------------------- APITER
def apiter(s: State, ispec: int, icon: int) -> None:
    """Port of SUBROUTINE APITER(ISPEC,ICON): Newton iteration on circulation, ADV, ADW and
    blade-angle change, up to NITERA iterations, to RMS(dGAM) <= 1e-7. Sets CONV on success.

    Not ported: the 0.75R "ADVFACT" correction, which 7.55 computes and then overwrites
    with 1.0 (the computation has no side effects).
    """
    ii = s.ii
    k1, k2, k3 = ii, ii + 1, ii + 2                   # Fortran K1, K2, K3 (0-based rows)
    clmax = np.zeros(ii)
    clmin = np.zeros(ii)
    dclstall = np.zeros(ii)
    q = np.zeros((IQ, IQ))                            # local Q(IQ,IQ)
    eps = 1.0e-07
    advfact = 1.0
    rhov2r2 = s.rho * (s.vel * s.vel) * (s.rad * s.rad)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for it in range(1, max(s.nitera, 1) + 1):
            if s.free or it == 1:
                induced_model(s)
            vcalc(s)
            setxw(s)
            tpq(s, 1)
            secs = sections(s)

            s.gresmx = 0.0
            s.fresmx = 0.0
            s.aresmx = 0.0
            q[k2, :k3 + 1] = 0.0

            if s.free:
                # ---- wake advance ratio from the equivalent prop's thrust and power
                s.dq[k2] = div(s.adwfctr * s.adw * s.twak, s.pwak) - s.adv * advfact
                z_tw = div(s.adwfctr * s.adw, s.pwak)
                z_pw = -div(s.adwfctr * s.adw * s.twak, s.pwak * s.pwak)
                q[k2, :ii] = z_tw * s.tw_gam[:ii] + z_pw * s.pw_gam[:ii]
                q[k2, k1] = z_tw * s.tw_adv + z_pw * s.pw_adv - advfact
                q[k2, k2] = z_tw * s.tw_adw + z_pw * s.pw_adw + div(s.adwfctr * s.twak, s.pwak)
                s.aresmx = fm.fmax(s.aresmx, abs(div(float(s.dq[k2]), s.adv)))
            else:
                s.dq[k2] = 0.0
                q[k2, k2] = 1.0

            # ---- Gamma-CL relation at every station of the real prop
            w = 0.0
            for i in range(ii):
                v = cscalc(s, i)
                w = v.w
                alfa = float(s.beta[i]) - v.phi
                al_dbe = 1.0
                al_p = -1.0
                ch = float(s.ch[i])
                rey = _rey(s, ch, w)
                r = _station_aero(s, secs, i, alfa, w, rey)
                clmax[i], clmin[i], dclstall[i] = r.clmax, r.clmin, r.dcl_stall
                cli = float(s.cl[i])

                s.dq[i] = ch * cli * w - 2.0 * float(s.gam[i])
                z_cl = ch * w
                z_w = ch * cli
                z_gi = -2.0
                z_vt = z_cl * (r.cl_alf * al_p * v.p_vt + r.cl_w * v.w_vt) + z_w * v.w_vt
                z_va = z_cl * (r.cl_alf * al_p * v.p_va + r.cl_w * v.w_va) + z_w * v.w_va
                z_adv = z_cl * (r.cl_alf * al_p * v.p_adv + r.cl_w * v.w_adv) + z_w * v.w_adv
                z_dbe = z_cl * (r.cl_alf * al_dbe)

                q[i, :ii] = z_vt * s.vind_gam[2, i, :ii] + z_va * s.vind_gam[0, i, :ii]
                q[i, i] = q[i, i] + z_gi
                q[i, k1] = z_adv
                q[i, k2] = z_vt * v.vt_adw + z_va * v.va_adw
                q[i, k3] = z_dbe
                s.gresmx = fm.fmax(s.gresmx, abs(div(float(s.dq[i]), 0.1 * w)))

            # ---- target equation
            if ispec == 1:
                t_spec = div(s.tspec, rhov2r2)
                s.dq[k1] = s.twak + s.tvis - t_spec
                q[k1, :ii] = s.tw_gam[:ii] + s.tv_gam[:ii]
                q[k1, k1] = s.tw_adv + s.tv_adv
                q[k1, k2] = s.tw_adw + s.tv_adw
                q[k1, k3] = s.tv_dbe
                s.fresmx = fm.fmax(s.fresmx, abs(float(s.dq[k1])))
            elif ispec == 2:
                q_spec = div(s.qspec, s.rho * (s.vel * s.vel) * (s.rad * s.rad * s.rad))
                s.dq[k1] = (s.pwak + s.pvis) * s.adv - q_spec
                q[k1, :ii] = (s.pw_gam[:ii] + s.pv_gam[:ii]) * s.adv
                q[k1, k1] = (s.pw_adv + s.pv_adv) * s.adv + s.pwak + s.pvis
                q[k1, k2] = (s.pw_adw + s.pv_adw) * s.adv
                q[k1, k3] = s.pv_dbe * s.adv
                s.fresmx = fm.fmax(s.fresmx, abs(float(s.dq[k1])))
            elif ispec == 3:
                p_spec = div(s.pspec, s.rho * (s.vel * s.vel * s.vel) * (s.rad * s.rad))
                s.dq[k1] = s.pwak + s.pvis - p_spec
                q[k1, :ii] = s.pw_gam[:ii] + s.pv_gam[:ii]
                q[k1, k1] = s.pw_adv + s.pv_adv
                q[k1, k2] = s.pw_adw + s.pv_adw
                q[k1, k3] = s.pv_dbe
                s.fresmx = fm.fmax(s.fresmx, abs(float(s.dq[k1])))
            elif ispec == 4:
                s.dq[k1] = 0.0
                q[k1, :ii] = 0.0
                q[k1, k1] = 1.0
                q[k1, k2] = 0.0
                q[k1, k3] = 0.0
            elif ispec == 5:
                rhov3r2 = s.rho * (s.vel * s.vel * s.vel) * (s.rad * s.rad)
                p_spec = div(s.pspec, rhov3r2)
                p_spec_adv = 0.0
                if s.lpwrvar:
                    rpm = div(s.vel, s.rad * s.adv * PI / 30.0)
                    rpm_adv = -div(rpm, s.adv)
                    n = s.npwrvar
                    s.pspec = seval(rpm, s.pwrvar, s.xpwrvar, s.rpmvar, n)
                    pspec_rpm = deval(rpm, s.pwrvar, s.xpwrvar, s.rpmvar, n)
                    pspec_adv = pspec_rpm * rpm_adv
                    p_spec = div(s.pspec, rhov3r2)
                    p_spec_adv = div(pspec_adv, rhov3r2)
                s.dq[k1] = s.pwak + s.pvis - p_spec
                q[k1, :ii] = s.pw_gam[:ii] + s.pv_gam[:ii]
                q[k1, k1] = s.pw_adv + s.pv_adv - p_spec_adv
                q[k1, k2] = s.pw_adw + s.pv_adw
                q[k1, k3] = s.pv_dbe
                s.fresmx = fm.fmax(s.fresmx, abs(float(s.dq[k1])))
            else:
                raise ValueError(f"APITER: ISPEC must be 1-5, not {ispec}")

            # ---- constraint
            s.dq[k3] = 0.0
            q[k3, :k3 + 1] = 0.0
            if icon == 1:
                q[k3, k1] = 1.0
            if icon == 2:
                q[k3, k3] = 1.0

            # ---- solve the linearized Newton system
            gauss(q, s.dq, k3 + 1)

            rlx = 1.0
            if it <= 2:
                rlx = 0.2
            # ---- limiters on the Newton update
            for i in range(ii):
                s.dgam[i] = -s.dq[i]
                dgam = float(s.dgam[i])
                cli = float(s.cl[i])
                # 7.55: W here is the last station's W from the loop above, for every I
                dcl = div(2.0 * dgam, float(s.ch[i]) * w)
                dclmin = fm.fmax(1.5 * dclstall[i], abs(cli - clmin[i]))
                dclmax = fm.fmax(1.5 * dclstall[i], abs(clmax[i] - cli))
                # MIN(0.5,DCLMIN,DCLMAX) compiles as MINSD(MINSD(DCLMAX,DCLMIN),0.5)
                dcllim = fm.fmin2(fm.fmin2(dclmax, dclmin), 0.5)
                dcllim = fm.fmax2(dcllim, 0.01)
                if rlx * abs(dcl) > dcllim:
                    rlx = fm.fmin2(rlx, div(dcllim, abs(dcl)))
                # ---- limit GAM changes that change sign
                if dgam * float(s.dgamold[i]) < 0.0:
                    if abs(dgam) > 0.2 * abs(float(s.dgamold[i])):
                        rlx = fm.fmin2(rlx, 0.2)

            s.dadv = -float(s.dq[k1])
            s.dadw = -float(s.dq[k2])
            s.dbet = -float(s.dq[k3])

            if s.nitera == 0:
                rlx = 0.0

            # ---- limit blade-angle change to 0.05 rad, and advance-ratio changes
            dbet, dadv, dadw = s.dbet, s.dadv, s.dadw
            if rlx * dbet > 0.05:
                rlx = fm.fmin(rlx, div(0.05, dbet))
            if rlx * dbet < -0.05:
                rlx = fm.fmin(rlx, div(-0.05, dbet))
            if rlx * dadv > 0.5 * s.adv:
                rlx = fm.fmin(rlx, div(0.5 * s.adv, dadv))
            if rlx * dadv < -0.3 * s.adv:
                rlx = fm.fmin(rlx, div(-0.3 * s.adv, dadv))
            if rlx * dadw > 0.5 * s.adw:
                rlx = fm.fmin(rlx, div(0.5 * s.adw, dadw))
            if rlx * dadw < -0.3 * s.adw:
                rlx = fm.fmin(rlx, div(-0.3 * s.adw, dadw))
            s.rlx = rlx

            # ---- update circulation and blade angles
            rms = 0.0
            for i in range(ii):
                dgam = float(s.dgam[i])
                s.gam[i] = s.gam[i] + rlx * dgam
                s.beta[i] = s.beta[i] + rlx * dbet
                s.beta0[i] = s.beta0[i] + rlx * dbet
                rms = rms + div(dgam * dgam, 1.0 + div(1.0, s.adv * s.adv))
                s.dgamold[i] = dgam
            s.dbeta = s.dbeta + rlx * dbet
            s.adv = s.adv + rlx * dadv
            s.adw = s.adw + rlx * dadw
            s.rms = fm.sqrt(div(rms, float(ii)))

            # ---- smooth GAM after heavily under-relaxed steps
            if rlx < 0.2:
                filter_(s.gam, 0.2 * ii, ii)

            # ---- converged: final update of everything for this solution
            if s.rms <= eps:
                if s.free:
                    induced_model(s)
                vcalc(s)
                tpq(s, 1)
                s.conv = True
                return


def aper(s: State, ispec: int, icon: int, linit: bool) -> None:
    """Port of SUBROUTINE APER(ISPEC,ICON,LINIT): optional initialization, then APITER.
    On failure CONV stays False (XROTOR prints 'Iteration limit exceeded')."""
    if linit:
        apinit(s)
    apiter(s, ispec, icon)
