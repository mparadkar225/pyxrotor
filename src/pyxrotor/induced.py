"""Induced-velocity models: circulation -> induced velocity influence coefficients.

Ported routines:
    VCALC      vcalc()      induced velocities VIND from VIND_GAM and GAM
    GRADMO     gradmo()     graded momentum with Prandtl's tip factor (FAST)
    HELICO     helico()     potential (Goldstein) formulation, spectral/finite difference
    VRTXCO     vrtxco()     discrete helical vortex wake (VRTX)
    VORSEGVEL  vorsegvel()  velocity of a straight vortex segment (vectorized)

Each writes VIND_GAM(3,I,J) (velocity at station I per unit GAM(J)) and VIND_ADW(3,I)
into the State, exactly as XROTOR does from APINIT/APITER/DESI.

Note: GRADMO, HELICO and VRTXCO do not include XROTOR.INC and compute their own
PI = 4*ATAN(1), the exact value, unlike the rest of XROTOR (pyxrotor.constants.PI).

Vectorization: VRTXCO sums ~10^3 segment velocities per station pair. Segment velocities
are computed with NumPy element-wise (bit-identical to scalar code); the running sums use
np.add.accumulate, which adds strictly left to right like the Fortran loop (np.sum does not).
"""
from __future__ import annotations

import math

import numpy as np

from . import fmath as fm
from .fmath import div
from .state import State

PI_EXACT = 4.0 * math.atan(1.0)          # PI = 4.0*ATAN(1.0) in GRADMO/HELICO/VRTXCO
PI4 = 12.56637062                         # DATA PI4 / 12.56637062 / in VORSEGVEL


def _seqsum(values: np.ndarray) -> float:
    """0 + v[0] + v[1] + ... in that order (Fortran's accumulation loop). The explicit
    leading 0.0 matters: 0.0 + (-0.0) is +0.0, so an all-(-0.0) sum is +0.0 as in Fortran."""
    return float(np.add.accumulate(np.concatenate(([0.0], values)))[-1])


def vcalc(s: State) -> None:
    """Port of SUBROUTINE VCALC: VIND(k,I) = sum over J of VIND_GAM(k,I,J)*GAM(J)."""
    ii = s.ii
    vsum = np.zeros((3, ii))
    for j in range(ii):                       # sequential over J, as in the Fortran
        vsum = vsum + s.vind_gam[:, :ii, j] * s.gam[j]
    s.vind[:, :ii] = vsum


def _duct_swirl(s: State, blds: float) -> None:
    """Ducted rotor: simple mean swirl at the blade (shared by GRADMO and VRTXCO)."""
    ii, adw = s.ii, s.adw
    for i in range(ii):
        s.vind_gam[:, i, :ii] = 0.0
        xi = float(s.xi[i])
        s.vind_gam[2, i, i] = div(blds, 4.0 * PI_EXACT * xi)
        s.vind_adw[2, i] = 0.0
        s.vind_adw[1, i] = 0.0
        s.vind_gam[0, i, i] = div(float(s.vind_gam[2, i, i]) * xi, adw)
        s.vind_adw[0, i] = -div(float(s.vind_gam[0, i, i]) * float(s.gam[i]), adw)


def gradmo(s: State) -> None:
    """Port of SUBROUTINE GRADMO: "graded momentum" swirl-circulation coefficients, with
    Prandtl's tip-loss factor F for a free tip."""
    ii, adw = s.ii, s.adw
    blds = float(s.nblds)
    xitip = float(s.xv[ii])
    if s.duct:
        _duct_swirl(s, blds)
        return

    sfac = fm.sqrt(1.0 + div(1.0, adw * adw))
    sf_adw = div(0.5, sfac) * -div(2.0, adw * adw * adw)
    for i in range(ii):
        s.vind_gam[:, i, :ii] = 0.0
        s.vind_adw[:, i] = 0.0
        xi = float(s.xi[i])
        gam = float(s.gam[i])
        tipf = 1.0 - div(xi, xitip)
        arg = fm.fmin(20.0, 0.5 * blds * tipf * sfac)
        ek = fm.exp(-arg)
        ek_adw = -ek * 0.5 * blds * tipf * sf_adw
        fk = fm.sqrt(1.0 - ek * ek)
        fk_adw = div(0.5, fk) * (-2.0 * ek * ek_adw)
        f = fm.atan2(fk, ek) * 2.0 / PI_EXACT
        f_adw = div(ek * fk_adw - fk * ek_adw, ek * ek + fk * fk) * 2.0 / PI_EXACT

        g3 = div(blds, 4.0 * PI_EXACT * f * xi)
        s.vind_gam[2, i, i] = g3
        s.vind_adw[2, i] = g3 * gam * -div(f_adw, f)
        g1 = div(g3 * xi, adw)
        s.vind_gam[0, i, i] = g1
        s.vind_adw[0, i] = div(float(s.vind_adw[2, i]) * xi, adw) - div(g1 * gam, adw)


def _gam_extended(s: State) -> np.ndarray:
    """GAM as HELICO sees it in memory: indices past IX run into the COMMON variables that
    follow GAM in /AERO_X/ (DTII, DPII, ...). HELICO reads GAM(I) up to I = IIMAX-1, beyond
    the II blade stations, so this reproduces exactly which values it reads."""
    return np.concatenate([s.gam, s.dtii, s.dpii, s.dtvi, s.dpvi, s.dtwi, s.dpwi])


def helico(s: State) -> None:
    """Port of SUBROUTINE HELICO: swirl-circulation influence coefficients by the potential
    (Goldstein) formulation. Each of 64 Fourier harmonics of the perturbation potential is
    solved by finite differences (a tridiagonal system), then the series is extrapolated.

    Written with Fortran indices: x[I] is X(I) (declared X(0:IDIM)), an_gam[I, J] is
    AN_GAM(I,J), sys_[k, I] is SYS(k+1, I); station I is s.xi[I-1].
    """
    idim = 150
    ii, nblds, adw = s.ii, s.nblds, s.adw
    gam = _gam_extended(s)                     # GAM(I) is gam[I-1]
    nn = 128

    x = np.zeros(idim + 1)
    x[1:ii + 1] = s.xi[:ii]
    xi0 = float(s.xv[0])
    xitip = float(s.xv[ii])
    x[0] = 2.0 * xi0 - float(s.xi[0])
    x[ii + 1] = 2.0 * xitip - float(s.xi[ii - 1])

    if s.duct:
        iimax = ii + 1
    else:
        # ---- position of outermost point (at "infinity")
        xinf = 4.0 * xitip
        x[ii + 2] = 2.0 * xitip - float(s.xi[ii - 2])
        x[ii + 3] = 2.0 * xitip - float(s.xi[ii - 3])
        x[ii + 4] = 2.0 * xitip - float(s.xi[ii - 4])
        xfac = div(x[ii + 4] - x[ii + 3], x[ii + 3] - x[ii + 2])
        dx = (x[ii + 4] - x[ii + 3]) * xfac
        for i in range(ii + 5, idim):
            x[i] = x[i - 1] + dx
            if x[i] >= xinf:
                break
            dx = dx * xfac
        else:
            i = idim                           # 'HELICO: Local array too small.'
        iimax = i

    s.vind_gam[:, :ii, :ii] = 0.0
    s.vind_adw[:, :ii] = 0.0

    # ==== tridiagonal system
    adwinv = div(1.0, adw * adw)
    qbsq = 0.25 * float(nblds * nblds)
    sys_ = np.zeros((4, idim + 1))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for i in range(1, iimax):
            xm, x0, xp = float(x[i - 1]), float(x[i]), float(x[i + 1])
            sys_[0, i] = div(x0 + xm, x0 - xm)
            sys_[2, i] = div(xp + x0, xp - x0)
            sys_[1, i] = qbsq * (div(1.0, x0) + x0 * adwinv) * (xp - xm)
            sys_[3, i] = qbsq * div(-2.0 * x0 * adwinv, adw) * (xp - xm)
        sys_[0, iimax] = -1.0 if s.duct else 1.0

        ainv = np.zeros(idim + 1)
        csav = np.zeros(idim + 1)
        an_gam = np.zeros((idim + 1, idim + 1))
        an_adw = np.zeros(idim + 1)
        vg3 = np.zeros((ii + 1, ii + 1))       # VIND_GAM(3,I,J) with Fortran I, J
        va3 = np.zeros(ii + 1)                 # VIND_ADW(3,I)
        diag = np.arange(1, ii + 1)

        # ==== loop over harmonics n = 2, 4, ..., NN
        for n in range(2, nn + 1, 2):
            rn = float(n)
            rn2 = rn * rn
            # ---- inner BC dAn/dx = 0, interior equations, outer BC
            ainv[0] = 1.0
            csav[0] = -1.0
            for i in range(1, iimax):
                b = sys_[0, i]
                a = -(sys_[0, i] + sys_[1, i] * rn2 + sys_[2, i])
                c = sys_[2, i]
                ainv[i] = div(1.0, float(a - b * csav[i - 1]))
                csav[i] = c * ainv[i]
            b = sys_[0, iimax]
            ainv[iimax] = div(1.0, float(1.0 - b * csav[iimax - 1]))

            # ==== An and dAn(i)/dGam(j): right-hand sides
            an_gam[:iimax + 1, :ii + 1] = 0.0
            an_gam[1:ii + 1, 0] = sys_[1, 1:ii + 1] * rn2 * gam[:ii]
            an_gam[diag, diag] = sys_[1, 1:ii + 1] * rn2

            # ---- forward elimination, only up to the nonzero elements
            for i in range(1, iimax + 1):
                jl = min(i, ii) + 1
                an_gam[i, :jl] = (an_gam[i, :jl] - sys_[0, i] * an_gam[i - 1, :jl]) * ainv[i]
            for i in range(iimax - 1, -1, -1):
                an_gam[i, :ii + 1] = an_gam[i, :ii + 1] - csav[i] * an_gam[i + 1, :ii + 1]

            # ==== dAn(i)/dAdw
            an_adw[0] = 0.0
            an_adw[1:iimax] = sys_[3, 1:iimax] * rn2 * (gam[:iimax - 1] + an_gam[1:iimax, 0])
            an_adw[iimax] = 0.0
            for i in range(1, iimax + 1):
                an_adw[i] = (an_adw[i] - sys_[0, i] * an_adw[i - 1]) * ainv[i]
            for i in range(iimax - 1, -1, -1):
                an_adw[i] = an_adw[i] - csav[i] * an_adw[i + 1]

            # ---- sum potential harmonics into the swirl-circulation coefficients
            vg3[1:, 1:] = vg3[1:, 1:] + an_gam[1:ii + 1, 1:ii + 1]
            vg3[diag, diag] = vg3[diag, diag] + 1.0
            va3[1:] = va3[1:] + an_adw[1:ii + 1]

        # ---- extrapolate the series assuming (An + Gam) ~ 1/n^2
        if not s.duct:
            fsum = 0.0
            for n in range(nn + 2, 4 * nn + 1, 2):
                t = float(nn) / float(n)
                fsum = fsum + t * t
            vg3[1:, 1:] = vg3[1:, 1:] + an_gam[1:ii + 1, 1:ii + 1] * fsum
            vg3[diag, diag] = vg3[diag, diag] + 1.0 * fsum
            va3[1:] = va3[1:] + an_adw[1:ii + 1] * fsum

        # ---- sawtooth self-influence term and scaling
        for i in range(1, ii + 1):
            bfac = div(float(nblds), 2.0 * PI_EXACT * float(x[i]))
            vg3[i, i] = vg3[i, i] + 0.5
            vg3[i, 1:] = vg3[i, 1:] * bfac
            va3[i] = va3[i] * bfac

        s.vind_gam[2, :ii, :ii] = vg3[1:, 1:]
        s.vind_adw[2, :ii] = va3[1:]

        # ---- other velocity components VX, VY from VZ
        for i in range(ii):
            xi = float(s.xi[i])
            s.vind_gam[0, i, :ii] = s.vind_gam[2, i, :ii] * xi / adw
            s.vind_gam[1, i, :ii] = 0.0
            vsum = 0.0
            for j in range(ii):
                vsum = vsum + float(s.gam[j]) * float(s.vind_gam[2, i, j])
            s.vind_adw[0, i] = div(float(s.vind_adw[2, i]) * xi, adw) - div(vsum * xi, adw * adw)
            s.vind_adw[1, i] = 0.0


def vorsegvel(a: np.ndarray, b: np.ndarray, full: bool = True):
    """Port of SUBROUTINE VORSEGVEL(A,B,UVW,UVW_A,UVW_B), vectorized over segments.

    a, b: (3, M) segment end points relative to the field point. Returns uvw (3, M) and the
    sensitivities uvw_a, uvw_b: (3, 3, M) if full, else only their first column (3, M),
    which is all VRTXCO uses.
    """
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        asq = a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
        bsq = b[0] * b[0] + b[1] * b[1] + b[2] * b[2]
        amag = np.sqrt(asq)
        bmag = np.sqrt(bsq)
        live = (amag * bmag) != 0.0            # NaN counts as nonzero, as in the Fortran

        axb = np.array([a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0]])
        zero = np.zeros_like(asq)
        axb_a = np.array([[zero, b[2], -b[1]], [-b[2], zero, b[0]], [b[1], -b[0], zero]])
        axb_b = np.array([[zero, -a[2], a[1]], [a[2], zero, -a[0]], [-a[1], a[0], zero]])
        adb = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

        den = amag * bmag + adb
        den_asq = 0.5 * bmag / amag
        den_bsq = 0.5 * amag / bmag
        t = (1.0 / amag + 1.0 / bmag) / den
        t_adb = -t / den
        t_asq = -t / den * den_asq - 0.5 / (den * amag * asq)
        t_bsq = -t / den * den_bsq - 0.5 / (den * bmag * bsq)

        uvw = np.where(live, axb * t, 0.0) / PI4
        cols = range(3) if full else range(1)
        ua = np.zeros((3, len(cols)) + asq.shape)
        ub = np.zeros((3, len(cols)) + asq.shape)
        for k in range(3):
            for l in cols:
                ua[k, l] = np.where(live, axb[k] * t_asq * a[l] * 2.0 + axb[k] * t_adb * b[l]
                                    + axb_a[k, l] * t, 0.0) / PI4
                ub[k, l] = np.where(live, axb[k] * t_bsq * b[l] * 2.0 + axb[k] * t_adb * a[l]
                                    + axb_b[k, l] * t, 0.0) / PI4
    if full:
        return uvw, ua, ub
    return uvw, ua[:, 0], ub[:, 0]


def theta_spacing(adw: float) -> np.ndarray:
    """VRTXCO's helix angle stations: fine near the blade, coarser downstream, to 50/ADW."""
    dth1 = PI_EXACT / 60.0
    thet1 = div(2.0, adw)
    dth2 = PI_EXACT / 20.0
    thet2 = div(4.0, adw)
    dth3 = PI_EXACT / 8.0
    thet3 = div(50.0, adw)
    n1 = int(2.0 * thet1 / (dth2 + dth1))               # IFIX truncates toward zero
    n2 = int(2.0 * (thet2 - thet1) / (dth3 + dth2))
    ddth1 = div(dth2 - dth1, float(n1 - 1))
    ddth2 = div(dth3 - dth2, float(n2 - 1))

    thetspc = []
    thet, dth = 0.0, dth1
    for _ in range(5000):                               # NTDIM
        if thet < thet1:
            thetspc.append(thet)
            thet = thet + dth
            dth = dth + ddth1
        elif thet < thet2:
            thetspc.append(thet)
            thet = thet + dth
            dth = dth + ddth2
        elif thet < thet3:
            thetspc.append(thet)
            dth = dth3
            thet = thet + dth
        else:
            break
    # (if 5000 is reached XROTOR prints 'Too many vortex segments for spacing array')
    return np.array(thetspc)


def vrtxco(s: State) -> None:
    """Port of SUBROUTINE VRTXCO: influence coefficients from a discrete helical vortex
    wake (one trailing leg from each vortex position XV(J) of each blade)."""
    ii, nblds, adw = s.ii, s.nblds, s.adw
    blds = float(nblds)
    xitip = float(s.xv[ii])
    tanrak = math.tan(s.rake)

    s.vind_gam[:, :ii, :ii] = 0.0
    s.vind_adw[:, :ii] = 0.0
    thetspc = theta_spacing(adw)
    nthet = thetspc.size

    if s.duct:
        _duct_swirl(s, blds)
        return

    dtbld = 2.0 * PI_EXACT / float(nblds)
    # ---- trailing-leg points for every leg J and blade N (independent of station I)
    legs = []
    for j in range(ii + 1):
        rv = float(s.xv[j])
        xxv = rv * tanrak
        px, py, pz, padw = [], [], [], []
        thetoff = 0.0
        for _ in range(nblds):
            ang = [th + thetoff for th in thetspc]
            px.append(xxv + thetspc * xitip * adw)
            py.append(np.array([rv * math.cos(t) for t in ang]))
            pz.append(np.array([rv * math.sin(t) for t in ang]))
            padw.append(thetspc * xitip)
            thetoff = thetoff + dtbld
        # segments L = 1..NTHET-1 of each blade, in loop order (blade, then L)
        r1 = np.array([np.concatenate([p[:-1] for p in q]) for q in (px, py, pz)])
        r2 = np.array([np.concatenate([p[1:] for p in q]) for q in (px, py, pz)])
        r1_adw = np.concatenate([p[:-1] for p in padw])
        r2_adw = np.concatenate([p[1:] for p in padw])
        legs.append((r1, r2, r1_adw, r2_adw))

    with np.errstate(invalid="ignore", over="ignore"):
        for i in range(ii):
            xi = float(s.xi[i])
            r0 = np.array([xi * tanrak, xi, 0.0])[:, None]
            for j in range(ii + 1):
                r1, r2, r1_adw, r2_adw = legs[j]
                uvw, ua1, ub1 = vorsegvel(r1 - r0, r2 - r0, full=False)
                vsum = np.array([_seqsum(uvw[k]) for k in range(3)])
                # VADW = VADW + UVW_A(k,1)*R1_ADW + UVW_B(k,1)*R2_ADW, step by step
                vadw = np.array([_seqsum(np.ravel(np.column_stack([ua1[k] * r1_adw, ub1[k] * r2_adw])))
                                 for k in range(3)])
                vsum[2] = -vsum[2]
                vadw[2] = -vadw[2]
                # ---- open wake: alternate +- influence of each vortex line
                if j < ii:
                    s.vind_gam[:, i, j] = -vsum
                    s.vind_adw[:, i] = s.vind_adw[:, i] - vadw * s.gam[j]
                if j > 0:
                    s.vind_gam[:, i, j - 1] = s.vind_gam[:, i, j - 1] + vsum
                    s.vind_adw[:, i] = s.vind_adw[:, i] + vadw * s.gam[j - 1]
