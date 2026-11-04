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

VCALC and HELICO are compiled (_speedups.pyx); the functions here call it.

Vectorization: VRTXCO sums ~10^3 segment velocities per station pair. Segment velocities
are computed with NumPy element-wise (bit-identical to scalar code); the running sums use
np.add.accumulate, which adds strictly left to right like the Fortran loop (np.sum does not).
"""
from __future__ import annotations

import math

import numpy as np

from . import _speedups
from . import fmath as fm
from .fmath import div
from .state import State

PI_EXACT = 4.0 * math.atan(1.0)          # PI = 4.0*ATAN(1.0) in GRADMO/HELICO/VRTXCO
PI4 = 12.56637062                         # DATA PI4 / 12.56637062 / in VORSEGVEL
_VRTXCO_CHUNK = 4096                       # segments per VORSEGVEL call in VRTXCO (speed only)
# d(A x B)(k)/dA(l) = sign * B(m) (and d/dB(l) = -sign * A(m)), for k != l
_CROSS = {(0, 1): (1, 2), (0, 2): (-1, 1), (1, 0): (-1, 2), (1, 2): (1, 0),
          (2, 0): (1, 1), (2, 1): (-1, 0)}


def vcalc(s: State) -> None:
    """Port of SUBROUTINE VCALC: VIND(k,I) = sum over J of VIND_GAM(k,I,J)*GAM(J)."""
    _speedups.vcalc(s)


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
    _speedups.helico(s, _gam_extended(s))


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
        # AXB_A(k,l) = d(AXB(k))/dA(l), AXB_B likewise: only the entries asked for are built
        def axb_a(k, l):
            if k == l:
                return zero
            sign, m = _CROSS[k, l]
            return b[m] if sign > 0 else -b[m]

        def axb_b(k, l):
            if k == l:
                return zero
            sign, m = _CROSS[k, l]
            return -a[m] if sign > 0 else a[m]
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
                                    + axb_a(k, l) * t, 0.0) / PI4
                ub[k, l] = np.where(live, axb[k] * t_bsq * b[l] * 2.0 + axb[k] * t_adb * a[l]
                                    + axb_b(k, l) * t, 0.0) / PI4
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
    # IFIX truncates toward zero; a non-finite or huge THET (ADW NaN, 0 or tiny) gives -2**31,
    # and N1-1 then wraps to +2**31-1 in the compiled reference
    n1 = fm.ifix(2.0 * thet1 / (dth2 + dth1))
    n2 = fm.ifix(2.0 * (thet2 - thet1) / (dth3 + dth2))
    ddth1 = div(dth2 - dth1, float(fm.iwrap(n1 - 1)))
    ddth2 = div(dth3 - dth2, float(fm.iwrap(n2 - 1)))

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

    # legs' segments (leg J, then blade, then L) grouped into chunks of several legs, so a
    # station I needs a few VORSEGVEL calls of a cache-friendly size instead of one per leg;
    # the per-leg sums below stay sequential in the Fortran's order
    nseg = legs[0][0].shape[1]
    per = max(1, _VRTXCO_CHUNK // max(nseg, 1))         # nseg is 0 when ADW is NaN
    chunks = []
    for j0 in range(0, ii + 1, per):
        grp = legs[j0:j0 + per]
        chunks.append((j0, len(grp),
                       np.concatenate([leg[0] for leg in grp], axis=1),
                       np.concatenate([leg[1] for leg in grp], axis=1),
                       np.concatenate([leg[2] for leg in grp]),
                       np.concatenate([leg[3] for leg in grp])))

    with np.errstate(invalid="ignore", over="ignore"):
        for i in range(ii):
            xi = float(s.xi[i])
            r0 = np.array([xi * tanrak, xi, 0.0])[:, None]
            vsums = np.empty((3, ii + 1))
            vadws = np.empty((3, ii + 1))
            for j0, n, r1a, r2a, r1a_adw, r2a_adw in chunks:
                uvw, ua1, ub1 = vorsegvel(r1a - r0, r2a - r0, full=False)
                lead = np.zeros((3, n, 1))
                # VSUM = VSUM + UVW(k) over the leg's segments: 0 + v1 + v2 + ... per leg J.
                # np.add.accumulate is sequential; the leading 0.0 makes an all-(-0.0) sum +0.0
                uvw = uvw.reshape(3, n, nseg)
                vsums[:, j0:j0 + n] = np.add.accumulate(
                    np.concatenate([lead, uvw], axis=2), axis=2)[:, :, -1]
                # VADW = VADW + UVW_A(k,1)*R1_ADW + UVW_B(k,1)*R2_ADW, step by step
                terms = np.stack([ua1 * r1a_adw, ub1 * r2a_adw], axis=-1).reshape(3, n, 2 * nseg)
                vadws[:, j0:j0 + n] = np.add.accumulate(
                    np.concatenate([lead, terms], axis=2), axis=2)[:, :, -1]
            for j in range(ii + 1):
                vsum = vsums[:, j].copy()
                vadw = vadws[:, j].copy()
                vsum[2] = -vsum[2]
                vadw[2] = -vadw[2]
                # ---- open wake: alternate +- influence of each vortex line
                if j < ii:
                    s.vind_gam[:, i, j] = -vsum
                    s.vind_adw[:, i] = s.vind_adw[:, i] - vadw * s.gam[j]
                if j > 0:
                    s.vind_gam[:, i, j - 1] = s.vind_gam[:, i, j - 1] + vsum
                    s.vind_adw[:, i] = s.vind_adw[:, i] + vadw * s.gam[j - 1]
