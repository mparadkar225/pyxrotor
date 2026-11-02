"""Radial discretization of the blade."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import fmath as fm
from .constants import PI


@dataclass
class RadialGrid:
    dt: float            # DT:  spacing in the cosine-spacing parameter t
    t: np.ndarray        # T:   t at each station, length II
    xi: np.ndarray       # XI:  r/R at each station (vortex-segment midpoints), length II
    dxi: np.ndarray      # DXI: radial width of each station, length II
    xv: np.ndarray       # XV:  r/R of the trailing-vortex positions, length II+1


def setx(ii: int, xi0: float, xitip: float = 1.0, ixspac: int = 2) -> RadialGrid:
    """Port of SUBROUTINE SETX (xrotor.f): radial stations from hub (xi0) to tip (xitip).

    ixspac = 2 gives sine spacing (XROTOR's default), anything else cosine spacing.
    Scalar math is used on purpose: it calls the same libm sin/cos as the Fortran, so the
    result matches the reference bit for bit. (NumPy's vectorized trig can differ by an ulp.)
    """
    dt = 0.5 * PI / float(ii)
    t = np.empty(ii)
    xi = np.empty(ii)
    dxi = np.empty(ii)
    xv = np.empty(ii + 1)

    xm = xi0
    xv[0] = xi0
    for i in range(1, ii + 1):
        ti = dt * (float(i) - 0.5)
        tp = dt * float(i)
        if ixspac == 2:
            # XROTOR also computes XI(I) from T(I) here, then overwrites it below.
            s, c = math.sin(tp), math.cos(tp)
            xp = math.sqrt(xitip * (s * s) + (xi0 * c) * (xi0 * c))
        else:
            xp = 0.5 * (1.0 - math.cos(2.0 * tp)) * (xitip - xi0) + xi0
        t[i - 1] = ti
        xi[i - 1] = (xp + xm) * 0.5
        dxi[i - 1] = xp - xm
        xm = xp
        xv[i] = xp
    xv[ii] = xitip
    return RadialGrid(dt, t, xi, dxi, xv)


def tinvrt(xx: float, xi0: float, xitip: float = 1.0, ixspac: int = 2) -> float:
    """Port of FUNCTION TINVRT (xrotor.f): the spacing parameter t at which SETX's
    stretching function gives r/R = xx (Newton, up to 40 iterations, |residual| <= 1e-5)."""
    xx = float(xx)
    tt = fm.atan2(xx, fm.sqrt(xitip * xitip - xx * xx))
    for _ in range(40):
        st = fm.sin(tt)
        ct = fm.cos(tt)
        if ixspac == 2:
            # ---- usual sine stretching, adjusted for nonzero root radius
            a, b = xitip * st, xi0 * ct
            rez = a * a + b * b - xx * xx
            drdt = 2.0 * (xitip * xitip - xi0 * xi0) * st * ct
        else:
            # ---- cosine stretching for more root resolution (also in SETX)
            rez = 0.5 * (1.0 - fm.cos(2.0 * tt)) * (xitip - xi0) + xi0 - xx
            drdt = (xitip - xi0) * fm.sin(2.0 * tt)
        tt = tt - fm.div(rez, drdt)
        if abs(rez) <= 1.0e-5:
            break
    # (else XROTOR prints 'TINVRT: Convergence failed')
    return tt
