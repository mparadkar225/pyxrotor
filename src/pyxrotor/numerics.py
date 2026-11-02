"""Numerical utilities used by the solver (xutils.f, spline.f, xoper.f FILTER).

Ported routines:
    GAUSS   gauss()    dense linear solve with partial pivoting (in place)
    TRISOL  trisol()   tridiagonal solve (in place)
    SPLINE  spline()   cubic-spline derivatives, zero-2nd-derivative ends
    SEVAL   seval()    spline evaluation
    DEVAL   deval()    spline derivative
    SEVLIN  sevlin()   piecewise-linear evaluation and slope
    FILTER  filter_()  implicit smoothing of an array (in place)

Element-wise row operations are vectorized with NumPy; this is bit-identical because each
element gets the same single multiply and subtract as in the Fortran. Reductions (sums along
a row) stay as sequential loops to keep the Fortran's summation order.
"""
from __future__ import annotations

import numpy as np

from .fmath import div


def gauss(z: np.ndarray, r: np.ndarray, nn: int) -> None:
    """Port of SUBROUTINE GAUSS(NSIZ,NN,Z,R,NRHS). Solves Z x = R for the leading NN x NN block.

    z (NSIZ x NSIZ) is destroyed; r (NSIZ x NRHS, or a 1-D NSIZ vector) is replaced by the
    solution(s). Like the Fortran, a singular pivot produces inf/NaN rather than an error.
    """
    if r.ndim == 1:
        r = r[:, None]                                 # a view: updates the caller's array
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for np_ in range(nn - 1):                       # NP = np_+1
            np1 = np_ + 1
            # ---- find max pivot index NX (arithmetic IF: only '>' moves it; NaN -> moves)
            nx = np_
            for n in range(np1, nn):
                d = abs(z[n, np_]) - abs(z[nx, np_])
                if not (d < 0.0 or d == 0.0):
                    nx = n
            pivot = div(1.0, float(z[nx, np_]))
            # ---- switch pivots
            z[nx, np_] = z[np_, np_]
            # ---- switch rows & normalize pivot row
            temp = z[nx, np1:nn] * pivot
            z[nx, np1:nn] = z[np_, np1:nn]
            z[np_, np1:nn] = temp
            temp = r[nx, :] * pivot
            r[nx, :] = r[np_, :]
            r[np_, :] = temp
            # ---- forward eliminate everything
            for k in range(np1, nn):
                ztmp = z[k, np_]
                z[k, np1:nn] = z[k, np1:nn] - ztmp * z[np_, np1:nn]
                r[k, :] = r[k, :] - ztmp * r[np_, :]
        # ---- solve for last row
        r[nn - 1, :] = r[nn - 1, :] / z[nn - 1, nn - 1]
        # ---- back substitute everything (sequential sum over K, vectorized over RHS)
        for np_ in range(nn - 2, -1, -1):
            for k in range(np_ + 1, nn):
                r[np_, :] = r[np_, :] - z[np_, k] * r[k, :]


def trisol(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray, kk: int) -> None:
    """Port of SUBROUTINE TRISOL. A diagonal, B sub-, C super-diagonal; D replaced by the
    solution. A and C are destroyed."""
    for k in range(1, kk):
        km = k - 1
        c[km] = div(float(c[km]), float(a[km]))
        d[km] = div(float(d[km]), float(a[km]))
        a[k] = a[k] - b[k] * c[km]
        d[k] = d[k] - b[k] * d[km]
    d[kk - 1] = div(float(d[kk - 1]), float(a[kk - 1]))
    for k in range(kk - 2, -1, -1):
        d[k] = d[k] - c[k] * d[k + 1]


def spline(x, s, n: int) -> np.ndarray:
    """Port of SUBROUTINE SPLINE(X,XS,S,N): dX/dS at each point, zero second derivative at
    both ends. Returns XS."""
    xs = np.zeros(n)
    if n == 1:
        return xs
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    for i in range(1, n - 1):
        dsm = s[i] - s[i - 1]
        dsp = s[i + 1] - s[i]
        b[i] = dsp
        a[i] = 2.0 * (dsm + dsp)
        c[i] = dsm
        xs[i] = 3.0 * (div((x[i + 1] - x[i]) * dsm, dsp) + div((x[i] - x[i - 1]) * dsp, dsm))
    # ---- zero second derivative end conditions
    a[0] = 2.0
    c[0] = 1.0
    xs[0] = div(3.0 * (x[1] - x[0]), s[1] - s[0])
    b[n - 1] = 1.0
    a[n - 1] = 2.0
    xs[n - 1] = div(3.0 * (x[n - 1] - x[n - 2]), s[n - 1] - s[n - 2])
    trisol(a, b, c, xs, n)
    return xs


def _interval(ss: float, s, n: int) -> int:
    """SEVAL/SEVLIN bisection: the (0-based) upper point I of the interval holding ss."""
    ilow = 0
    i = n - 1
    while i - ilow > 1:
        imid = (i + ilow + 2) // 2 - 1                # (I+ILOW)/2 with 1-based I, ILOW
        if ss < s[imid]:
            i = imid
        else:
            ilow = imid
    return i


def seval(ss: float, x, xs, s, n: int) -> float:
    """Port of FUNCTION SEVAL(SS,X,XS,S,N)."""
    if n == 1:
        return float(x[0])
    ss = float(ss)
    i = _interval(ss, s, n)
    ds = float(s[i] - s[i - 1])
    t = div(ss - float(s[i - 1]), ds)
    cx1 = ds * float(xs[i - 1]) - float(x[i]) + float(x[i - 1])
    cx2 = ds * float(xs[i]) - float(x[i]) + float(x[i - 1])
    return t * float(x[i]) + (1.0 - t) * float(x[i - 1]) + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)


def sevlin(ss: float, x, s, n: int) -> tuple[float, float]:
    """Port of SUBROUTINE SEVLIN(SS,X,S,N,XX,XXS): linear interpolation. Returns (XX, XXS)."""
    ss = float(ss)
    i = _interval(ss, s, n)
    ds = float(s[i] - s[i - 1])
    t = div(ss - float(s[i - 1]), ds)
    xx = t * float(x[i]) + (1.0 - t) * float(x[i - 1])
    xxs = div(float(x[i] - x[i - 1]), ds)
    return xx, xxs


def filter_(q: np.ndarray, smlen: float, n: int) -> None:
    """Port of SUBROUTINE FILTER(Q,SMLEN,N): smooth Q in place over about SMLEN points."""
    if n > 500:
        return                                        # 'FILTER: Array overflow. No action taken'
    con = smlen * smlen
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    a[0] = 1.0
    c[0] = 0.0
    for i in range(1, n - 1):
        b[i] = -con
        a[i] = 2.0 * con + 1.0
        c[i] = -con
    a[n - 1] = 1.0
    b[n - 1] = 0.0
    trisol(a, b, c, q, n)


def deval(ss: float, x, xs, s, n: int) -> float:
    """Port of FUNCTION DEVAL(SS,X,XS,S,N): dX/dS of the spline at SS."""
    if n == 1:
        return float(xs[0])
    ss = float(ss)
    i = _interval(ss, s, n)
    ds = float(s[i] - s[i - 1])
    t = div(ss - float(s[i - 1]), ds)
    cx1 = ds * float(xs[i - 1]) - float(x[i]) + float(x[i - 1])
    cx2 = ds * float(xs[i]) - float(x[i]) + float(x[i - 1])
    d = float(x[i]) - float(x[i - 1]) + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2
    return div(d, ds)
