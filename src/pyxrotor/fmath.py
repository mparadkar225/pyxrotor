"""Scalar math with gfortran's IEEE behavior.

Python's float operations raise where compiled Fortran quietly returns inf or NaN
(x/0, overflow in pow/exp, sqrt or log of a negative number, sin of inf). XROTOR's Newton
iterations can diverge through exactly those values, and the port must follow the Fortran
there too, so ported code uses these functions instead of bare operators where a zero
divisor, overflow, or invalid argument is possible.

MIN/MAX with a NaN argument are unspecified in gfortran (GCC 9 release notes) and depend on
how each call site compiles. In the reference library, every site tested so far returns the
*first* argument when either is NaN: MAX(0,NaN) = 0, MIN(NaN,x) = NaN. That is exactly
Python's built-in min/max (keep the first value unless a later one compares smaller/larger),
which fmin/fmax spell out. A standalone probe program behaved differently, so only tests
against the compiled routine itself (test_clcdcm_non_finite) count as evidence.
"""
from __future__ import annotations

import math

inf = math.inf
nan = math.nan


def div(a: float, b: float) -> float:
    """a / b: x/0 -> +-inf, 0/0 and nan/0 -> nan."""
    try:
        return a / b
    except ZeroDivisionError:
        if a == 0.0 or a != a:
            return nan
        return math.copysign(inf, a) * math.copysign(1.0, b)


def pow(x: float, y: float) -> float:
    """x**y with a REAL exponent (libm pow): overflow -> inf, invalid -> nan."""
    try:
        return math.pow(x, y)
    except OverflowError:
        return math.copysign(inf, x) if (x < 0 and float(y).is_integer() and int(y) % 2) else inf
    except ValueError:
        return nan


def sqrt(x: float) -> float:
    return math.sqrt(x) if x >= 0.0 else nan          # also nan for x = nan


def exp(x: float) -> float:
    try:
        return math.exp(x)
    except OverflowError:
        return inf


def log(x: float) -> float:
    if x > 0.0:
        return math.log(x)
    if x == 0.0:
        return -inf
    return nan                                       # negative or nan


def sin(x: float) -> float:
    return math.sin(x) if math.isfinite(x) else nan


def cos(x: float) -> float:
    return math.cos(x) if math.isfinite(x) else nan


atan2 = math.atan2
atan = math.atan


def fmax(a: float, b: float) -> float:
    """Fortran MAX(A,B) as compiled in the reference: A unless B > A (NaN -> A)."""
    return b if b > a else a


def fmin(a: float, b: float) -> float:
    """Fortran MIN(A,B) as compiled in the reference: A unless B < A (NaN -> A)."""
    return b if b < a else a


def sign(a: float, b: float) -> float:
    """Fortran SIGN(A,B): |A| with the sign of B (SIGN(1,-0.0) is -1)."""
    return math.copysign(abs(a), b)
