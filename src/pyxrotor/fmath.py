"""Scalar math with gfortran's IEEE behavior.

Python's float operations raise where compiled Fortran quietly returns inf or NaN
(x/0, overflow in pow/exp, sqrt or log of a negative number, sin of inf). XROTOR's Newton
iterations can diverge through exactly those values, and the port must follow the Fortran
there too, so ported code uses these functions instead of bare operators where a zero
divisor, overflow, or invalid argument is possible.

MIN/MAX with a NaN argument are unspecified in gfortran. At -O1 on x86-64 each call
compiles to MINSD/MAXSD, which return their *source* operand when either operand is NaN,
and which Fortran argument becomes the source varies by call site: a constant or a value
in memory always does, otherwise it depends on register allocation. So each site is
checked in the compiled reference with tools/minmax_sites.py and ported with
    fmin/fmax    NaN -> first argument   (same as Python's built-in min/max)
    fmin2/fmax2  NaN -> second argument
e.g. APINIT's MAX(VHSQ,-0.25) returns -0.25 for a NaN VHSQ (fmax2), while CLCDCM's
MAX(0.0, x) returns 0.0 (fmax). Pin new sites with a non-finite test.
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


def _odd_integer(y: float) -> bool:
    return math.isfinite(y) and float(y).is_integer() and int(y) % 2 == 1


def pow(x: float, y: float) -> float:
    """x**y with a REAL exponent, as libm's pow (C99) returns it: overflow -> +-inf; a zero
    base with a negative exponent (a pole) -> +-inf, negative only for -0 and an odd
    integer exponent; a negative base with a non-integer exponent -> nan. Python's
    math.pow raises for all three. (The pole is reached in CLCDCM when REY/REREF
    underflows to 0 with a negative REXP; found by the compiled-vs-pure comparison.)"""
    try:
        return math.pow(x, y)
    except OverflowError:
        return math.copysign(inf, x) if (x < 0 and _odd_integer(y)) else inf
    except ValueError:
        if x == 0.0:                                  # y < 0: pole
            return math.copysign(inf, x) if _odd_integer(y) else inf
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


def fmax2(a: float, b: float) -> float:
    """Fortran MAX(A,B) compiled as MAXSD with B as the source operand: NaN -> B."""
    return a if a > b else b


def fmin2(a: float, b: float) -> float:
    """Fortran MIN(A,B) compiled as MINSD with B as the source operand: NaN -> B."""
    return a if a < b else b


def sign(a: float, b: float) -> float:
    """Fortran SIGN(A,B): |A| with the sign of B (SIGN(1,-0.0) is -1)."""
    return math.copysign(abs(a), b)


INT32_MIN = -2**31


def ifix(x: float) -> int:
    """IFIX / INT of a REAL into INTEGER*4: truncation toward zero. On x86-64 gfortran
    (cvttsd2si) NaN, +-inf and anything outside the INTEGER*4 range give -2**31."""
    if x != x or not (-2147483649.0 < x < 2147483648.0):
        return INT32_MIN
    return int(x)


def iwrap(n: int) -> int:
    """INTEGER*4 arithmetic result wrapped to 32 bits (two's complement, as compiled)."""
    return (n + 2**31) % 2**32 - 2**31
