"""Fortran formatted I/O as gfortran does it, for files that must match XROTOR byte for byte.

Output:  fmt_f (Fw.d), fmt_e (Ew.d), fmt_g (Gw.d), fmt_i (Iw), fmt_l (Lw)
Input:   read_list (list-directed READ of reals, integers, logicals), getflt (userio.f GETFLT)

Checked against the compiled reference (tests/differential/test_fortio.py). gfortran rounds
the exact binary value to nearest, ties to even (1.03125 in F6.4 is 1.0312), decides G's
F-or-E choice after rounding (99999.5 in G12.5 is 0.10000E+06), writes three-digit
exponents without the letter (0.10000+101), and shortens Infinity to Inf (or -Inf) when the
field is narrow.
"""
from __future__ import annotations

import math
import re
from decimal import ROUND_HALF_EVEN, Context, Decimal, localcontext

# every double is exact in at most 767 significant digits
_CTX = Context(prec=1200, Emax=10**6, Emin=-10**6)


def _special(x: float, w: int) -> str | None:
    if x != x:
        return "NaN".rjust(w) if w >= 3 else "*" * w
    if math.isinf(x):
        neg = x < 0
        for text in (("-Infinity", "-Inf") if neg else ("Infinity", "Inf")):
            if len(text) <= w:
                return text.rjust(w)
        return "*" * w
    return None


def fmt_f(x: float, w: int, d: int) -> str:
    """Fw.d"""
    with localcontext(_CTX):
        return _fmt_f(float(x), w, d)


def _fmt_f(x: float, w: int, d: int) -> str:
    sp = _special(x, w)
    if sp is not None:
        return sp
    neg = math.copysign(1.0, x) < 0
    q = Decimal(abs(x)).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_EVEN)
    digits = f"{q:f}"
    if d == 0:
        digits += "."                             # F always writes the decimal point
    if digits.startswith("0") and digits != "0.":
        digits = digits[1:]                      # ".1234", leading zero optional
    text = ("-" if neg else "") + digits
    if len(text) < w and digits.startswith("."):
        text = ("-" if neg else "") + "0" + digits
    if len(text) > w:
        return "*" * w
    return text.rjust(w)


def _sig(x: float, d: int) -> tuple[int, int]:
    """abs(x) rounded to d significant digits: (integer mantissa with d digits, exponent e)
    with value = mantissa * 10**(e-d), i.e. 0.MMMMM x 10**e. x must be finite and nonzero."""
    with localcontext(_CTX):
        return _sig_exact(x, d)


def _sig_exact(x: float, d: int) -> tuple[int, int]:
    dec = Decimal(abs(x))
    e = dec.adjusted() + 1                       # 10**(e-1) <= |x| < 10**e
    m = (dec.scaleb(d - e)).quantize(Decimal(1), rounding=ROUND_HALF_EVEN)
    if m >= 10 ** d:                             # rounded up to the next decade
        m = m.scaleb(-1).quantize(Decimal(1), rounding=ROUND_HALF_EVEN)
        e += 1
    return int(m), e


def fmt_e(x: float, w: int, d: int) -> str:
    """Ew.d (exponent of two digits as E+nn, three as +nnn)"""
    x = float(x)
    sp = _special(x, w)
    if sp is not None:
        return sp
    neg = math.copysign(1.0, x) < 0
    if x == 0.0:
        m, e = 0, 0
    else:
        m, e = _sig(x, d)
    mant = "." + str(m).rjust(d, "0")
    if abs(e) <= 99:
        exp = "E" + ("-" if e < 0 else "+") + f"{abs(e):02d}"
    elif abs(e) <= 999:
        exp = ("-" if e < 0 else "+") + f"{abs(e):03d}"
    else:
        return "*" * w
    text = ("-" if neg else "") + mant + exp
    if len(text) < w:
        text = ("-" if neg else "") + "0" + mant + exp
    if len(text) > w:
        return "*" * w
    return text.rjust(w)


def fmt_g(x: float, w: int, d: int) -> str:
    """Gw.d for reals: F(w-4).(d-k) plus 4 blanks when 0.1 <= |x| < 10**d after rounding
    to d significant digits (k the number of integer digits), otherwise Ew.d."""
    x = float(x)
    sp = _special(x, w)
    if sp is not None:
        return sp
    if x == 0.0:
        return fmt_f(x, w - 4, d - 1) + "    "
    m, e = _sig(x, d)
    if 0 <= e <= d:
        return fmt_f(x, w - 4, d - e) + "    "
    return fmt_e(x, w, d)


def fmt_i(n: int, w: int) -> str:
    text = str(int(n))
    return "*" * w if len(text) > w else text.rjust(w)


def fmt_l(v: bool, w: int) -> str:
    return ("T" if v else "F").rjust(w)


# ------------------------------------------------------------------------------- input
_REAL = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([EeDdQq]?[+-]?\d+)?$")


def _real(tok: str) -> float:
    t = tok.strip()
    if not _REAL.match(t):
        low = t.lower().lstrip("+-")
        if low in ("nan", "inf", "infinity"):
            return float(t)
        raise ValueError(f"bad real {tok!r}")
    mant, exp = re.match(r"^([+-]?(?:\d+\.?\d*|\.\d+))(.*)$", t).groups()
    if exp:
        exp = exp.lstrip("EeDdQq")
        return float(f"{mant}e{exp}")
    return float(mant)


def _logical(tok: str) -> bool:
    t = tok.strip().upper()
    if t.startswith("."):
        t = t[1:]
    if t[:1] == "T":
        return True
    if t[:1] == "F":
        return False
    raise ValueError(f"bad logical {tok!r}")


def _integer(tok: str) -> int:
    t = tok.strip()
    if not re.match(r"^[+-]?\d+$", t):
        raise ValueError(f"bad integer {tok!r}")
    return int(t)


def _tokens(record: str):
    """List-directed values: separators are blanks and one comma, '/' ends the read; r*c
    repeats c, r* (or an empty field between commas) is a null value (None)."""
    out = []
    s = record
    i, n = 0, len(s)
    while i < n:
        while i < n and s[i] == " ":
            i += 1
        if i >= n:
            break
        if s[i] == "/":
            break
        if s[i] == ",":                          # null value
            out.append(None)
            i += 1
            continue
        j = i
        while j < n and s[j] not in " ,/":
            j += 1
        tok = s[i:j]
        k = j
        while k < n and s[k] == " ":
            k += 1
        if k < n and s[k] == ",":
            k += 1
        i = k
        m = re.match(r"^(\d+)\*(.*)$", tok)
        if m:
            r, c = int(m.group(1)), m.group(2)
            out.extend([c if c else None] * r)
        else:
            out.append(tok)
    return out


def read_list(record: str, kinds: str, current=None) -> list:
    """READ(record,*) of len(kinds) items; kinds is a string of 'r' (real), 'i' (integer),
    'l' (logical). A null value or a '/' leaves the item at current[k] (None if not given).
    Raises ValueError on a bad value (Fortran's ERR condition) and EOFError when the record
    has too few values (END)."""
    toks = _tokens(record)
    slash = "/" in record and len(toks) < len(kinds)
    vals = []
    for k, kind in enumerate(kinds):
        if k >= len(toks):
            if slash:
                vals.append(None if current is None else current[k])
                continue
            raise EOFError("end of record")
        t = toks[k]
        if t is None:
            vals.append(None if current is None else current[k])
            continue
        vals.append(_real(t) if kind == "r" else _integer(t) if kind == "i" else _logical(t))
    return vals


def getflt(text: str, nmax: int = 0, limit: int = 128) -> tuple[list[float], bool]:
    """Port of SUBROUTINE GETFLT (userio.f): (values, error). Counts numbers in the first
    128 characters the way 7.55 does, then list-directed READs that many (at most nmax).

    7.55: the count treats a comma after a blank ("1,2 ,3") as one more number, and the
    READ has ERR= but no END=, so XROTOR stops with 'Fortran runtime error: End of file'.
    That case raises EOFError here."""
    ilen = min(len(text), limit)                   # 7.55: only 128 characters parsed
    rec = text[:ilen].ljust(ilen) + " ,"
    k = rec.find("!")
    if k >= 0:
        rec = rec[:k].ljust(ilen) + rec[ilen:]
    ilenp = ilen + 2
    n = 0
    k = 1                                          # 1-based, as in the Fortran
    for _ in range(ilen):
        sub = rec[k - 1:ilenp]
        kspace = sub.find(" ") + 1 + k - 1 if " " in sub else k - 1
        kcomma = sub.find(",") + 1 + k - 1 if "," in sub else k - 1
        if k == kspace:
            k += 1
        elif k == kcomma:
            n += 1
            k += 1
        else:
            n += 1
            k = min(kspace, kcomma) + 1
        if k >= ilen:
            break
    if nmax > 0:
        n = min(n, nmax)
    try:
        vals = read_list(rec[:ilen], "r" * n)
    except EOFError:
        raise                                    # 7.55 aborts (no END= on the READ)
    except ValueError:
        return [], True
    return [0.0 if v is None else v for v in vals], False
