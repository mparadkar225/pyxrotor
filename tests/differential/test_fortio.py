"""Fortran formatted output and list-directed input against the compiled reference (gfortran)."""
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor.fortio import fmt_e, fmt_f, fmt_g, fmt_i, getflt, read_list

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)


def ref_fmt(ref, kind, x, w, d=0):
    ref.call("REFFMT", kind, float(x), w, d)
    return ref["FNAME"][()].decode()[:w]          # raw: trailing blanks matter


# values near the decade and rounding boundaries G and F have to get right
boundaries = st.sampled_from([0.0, -0.0, 0.1, 0.0999995, 0.09999949999, 99999.4, 99999.5, 1e5,
                              1.03125, 1.03135, 2.5e-5, 1e100, -1e-100, 1e-310, 5e-324,
                              float("nan"), float("inf"), -float("inf"), 999.0, 7.55, 0.5, 9.99995])
reals = st.one_of(boundaries,
                  st.floats(allow_nan=False, allow_infinity=False),
                  st.floats(-1e7, 1e7, allow_nan=False),
                  st.builds(lambda m, e: m * 10.0 ** e, st.integers(-10**6, 10**6), st.integers(-12, 8)))


@SETTINGS
@given(x=reals, w=st.integers(8, 20), d=st.integers(1, 7))
def test_fmt_g(ref_module, x, w, d):
    if d > w - 6:
        d = max(1, w - 6)
    assert fmt_g(x, w, d) == ref_fmt(ref_module, 1, x, w, d)


@SETTINGS
@given(x=reals, w=st.integers(8, 20), d=st.integers(1, 7))
def test_fmt_e(ref_module, x, w, d):
    assert fmt_e(x, w, d) == ref_fmt(ref_module, 2, x, w, d)


@SETTINGS
@given(x=reals, w=st.integers(3, 20), d=st.integers(0, 7))
def test_fmt_f(ref_module, x, w, d):
    assert fmt_f(x, w, d) == ref_fmt(ref_module, 3, x, w, d)


@SETTINGS
@given(n=st.integers(-10**8, 10**8), w=st.integers(1, 12))
def test_fmt_i(ref_module, n, w):
    assert fmt_i(n, w) == ref_fmt(ref_module, 4, float(n), w)


def _record(draw_vals, seps):
    out = []
    for v, sep in zip(draw_vals, seps):
        out.append(v + sep)
    return "".join(out)


real_tokens = st.one_of(
    reals.filter(math.isfinite).map(lambda x: fmt_g(x, 12, 5).strip()),
    reals.filter(math.isfinite).map(repr),
    st.sampled_from(["1", "-2", "+3.", ".5", "1.0D0", "2.5d-3", "1.0+07", "1E5", "0.30000E+07",
                     "3*1.5", "2*", "-0.0", "1.e-2"]))


@SETTINGS
@given(toks=st.lists(real_tokens, min_size=1, max_size=8),
       seps=st.lists(st.sampled_from([" ", "  ", ",", " , ", "\t"]), min_size=8, max_size=8),
       tail=st.sampled_from(["", " ", " /", "/ junk", ",,"]), n=st.integers(1, 10))
def test_read_list_reals(ref_module, toks, seps, tail, n):
    """READ(FNAME,*) (A(I),I=1,N): same values, same error/end conditions."""
    rec = _record(toks, seps) + tail
    if "\t" in rec:
        rec = rec.replace("\t", " ")
    if len(rec) > 80:
        return
    a = np.full(16, -12345.0)
    ref_module["FNAME"] = rec
    _, out, ierr = ref_module.call("REFRDL", n, a, 0)
    try:
        got = read_list(rec, "r" * n, current=[-12345.0] * n)
        err = 0
    except (ValueError, EOFError):
        err = 1
    assert err == int(ierr)
    if not err:
        np.testing.assert_array_equal(np.array(got), out[:n])


@SETTINGS
@given(toks=st.lists(real_tokens, max_size=14),
       seps=st.lists(st.sampled_from([" ", "   ", ",", " ,"]), min_size=14, max_size=14),
       comment=st.sampled_from(["", " ! note 1 2", "!"]), nmax=st.sampled_from([0, 0, 3, 13]))
def test_getflt(ref_module, toks, seps, comment, nmax):
    rec = (_record(toks, seps) + comment)[:80]
    try:
        vals, err = getflt(rec.ljust(80), nmax)
    except EOFError:
        return                     # 7.55 aborts the process here: test_getflt_aborts
    a = np.zeros(40)
    ref_module["FNAME"] = rec
    out, n, ierr = ref_module.call("REFGETFLT", a, nmax, 0)
    assert err == bool(ierr)
    if not err:
        assert len(vals) == int(n)
        np.testing.assert_array_equal(np.array(vals), out[: int(n)])


@pytest.mark.parametrize("rec", ["1,2 ,3", "0.1 , 0.5 0.1", "4.0 ,"])
def test_getflt_aborts(rec):
    """GETFLT over-counts a comma after a blank; the reference dies with a Fortran runtime
    error (run in a subprocess), and the port raises EOFError."""
    import subprocess
    import sys
    from pathlib import Path
    with pytest.raises(EOFError):
        getflt(rec.ljust(80))
    code = ("import sys, numpy as np; sys.path[:0] = ['.']\n"
            "from reference.ref755 import Ref755\n"
            f"r = Ref755(); r['FNAME'] = {rec!r}; r.call('REFGETFLT', np.zeros(40), 0, 0)\n"
            "print('survived')")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       cwd=Path(__file__).resolve().parents[2])
    assert p.returncode != 0 and "End of file" in p.stderr and "survived" not in p.stdout


@SETTINGS
@given(vals=st.lists(reals, min_size=1, max_size=3))
def test_write_list(ref_module, vals):
    """List-directed WRITE of REAL*8 values (VPUT's file format), up to the trailing blanks."""
    from pyxrotor.fortio import write_list
    ref_module.call("REFLDW", len(vals), np.array(vals, dtype=float))
    got = (ref_module["SAVFIL"][()] + ref_module["FNAME"][()]).decode().rstrip()
    assert write_list(*vals).rstrip() == got
