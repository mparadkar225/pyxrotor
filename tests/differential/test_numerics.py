"""GAUSS, TRISOL, SPLINE, SEVAL, SEVLIN, FILTER, TINVRT against the reference."""
import numpy as np
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from harness import assert_close
from pyxrotor.geometry import tinvrt
from pyxrotor.numerics import deval, filter_, gauss, seval, sevlin, spline, trisol
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)

values = finite(-1e3, 1e3)


def vec(n, elements=values):
    return arrays(np.float64, n, elements=elements)


@st.composite
def increasing(draw, n):
    """n strictly increasing abscissae."""
    steps = draw(arrays(np.float64, n, elements=finite(0.01, 2.0)))
    return np.cumsum(steps) + draw(finite(-5.0, 5.0))


@SETTINGS
@given(nn=st.integers(1, 40), extra=st.integers(0, 3), nrhs=st.integers(1, 3), data=st.data())
def test_gauss(ref_module, nn, extra, nrhs, data):
    nsiz = nn + extra                                    # Fortran leading dimension >= NN
    z = np.asfortranarray(data.draw(arrays(np.float64, (nsiz, nsiz), elements=finite(-10, 10))))
    if data.draw(st.booleans()):
        z += np.eye(nsiz) * 30.0                          # well conditioned
    r = np.asfortranarray(data.draw(arrays(np.float64, (nsiz, nrhs), elements=values)))
    _, _, zf, rf, _ = ref_module.call("GAUSS", nsiz, nn, z.copy(order="F"), r.copy(order="F"), nrhs)
    zp, rp = z.copy(), r.copy()
    gauss(zp, rp, nn)
    assert_close(zp, zf, rtol=0, name="Z")
    assert_close(rp, rf, rtol=0, name="R")


def test_gauss_singular_gives_nonfinite_like_fortran(ref):
    z = np.zeros((3, 3), order="F")
    r = np.ones((3, 1), order="F")
    _, _, zf, rf, _ = ref.call("GAUSS", 3, 3, z.copy(order="F"), r.copy(order="F"), 1)
    zp, rp = z.copy(), r.copy()
    gauss(zp, rp, 3)
    assert_close(rp, rf, rtol=0, name="R")


@SETTINGS
@given(kk=st.integers(1, 60), data=st.data())
def test_trisol(ref_module, kk, data):
    a = data.draw(vec(kk, finite(2.0, 10.0)))
    b, c, d = (data.draw(vec(kk, finite(-1.0, 1.0))) for _ in range(3))
    af, bf, cf, df, _ = ref_module.call("TRISOL", a.copy(), b.copy(), c.copy(), d.copy(), kk)
    ap, cp, dp = a.copy(), c.copy(), d.copy()
    trisol(ap, b.copy(), cp, dp, kk)
    for name, got, want in (("A", ap, af), ("C", cp, cf), ("D", dp, df)):
        assert_close(got, want, rtol=0, name=name)


@SETTINGS
@given(n=st.integers(1, 40), data=st.data())
def test_spline_seval_sevlin(ref_module, n, data):
    s = data.draw(increasing(n))
    x = data.draw(vec(n))
    _, xsf, _, _ = ref_module.call("SPLINE", x.copy(), np.zeros(n), s.copy(), n)
    xs = spline(x, s, n)
    assert_close(xs, xsf, rtol=0, name="XS")
    for ss in data.draw(st.lists(finite(s[0] - 1.0, s[-1] + 1.0), min_size=1, max_size=5)):
        assert seval(ss, x, xs, s, n) == ref_module.func("SEVAL", ss, x.copy(), xsf.copy(), s.copy(), n)
        assert_close(deval(ss, x, xs, s, n), ref_module.func("DEVAL", ss, x.copy(), xsf.copy(), s.copy(), n),
                     rtol=0, name="DEVAL")
        if n >= 2:
            out = ref_module.call("SEVLIN", ss, x.copy(), s.copy(), n, 0.0, 0.0)
            assert sevlin(ss, x, s, n) == (out[4].item(), out[5].item())


@SETTINGS
@given(n=st.integers(2, 120), smlen=finite(0.0, 5.0), data=st.data())
def test_filter(ref_module, n, smlen, data):
    q = data.draw(vec(n))
    qf = ref_module.call("FILTER", q.copy(), smlen, n)[0]
    qp = q.copy()
    filter_(qp, smlen, n)
    assert_close(qp, qf, rtol=0, name="Q")


@SETTINGS
@given(xi0=finite(0.01, 0.5), ixspac=st.sampled_from([1, 2]), frac=finite(0.0, 1.0))
def test_tinvrt(ref_module, xi0, ixspac, frac):
    xx = xi0 + frac * (1.0 - xi0)
    ref_module["XI0"], ref_module["XITIP"], ref_module["IXSPAC"] = xi0, 1.0, ixspac
    assert_close(tinvrt(xx, xi0, 1.0, ixspac), ref_module.func("TINVRT", xx), rtol=0, name="T")
