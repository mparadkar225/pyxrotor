"""The compiled hot loops (pyxrotor._speedups) against their pure-Python counterparts: the
same results bit for bit (NaN payloads aside), including on non-finite inputs. The rest of
the differential suite checks whichever implementation is active against the Fortran; CI
also runs it with PYXROTOR_PURE=1."""
import copy
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from pyxrotor import _accel
from pyxrotor.aero import AeroSection, clcdcm
from pyxrotor.induced import helico
from pyxrotor.numerics import gauss
from strategies import finite
from test_induced import wake_inputs

pytestmark = [pytest.mark.differential,
              pytest.mark.skipif(not _accel.available, reason="compiled speedups not built")]
SETTINGS = settings(deadline=None)


def same(a, b) -> bool:
    """Equal bit patterns, or both NaN."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return False
    nan = np.isnan(a) & np.isnan(b)
    return bool(np.all(nan | (a.view(np.int64) == b.view(np.int64))))


def both(fn, *args):
    """fn(*args) compiled and pure (on deep copies of the arguments)."""
    a1 = copy.deepcopy(args)
    a2 = copy.deepcopy(args)
    _accel.use_compiled(True)
    try:
        r1 = fn(*a1)
    finally:
        _accel.use_compiled(True)
    with _accel.pure():
        r2 = fn(*a2)
    return (r1, a1), (r2, a2)


values = st.one_of(finite(-50.0, 50.0), st.sampled_from([0.0, -0.0, math.inf, -math.inf, math.nan]))


@st.composite
def sections(draw):
    pos = lambda lo, hi: draw(st.one_of(finite(lo, hi), st.sampled_from([0.0, math.nan])))
    return AeroSection(a0=draw(finite(-0.2, 0.1)), clmax=pos(0.5, 2.0), clmin=draw(finite(-1.5, 0.0)),
                       dclda=pos(2.0, 7.0), dclda_stall=pos(0.0, 0.5), dcl_stall=pos(0.0, 0.5),
                       cdmin=pos(0.0, 0.05), cldmin=draw(finite(-0.5, 1.0)), dcdcl2=pos(0.0, 0.05),
                       cmcon=draw(finite(-0.2, 0.1)), reref=draw(finite(1e4, 1e7)),
                       rexp=draw(finite(-1.0, 0.0)), mcrit=pos(0.3, 1.0))


@SETTINGS
@given(alf=values, w=st.one_of(finite(0.0, 15.0), st.sampled_from([0.0, math.inf, math.nan])),
       rey=st.one_of(finite(-1.0, 1e8), st.just(math.nan)), sec=sections())
def test_clcdcm(alf, w, rey, sec):
    (r1, _), (r2, _) = both(clcdcm, alf, w, rey, 120.0, 340.0, sec)
    assert type(r1) is type(r2)
    for k, (a, b) in enumerate(zip(r1, r2)):
        assert same(a, b), (r1._fields[k], a, b)


@SETTINGS
@given(n=st.integers(1, 40), nrhs=st.integers(1, 3), data=st.data())
def test_gauss(n, nrhs, data):
    el = st.one_of(finite(-5.0, 5.0), st.sampled_from([0.0, -0.0]))
    z = data.draw(arrays(np.float64, (n + 2, n + 2), elements=el))
    r = data.draw(arrays(np.float64, (n + 2, nrhs), elements=el))
    if data.draw(st.booleans()):
        z[data.draw(st.integers(0, n - 1)), data.draw(st.integers(0, n - 1))] = math.nan
    (_, (z1, r1, _)), (_, (z2, r2, _)) = both(gauss, z, r, n)
    assert same(z1, z2) and same(r1, r2)


def test_gauss_vector_rhs_and_non_contiguous_fallback():
    rng = np.random.default_rng(3)
    z = rng.standard_normal((12, 12))
    r = rng.standard_normal(12)
    (_, (z1, r1, _)), (_, (z2, r2, _)) = both(gauss, z, r, 10)
    assert same(z1, z2) and same(r1, r2)
    zf = np.asfortranarray(z)                  # not C-contiguous: the Python path runs
    (_, (z1, r1, _)), (_, (z2, r2, _)) = both(gauss, zf, r, 10)
    assert same(z1, z2) and same(r1, r2)


@SETTINGS
@given(s=wake_inputs(), stale=finite(-1.0, 1.0), duct=st.booleans(),
       adw=st.one_of(st.none(), st.sampled_from([0.0, math.nan, math.inf, 1e-300])))
def test_helico(s, stale, duct, adw):
    s.gam[s.ii:] = stale
    s.duct = duct
    if adw is not None:
        s.adw = adw
    (_, (t1,)), (_, (t2,)) = both(helico, s)
    assert same(t1.vind_gam, t2.vind_gam) and same(t1.vind_adw, t2.vind_adw)
