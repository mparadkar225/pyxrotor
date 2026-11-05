"""CR-176883's airfoil data banks (pyxrotor.databank) against the compiled transcription
(reference/cr176883, REAL*8): every output bit for bit, NaN where the Fortran gives NaN.

Inputs range past the banks' own limits (Mach to 1.3, t/c to 0.45, negative design CL)
so that the extrapolating and supersonic branches are covered too.
"""
import math

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor import databank
from reference import crdb

pytestmark = [pytest.mark.differential,
              pytest.mark.skipif(not crdb.available(), reason="reference/build.sh not run")]
SETTINGS = settings(deadline=None)


@pytest.fixture(scope="module")
def lib():
    return crdb.CRDB()


def same(got, want, name):
    got, want = float(got), float(want)
    if math.isnan(want):
        assert math.isnan(got), f"{name}: got {got!r}, want nan"
    else:
        assert got == want, f"{name}: got {got!r}, want {want!r} (diff {got - want:.3g})"


def fl(lo, hi):
    return st.floats(lo, hi, allow_nan=False)


altitudes = fl(-1000.0, 60000.0)


@SETTINGS
@given(alt=altitudes)
def test_atmosphere(lib, alt):
    out = lib.call("ATMCON", alt, 0.0, 0.0, 0.0, 0.0, 0.0)
    got = databank.atmosphere(alt)
    for name, want in zip(("temp", "pres", "den", "vis", "wa"), out[1:]):
        same(getattr(got, name), want, name)
    # the 4-digit bank's ATMCON (other argument order) is the same arithmetic
    out4 = lib.call("ATMCN4", alt, 0.0, 0.0, 0.0)
    for name, want in zip(("den", "vis", "wa"), out4[1:]):
        same(getattr(got, name), want, "ATMCN4 " + name)


@SETTINGS
@given(xm=st.one_of(fl(0.0, 1.3), st.sampled_from([0.95, 1.0, 0.0])),
       alpha=fl(-10.0, 15.0), icl=st.sampled_from([0, 1]), chord=fl(0.05, 5.0),
       tc=st.one_of(fl(0.0, 0.45), st.sampled_from([0.03, 0.031, 0.4])),
       cld=fl(-0.3, 1.0), xr=fl(0.0, 1.0), alt=altitudes)
def test_clcd16(lib, xm, alpha, icl, chord, tc, cld, xr, alt):
    atm = databank.atmosphere(alt)
    out = lib.call("CLCD16", xm, alpha, 0.0, 0.0, icl, 1, chord, tc, 0.0, cld, xr,
                   atm.wa, atm.vis, atm.den, 0.0, 0.0, 0.0, 0.0)
    got = databank.clcd16(xm, alpha, icl, chord, tc, cld, xr, atm.wa, atm.vis, atm.den)
    want = dict(cl=out[2], cd=out[3], re=out[8], cd_friction=out[14], cd_lift=out[15],
                cd_camber=out[16], cd_zero=out[17])
    for name, w in want.items():
        same(getattr(got, name), w, name)


@SETTINGS
@given(camp=st.sampled_from([0, 1]), alpha=st.one_of(fl(-12.0, 12.0), st.just(0.0)),
       mach=fl(0.0, 0.95), tc=fl(0.01, 0.3), cz=fl(0.0, 1.0), chord=fl(0.1, 5.0),
       alt=altitudes, re=st.one_of(st.just(0.0), fl(5e4, 2e7)), fix=st.booleans())
def test_naca4(lib, camp, alpha, mach, tc, cz, chord, alt, re, fix):
    cld, z = (cz, 0.0) if camp == 0 else (0.0, cz * 0.1)
    out = lib.call("NACA4", float(camp), tc, cld, z, alpha, mach, re, alt, chord, fix,
                   *[0.0] * 7)
    kw = dict(cld=cld) if camp == 0 else dict(z=z)
    got = databank.naca4(alpha, mach, tc, chord_ft=chord, altitude_ft=alt,
                         re=re if re > 0 else None, fix_x8=fix, **kw)
    want = dict(cld=out[2], z=out[3], re=out[6], cl=out[10], cd=out[11], clmax=out[12],
                clopt=out[13], x=out[14], cdmin=out[15], delcd=out[16])
    for name, w in want.items():
        same(getattr(got, name), w, name)


def test_x8_fix_matters(lib):
    """The listing's X8(13) = 1. changes CDL between X = 0 and 1 on cambered (Z > 0.04)
    sections; the fix interpolates the table as drawn."""
    a = databank.cdl(0.1, 0.07, 0.12)
    b = databank.cdl(0.1, 0.07, 0.12, fix_x8=True)
    assert a != b
    assert databank.cdl(0.1, 0.03, 0.12) == databank.cdl(0.1, 0.03, 0.12, fix_x8=True)
