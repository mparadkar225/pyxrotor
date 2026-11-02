"""Hypothesis strategies and reference set-up shared by the differential tests."""
from __future__ import annotations

import numpy as np
from hypothesis import strategies as st

from pyxrotor.aero import AERODATA_ROW, AeroSection


def finite(lo, hi):
    return st.floats(min_value=lo, max_value=hi, allow_nan=False, allow_infinity=False)


@st.composite
def sections_params(draw, xisect=None):
    """Section parameters spanning realistic and somewhat extreme airfoils."""
    return AeroSection(
        xisect=draw(finite(0.0, 1.0)) if xisect is None else xisect,
        a0=draw(finite(-0.15, 0.05)),
        clmax=draw(finite(0.6, 2.0)),
        clmin=draw(finite(-1.5, -0.1)),
        dclda=draw(finite(3.0, 7.0)),
        dclda_stall=draw(finite(0.0, 0.6)),
        dcl_stall=draw(finite(0.02, 0.5)),
        cdmin=draw(finite(0.004, 0.03)),
        cldmin=draw(finite(-0.2, 0.8)),
        dcdcl2=draw(finite(0.0, 0.02)),
        cmcon=draw(finite(-0.2, 0.05)),
        reref=draw(finite(1e5, 1e7)),
        rexp=draw(finite(-1.0, 0.0)),
        mcrit=draw(finite(0.5, 0.9)),
    )


@st.composite
def section_sets(draw, max_sections=5):
    """1 to max_sections sections, sorted by r/R, the first at r/R = 0 (as SETDEF leaves it)."""
    n = draw(st.integers(min_value=1, max_value=max_sections))
    xs = sorted([0.0] + draw(st.lists(finite(0.05, 1.0), min_size=n - 1, max_size=n - 1, unique=True)))
    return [draw(sections_params(xisect=x)) for x in xs]


# Angle of attack well past both stalls; local speed giving local Mach 0 to ~1.3 at VEL/VSO
# in [0.05, 0.6] (so the Mach-limiting branch is exercised); Re off (0) or 1e4..1e7.
alphas = finite(-0.8, 0.8)
speeds = finite(0.0, 6.0)
reynolds = st.one_of(st.just(0.0), finite(1e4, 1e7))
flight = st.tuples(finite(20.0, 200.0), finite(300.0, 400.0))      # (VEL, VSO)


def load_sections(ref, sections):
    """Write sections into the reference's XIAERO/AERODATA/NAERO, as PUTAERO would."""
    aerodata, xiaero = ref["AERODATA"], ref["XIAERO"]
    for n, sec in enumerate(sections):
        xiaero[n] = sec.xisect
        for name, row in AERODATA_ROW.items():
            aerodata[row - 1, n] = getattr(sec, name)
    ref["NAERO"] = len(sections)


def iaero_to_fortran(iaero):
    return np.asarray(iaero) + 1


def iaero_from_fortran(iaero):
    return np.asarray(iaero, dtype=np.int64) - 1
