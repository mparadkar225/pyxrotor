"""Realistic solver states for differential tests.

A converged operating point of the Rotol Spitfire propeller (ARC R&M 2357 geometry) is
computed once in the reference. Tests then perturb it with Hypothesis so routines see
varied but physically shaped inputs: induced velocities, sensitivities, ducts, added
inflow, flags.
"""
from __future__ import annotations

import functools
import json
from pathlib import Path

import numpy as np
from hypothesis import strategies as st

from pyxrotor.constants import PI
from pyxrotor.state import State
from reference.ref755 import Ref755
from statesync import state_from_ref
from strategies import finite

ROOT = Path(__file__).resolve().parents[2]
RADIUS = 10.75 * 0.3048 / 2
SECTIONS = [  # r/R, a0 (deg), dCL/da, CLmax, CLmin, CDmin, CL@CDmin, dCD/dCL^2, Mcrit
    (0.0, -5.0, 5.5, 1.10, -0.5, 0.020, 0.30, 0.010, 0.55),
    (0.4, -4.5, 5.9, 1.30, -0.6, 0.011, 0.35, 0.008, 0.66),
    (0.7, -3.5, 6.1, 1.35, -0.6, 0.0085, 0.40, 0.007, 0.74),
    (1.0, -2.5, 6.2, 1.30, -0.5, 0.0075, 0.35, 0.007, 0.79),
]


def geometry():
    return json.loads((ROOT / "validation" / "rm2357" / "geometry.json").read_text())["stations"]


def rotol_ref(theta07: float = 22.0, ii: int | None = None, nblds: int = 2) -> Ref755:
    """A reference instance with the Rotol blade loaded through ARBI (not yet solved)."""
    ref = Ref755()
    if ii is not None:
        ref["II"] = ii
        ref["IINF"] = ii + ii // 2
    _, vso, rho, rmu = ref.call("ATMO", 0.0, 0.0, 0.0, 0.0)
    ref["ALT"], ref["VSO"], ref["RHO"], ref["RMU"] = 0.0, float(vso), float(rho), float(rmu)
    for n, (x, a0, dcla, clmax, clmin, cdmin, cldmin, dcdcl2, mcrit) in enumerate(SECTIONS, 1):
        ref.call("PUTAERO", n, x, a0 * PI / 180, clmax, clmin, dcla, 0.1, 0.1,
                 cdmin, cldmin, dcdcl2, -0.1, mcrit, 2.0e6, -0.15)
    ref["NAERO"] = len(SECTIONS)
    geo = geometry()
    ref.call("REFARBI", nblds, 50.0, RADIUS, 0.2 * RADIUS, len(geo),
             np.array([g["r"] for g in geo]), np.array([g["cR"] for g in geo]),
             np.array([theta07 + g["dth"] for g in geo]))
    return ref


@functools.lru_cache(maxsize=None)
def _solved(theta07: float, j: float, mtip: float) -> State:
    ref = rotol_ref(theta07)
    n = mtip * ref.get("VSO") / (2 * PI * RADIUS)
    vel = j * n * 2 * RADIUS
    ref["VEL"] = vel
    ref["ADV"] = vel / (RADIUS * n * 60 * PI / 30.0)
    ref["CONV"] = False
    ref.call("APER", 4, 2, True)
    assert ref.get("CONV") == 1
    return state_from_ref(ref)


def solved_state(theta07: float = 22.0, j: float = 0.55, mtip: float = 0.41) -> State:
    """A fresh copy of a converged operating point."""
    return _solved(theta07, j, mtip).copy()


@st.composite
def perturbed_states(draw, duct=True, added_inflow=True):
    """A converged state with its velocities, sensitivities and options perturbed."""
    s = solved_state(draw(st.sampled_from([15.0, 22.0, 27.0])),
                     draw(st.sampled_from([0.3, 0.55, 0.8])))
    ii = s.ii
    rng = np.random.default_rng(draw(st.integers(0, 2**32 - 1)))
    amp = draw(finite(0.0, 0.5))
    s.adv *= draw(finite(0.4, 2.0))
    s.adw = s.adv * draw(finite(0.8, 1.4))
    for arr in (s.vind, s.vind_adw, s.vind_gam):
        arr *= 1.0 + amp * rng.standard_normal(arr.shape)
    if duct and draw(st.booleans()):
        s.duct = True
        s.urduct = draw(finite(0.7, 1.4))
    if draw(st.booleans()):
        s.ubody[:ii] = draw(finite(-0.1, 0.1)) * rng.standard_normal(ii)
    if added_inflow and draw(st.booleans()):
        n = draw(st.integers(2, 12))
        s.nadd = n
        r = np.sort(rng.uniform(0.0, 1.2 * s.rad, n))
        r[0] = min(r[0], 0.1 * s.rad)
        s.radd[:n] = r + np.arange(n) * 1e-3
        for arr in (s.uadd, s.vadd, s.uaddr, s.vaddr):
            arr[:n] = rng.normal(0.0, 0.2 * s.vel, n)
    return s
