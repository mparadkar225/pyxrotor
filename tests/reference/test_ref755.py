"""Checks on the XROTOR 7.55 reference library itself.

These establish that the library computes what the interactive program does, so that
differential tests against it mean something.
"""
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pyxrotor.constants import PI
from reference import ref755

pytestmark = pytest.mark.reference
ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.skipif(shutil.which("nm") is None, reason="needs nm")
def test_common_layout_matches_compiled_library(ref):
    """Every COMMON block parsed from XROTOR.INC has exactly the size the compiler gave it."""
    out = subprocess.run(["nm", "-S", str(ref755.LIB)], capture_output=True, text=True, check=True).stdout
    sizes = {f[3]: int(f[1], 16) for f in (line.split() for line in out.splitlines()) if len(f) == 4}
    for blk in ref.blocks.values():
        assert sizes.get(blk.symbol) == blk.size, blk.name
    assert len(ref.blocks) == 27


def test_startup_state(ref):
    assert ref.get("VERSION") == pytest.approx(7.55)
    assert ref.get("II") == 30 and ref.get("IXSPAC") == 2
    assert ref.get("NITERD") == 40 and ref.get("NITERA") == 40
    assert ref.get("XINF") == 3.0 and ref.get("GEE") == 9.81
    assert ref.get("NAERO") == 1


def test_instances_are_isolated(ref):
    other = ref755.Ref755()
    other["II"] = 7
    assert ref.get("II") == 30 and other.get("II") == 7


def test_atmosphere_matches_interactive_program(ref):
    """`atmo 7.62` in the interactive program printed 309.520 m/s, 0.54994 kg/m^3, 0.1530E-04."""
    _, vso, rho, rmu = ref.call("ATMO", 7.62, 0.0, 0.0, 0.0)
    assert float(vso) == pytest.approx(309.520, abs=5e-4)
    assert float(rho) == pytest.approx(0.54994, abs=5e-6)
    assert float(rmu) == pytest.approx(1.530e-5, abs=5e-9)


# --- A full analysis through the library, against the interactive program ---------------
# The Rotol Spitfire propeller from ARC R&M 2357, set up exactly as in the interactive runs
# of Nov 2026: four parametric sections, ARBI geometry, then `velo` + `rpm` in OPER.
SECTIONS = [  # r/R, a0 (deg), dCL/da, CLmax, CLmin, CDmin, CL@CDmin, dCD/dCL^2, Mcrit
    (0.0, -5.0, 5.5, 1.10, -0.5, 0.020, 0.30, 0.010, 0.55),
    (0.4, -4.5, 5.9, 1.30, -0.6, 0.011, 0.35, 0.008, 0.66),
    (0.7, -3.5, 6.1, 1.35, -0.6, 0.0085, 0.40, 0.007, 0.74),
    (1.0, -2.5, 6.2, 1.30, -0.5, 0.0075, 0.35, 0.007, 0.79),
]
# Printed by the interactive program at theta(0.7R) = 22 deg: (V m/s, rpm, Ct, Cp, eff)
INTERACTIVE = [(18.207, 805.301, 0.0966, 0.06013, 0.6652)]


def rotol_setup(ref, theta07):
    geo = json.loads((ROOT / "validation" / "rm2357" / "geometry.json").read_text())["stations"]
    radius = 10.75 * 0.3048 / 2
    _, vso, rho, rmu = ref.call("ATMO", 0.0, 0.0, 0.0, 0.0)
    ref["ALT"], ref["VSO"], ref["RHO"], ref["RMU"] = 0.0, float(vso), float(rho), float(rmu)
    for n, (x, a0, dcla, clmax, clmin, cdmin, cldmin, dcdcl2, mcrit) in enumerate(SECTIONS, start=1):
        ref.call("PUTAERO", n, x, a0 * PI / 180, clmax, clmin, dcla, 0.1, 0.1,
                 cdmin, cldmin, dcdcl2, -0.1, mcrit, 2.0e6, -0.15)
    ref["NAERO"] = len(SECTIONS)
    ref.call("REFARBI", 2, 50.0, radius, 0.2 * radius, len(geo),
             np.array([g["r"] for g in geo]), np.array([g["cR"] for g in geo]),
             np.array([theta07 + g["dth"] for g in geo]))
    return radius


@pytest.mark.parametrize("vel,rpm,ct,cp,eff", INTERACTIVE)
def test_rotol_operating_point_matches_interactive_program(ref, vel, rpm, ct, cp, eff):
    radius = rotol_setup(ref, 22.0)
    ref["VEL"] = vel                                     # OPER `velo`
    ref["ADV"] = vel / (radius * rpm * PI / 30.0)       # OPER `rpm`
    ref["CONV"] = False
    ref.call("APER", 4, 2, True)
    assert ref.get("CONV") == 1

    rho, n, d = ref.get("RHO"), rpm / 60.0, 2 * radius
    thrust = ref.get("TTOT") * rho * vel**2 * radius**2
    power = ref.get("PTOT") * rho * vel**3 * radius**2
    # The interactive build kept DOUBLE PRECISION at 16 bytes; this one uses 8, and the
    # printed values have 4-5 digits, so agreement is to print precision, not bitwise.
    assert thrust / (rho * n**2 * d**4) == pytest.approx(ct, rel=1e-3)
    assert power / (rho * n**3 * d**5) == pytest.approx(cp, rel=1e-3)
    assert ref.get("TTOT") / ref.get("PTOT") == pytest.approx(eff, abs=2e-4)
