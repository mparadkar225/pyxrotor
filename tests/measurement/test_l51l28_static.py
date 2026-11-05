"""XROTOR's compressibility model against static tests of three two-blade NACA propellers
with thin NACA 16-series sections, NACA RM L51L28 (Wood & Swihart 1952;
validation/naca_l51l28): the 10-(0)(066)-03, 10-(3)(066)-03 and 10-(5)(066)-03, design CL
0, 0.3 and 0.5, at blade angles 0-16 deg and tip Mach numbers up to 1.02.

The comparison is at the test's thrust: the model's blade angle is chosen to give the test's
CT, and its power is compared. At the test's blade angle the test's CT and CP rise 15-35% more
with tip Mach number than the model's, which the report puts partly on blade twist under load.
What these tests pin:
  * the digitized curves are physical (figure of merit <= 1, CT and CP grow with blade angle
    and with camber);
  * at low tip Mach number the model needs 25-42% too little power at 12 deg and on the
    symmetric blade at 8 deg (the inboard blade is stalled at zero advance, and XROTOR's post-stall drag,
    2 (excess angle)^2, is far below a stalled section's);
  * with 7.55's model and Mcrit = 0.85 - 1.2 t, the power at the test's thrust rises 1.4-2.8x
    too much from Mt 0.8 to 1.0, and already 1.1-1.6x at Mt 0.9;
  * with Mcrit + 0.20 the rise is within 10% of the test's in every case but the symmetric
    blade at 12 deg (+26%); Mcrit + 0.15, the P-47 correction, is within 8% for the cambered
    blades at 8-12 deg and 14-52% high where the sections work far from their design CL.
"""
import csv
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
HERE = Path(__file__).resolve().parents[2] / "validation" / "naca_l51l28"
_spec = importlib.util.spec_from_file_location("naca_l51l28_compare", HERE / "compare.py")
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)
model = compare.model

CASES = [(0.0, 4.0), (0.0, 8.0), (0.0, 12.0), (0.3, 4.0), (0.3, 8.0), (0.3, 12.0),
         (0.5, 4.0), (0.5, 8.0), (0.5, 12.0)]


@pytest.fixture(scope="module")
def matched():
    compare.SETTINGS = {k: compare.SETTINGS[k]
                        for k in ("7.55", "mcrit+0.15", "mcrit+0.20", "incompressible")}
    return compare.at_test_thrust()


def test_digitized_curves_are_physical():
    """Figure of merit CT^1.5 / (sqrt(pi/2) CP) below 0.8 (it is up to 0.68); CT and CP grow with blade angle at
    each tip Mach number, and CT with design CL at each blade angle."""
    rows = model.tests()
    for r in rows:
        ct, cp = float(r["CT"]), float(r["CP"])
        assert ct ** 1.5 / (math.sqrt(math.pi / 2.0) * cp) < 0.8, r
        assert abs(ct / cp / float(r["CT_over_CP"]) - 1.0) < 0.01     # CT, CP rounded to 4 places
    key = {(r["design_cl"], float(r["beta75"]), r["Mt"]): r for r in rows}
    for (cld, beta, mt), r in key.items():
        up = key.get((cld, beta + 4.0, mt))
        if up:
            assert float(up["CT"]) > float(r["CT"]) and float(up["CP"]) > float(r["CP"])
        more = key.get(({"0.0": "0.3", "0.3": "0.5"}.get(cld), beta, mt))
        if more:
            assert float(more["CT"]) > float(r["CT"])


def test_static_limit():
    """J_STATIC = 0.002 is the static limit to 0.3%: halving it changes CT and CP by < 0.2%."""
    for cld, beta, mt in ((0.3, 8.0, 0.6), (0.0, 12.0, 0.9), (0.5, 4.0, 1.0)):
        a = model.static(cld, beta, mt)
        b = model.static(cld, beta, mt, j=0.5 * model.J_STATIC)
        assert a["converged"] and b["converged"]
        assert abs(a["CT"] / b["CT"] - 1.0) < 0.002 and abs(a["CP"] / b["CP"] - 1.0) < 0.002


def ratio(matched, cld, beta, mt, name="incompressible"):
    r = [r for r in matched if (r["design_cl"], r["beta75"], r["Mt"]) == (cld, beta, mt)][0]
    return r[name]["CP"] / r["CP_test"]


def test_static_power_too_low_at_high_blade_angles(matched):
    """At Mt 0.4 the model needs too little power for the test's thrust: 38-42% on the
    symmetric blade at 8 and 12 deg, 25-31% on the cambered ones at 12 deg; within 6% at
    4 deg on the cambered blades."""
    for cld, beta in CASES:
        r = ratio(matched, cld, beta, 0.4)
        if beta >= 12.0 or (beta >= 8.0 and cld == 0.0):
            assert r < 0.85, (cld, beta, r)
        if beta == 4.0 and cld > 0.0:
            assert 0.9 < r < 1.0, (cld, beta, r)


def test_755_power_rise_far_too_large(matched):
    """7.55 with Mcrit = 0.85 - 1.2 t: from Mt 0.8 the power at the test's thrust grows 1.4x
    or more too much by Mt 1.0 (or no solution reaches the test's thrust), and 1.1x by 0.9."""
    for cld, beta in CASES:
        d = compare.drift(matched, cld, beta, "7.55")
        assert d[0.9] > 1.09, (cld, beta, d)
        assert 1.0 not in d or d[1.0] > 1.35, (cld, beta, d)


def test_mcrit_shift(matched):
    """Mcrit + 0.20: the power rise from Mt 0.8 to 1.0 within 10% of the test's except on the
    symmetric blade at 12 deg (+26%). Mcrit + 0.15 (the P-47 correction): within 8% for the
    cambered blades at 8 and 12 deg, too high at 4 deg (local CL far below the design CL,
    where XROTOR lowers Mcrit by 0.25 |CL - CLDMIN|) and on the symmetric blade."""
    for cld, beta in CASES:
        d20 = compare.drift(matched, cld, beta, "mcrit+0.20")
        d15 = compare.drift(matched, cld, beta, "mcrit+0.15")
        if (cld, beta) == (0.0, 12.0):
            assert 1.15 < d20[1.0] < 1.35 and d15[1.0] > 1.4
            continue
        assert abs(d20[1.0] - 1.0) < 0.1, (cld, beta, d20)
        if cld > 0.0 and beta >= 8.0:
            assert abs(d15[1.0] - 1.0) < 0.08, (cld, beta, d15)
        if beta == 4.0:
            assert d15[1.0] > 1.1, (cld, beta, d15)


def test_incompressible_misses_the_low_angle_rise(matched):
    """Without compressibility the power at the test's thrust falls 3-12% behind the test's by
    Mt 1.0 at 4 deg: the test has a compressibility power rise there that the model with
    Mcrit + 0.20 matches."""
    for cld, beta in [(c, b) for c, b in CASES if b == 4.0]:
        d = compare.drift(matched, cld, beta, "incompressible")
        assert d[1.0] < 0.98, (cld, d)
