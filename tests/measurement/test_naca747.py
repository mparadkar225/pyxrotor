"""NACA Report 747: four-blade single- and dual-rotating 10-ft propellers in the Langley
propeller-research tunnel (Hamilton Standard 3155-6 / 3156-6 blades, Clark Y sections).

The tested "single" propeller is two 2-blade hubs in tandem, co-rotating 9.7 in apart; the
dual one is the same pair counter-rotating, with the rear blades set 1.1-1.6 deg lower.
Both are modelled as pyxrotor.pair.RotorPair (validation/naca747/model.py), so the dual
vs. single difference isolates the counter-rotation effect, while the estimated Clark Y
polars largely cancel. The full comparison (220 points, ~10 min) is
validation/naca747/compare.py + report.py; this file checks the data and a few points.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
DATA = Path(__file__).resolve().parents[2] / "validation" / "naca747"
sys.path.insert(0, str(DATA))
import model  # noqa: E402

TESTS = json.loads((DATA / "tests.json").read_text())
CURVES = {(c["rotation"], c["quantity"], c["beta_front"]): c for c in TESTS["curves"]}


def measured(rot, q, beta, J):
    c = CURVES[(rot, q, beta)]
    assert c["J"][0] - 1e-9 <= J <= c["J"][-1] + 1e-9
    return float(np.interp(J, c["J"], c["value"]))


@pytest.mark.parametrize("rot", ["single", "dual"])
def test_digitized_curves_agree_with_the_reports_efficiency_curves(rot):
    """eta = CT J / CP from the digitized CT and CP figures never exceeds the report's own
    efficiency curves (figures 7, 10) by more than reading error, and comes close to
    their peak. Near zero thrust a 0.002 reading error is several percent, so CT > 0.07.
    This caught three CP crops whose lowest grid line was not the axis."""
    peaks = TESTS["peak_efficiency"][rot]
    for b in (20, 25, 30, 35, 40, 45):
        t, p = CURVES[(rot, "CT", b)], CURVES[(rot, "CP", b)]
        J = np.linspace(max(t["J"][0], p["J"][0]), min(t["J"][-1], p["J"][-1]), 25)
        ct = np.interp(J, t["J"], t["value"])
        eta = (ct * J / np.interp(J, p["J"], p["value"]))[ct > 0.07]   # reading error grows near 0
        assert eta.max() < peaks[str(b)] + 0.03, b
        assert eta.max() > peaks[str(b)] - 0.06, b


def test_operating_points_follow_the_test_schedule():
    v, rpm = model.operating_point(0.5)
    assert rpm == pytest.approx(550.0) and v == pytest.approx(0.5 * 550 / 60 * model.D)
    v, rpm = model.operating_point(2.5)
    assert v == pytest.approx(110 * 0.44704) and rpm < 550.0


@pytest.fixture(scope="module")
def points():
    """single and dual at two blade angles, near peak efficiency (attached flow)."""
    out = {}
    for beta, J in ((25, 0.9), (45, 2.2)):
        out[("single", beta)] = model.single(beta, J)
        out[("dual", beta)] = model.dual(beta, J)
        out[("J", beta)] = J
    return out


@pytest.mark.parametrize("beta", [25, 45])
def test_single_rotation_thrust_and_power(points, beta):
    """The co-rotating tandem pair: CT and CP within 8% near peak efficiency."""
    J, r = points[("J", beta)], points[("single", beta)]
    assert r["converged"]
    assert abs(r["CT"] / measured("single", "CT", beta, J) - 1) < 0.08
    assert abs(r["CP"] / measured("single", "CP", beta, J) - 1) < 0.08


@pytest.mark.parametrize("beta", [25, 45])
def test_counter_rotation_efficiency_gain(points, beta):
    """Dual minus single efficiency at the same front angle and J: measured +0.008 (25 deg)
    and +0.023 (45 deg); the model gives +0.015 and +0.020."""
    J = points[("J", beta)]
    meas = {rot: measured(rot, "CT", beta, J) * J / measured(rot, "CP", beta, J)
            for rot in ("single", "dual")}
    gain_m = meas["dual"] - meas["single"]
    gain_p = points[("dual", beta)]["eta"] - points[("single", beta)]["eta"]
    assert gain_p > 0
    assert abs(gain_p - gain_m) < 0.012


def test_known_discrepancy_dual_loading(points):
    """Known: at equal nominal blade angles the model's dual propeller is more heavily loaded
    than the tested one (CT 10-20% high near zero thrust, about 0.6 deg of effective
    pitch), while the single one matches. Figure 12 shows the deficit is mostly the front
    rotor's. Pinned so a model change that moves it is noticed."""
    J, r = points[("J", 45)], points[("dual", 45)]
    ratio = r["CT"] / measured("dual", "CT", 45, J)
    assert 1.10 < ratio < 1.25
    # figure 12: the excess is mostly the front rotor's (model 1.20x, rear 1.09x at J = 2.2)
    front = r["CP_front"] / measured("dual", "CP_front", 45, J)
    rear = r["CP_rear"] / measured("dual", "CP_rear", 45, J)
    assert 1.12 < front < 1.28 and 1.02 < rear < 1.15 and front > rear + 0.05


def test_individual_powers_add_up_to_the_total():
    """Figure 12's front and rear CP (digitized separately: solid and dashed lines) sum to
    figure 13's total, which was digitized independently."""
    for b in (25, 30, 35, 40, 45, 50):
        f, r, t = (CURVES[("dual", q, b)] for q in ("CP_front", "CP_rear", "CP"))
        J = np.array([j for j in r["J"] if t["J"][0] <= j <= t["J"][-1] and f["J"][0] <= j <= f["J"][-1]])
        assert len(J) >= 5, b
        err = (np.interp(J, f["J"], f["value"]) + np.interp(J, r["J"], r["value"])
               - np.interp(J, t["J"], t["value"]))
        assert abs(err.mean()) < 0.003 and abs(err).max() < 0.012, b   # worst: a curve end
