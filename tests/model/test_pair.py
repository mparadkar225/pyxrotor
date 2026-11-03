"""The contra-rotating pair solver (pyxrotor.pair): limits it must reach exactly or nearly,
and its operating modes' targets. Each limit ties the pair back to something independent:
the isolated XROTOR rotor, XROTOR's own VPUT/VGET procedure, a single rotor with twice the
blades, and angular-momentum conservation."""
import numpy as np
import pytest

from pyxrotor.constants import PI
from pyxrotor.pair import PairConvergenceError, RotorPair
from test_rotor import rotol_both

pytestmark = pytest.mark.model
RPM = 1100.0


def rotol(theta=22.0, nblds=2, model=None):
    _, r = rotol_both(theta)
    if nblds != 2:
        r.s.nblds = nblds
    if model:
        r.design_formulation(model)
    return r


@pytest.fixture(scope="module")
def isolated():
    r = rotol()
    assert r.solve_rpm(RPM)
    return r


def rel(a, b):
    return abs(a / b - 1.0)


def test_far_apart_front_is_isolated_and_rear_is_xrotors_vget(isolated):
    """At 50 R spacing the rear's upstream field is ~1e-4 of its disk value: the front is
    the isolated rotor, and the rear sees the front's fully developed slipstream, which is
    exactly what XROTOR's VPUT then VGET (weights 1, -1) imposes."""
    R = isolated.s.rad
    res = RotorPair(rotol(), rotol(), 50.0 * R).solve_rpm(RPM)
    assert res.converged
    iso = isolated.performance()
    for k in ("thrust", "power", "torque"):
        assert rel(res.front[k], iso[k]) < 2e-4, k
    # XROTOR's procedure by hand
    p = isolated.slipstream_profile()
    rear = rotol()
    rear.set_slipstream(p["r"], p["u"], p["v"], axial_weight=1.0, swirl_weight=-1.0)
    assert rear.solve_rpm(RPM)
    hand = rear.performance()
    for k in ("thrust", "power", "torque"):
        assert rel(res.rear[k], hand[k]) < 1e-3, k


def test_unloaded_rear_leaves_front_isolated(isolated):
    """A rear rotor carrying almost no circulation (chord x 1e-5) changes nothing upstream."""
    rear = rotol()
    rear.s.ch[: rear.s.ii] *= 1e-5
    res = RotorPair(rotol(), rear, 0.1 * isolated.s.rad).solve_rpm(RPM)
    iso = isolated.performance()
    for k in ("thrust", "power", "torque"):
        assert rel(res.front[k], iso[k]) < 1e-4, k
    assert abs(res.rear["thrust"]) < 1e-3 * iso["thrust"]


@pytest.mark.parametrize("model", ["pot", "grad"])
def test_coincident_corotating_rows_approach_twice_the_blades(model):
    """Two identical 2-blade rows co-rotating at (almost) zero spacing are nearly one
    4-blade rotor. Not exact: each row's own solve applies a 2-blade tip loss (stronger than
    a 4-blade one), which makes the pair 2-3% low."""
    four = rotol(nblds=4, model=model)
    assert four.solve_rpm(RPM)
    p4 = four.performance()
    res = RotorPair(rotol(model=model), rotol(model=model), 1e-6 * four.s.rad,
                    counter_rotating=False).solve_rpm(RPM)
    assert res.converged
    assert -0.04 < res.thrust / p4["thrust"] - 1.0 < 0.0
    assert -0.04 < res.power / p4["power"] - 1.0 < 0.0


def test_equal_torque_counter_rotation_recovers_the_swirl(isolated):
    """With equal and opposite torques the pair's far wake carries no net angular momentum,
    and its swirl kinetic-energy loss nearly vanishes; co-rotating, the swirl adds up."""
    R = isolated.s.rad
    s = isolated.s
    f = isolated.induced_field()
    r = np.linspace(0, R, 401)[1:]
    u, _, v = f.velocity(r, 1e3 * R)
    iso_loss = np.sum(0.5 * s.rho * v * v * (s.vel + u) * 2 * PI * r * (r[1] - r[0])) / isolated.performance()["power"]
    assert 0.005 < iso_loss < 0.02

    cr = RotorPair(rotol(), rotol(), 0.2 * R).solve_power(170e3, RPM, torque_split=0.5)
    assert cr.converged and abs(cr.net_torque) < 1e-9 * abs(cr.front["torque"])
    assert cr.swirl_power_fraction < 0.05 * iso_loss
    # angular momentum flux of the far wake: the residual swirl is small and changes sign
    rr, vv = cr.residual_swirl["r"], cr.residual_swirl["v"]
    flux = np.sum(vv * rr * rr)
    ref = np.sum(np.abs(v) * r * r)
    assert abs(flux) < 0.05 * ref

    co = RotorPair(rotol(), rotol(), 0.2 * R, counter_rotating=False).solve_rpm(RPM)
    assert co.swirl_power_fraction > iso_loss
    assert cr.efficiency > co.efficiency + 0.01


def test_operating_modes_meet_their_targets(isolated):
    R = isolated.s.rad
    p = RotorPair(rotol(), rotol(), 0.2 * R).solve_power(170e3, RPM, torque_split=0.6)
    assert p.converged and rel(p.power, 170e3) < 1e-9 and abs(p.torque_split - 0.6) < 1e-9
    assert p.front["rpm"] == pytest.approx(RPM) and p.rear["rpm"] == pytest.approx(RPM)

    t = RotorPair(rotol(), rotol(), 0.2 * R).solve_thrust(2800.0, RPM, torque_split=0.5)
    assert t.converged and rel(t.thrust, 2800.0) < 1e-6 and abs(t.torque_split - 0.5) < 1e-9

    fp = RotorPair(rotol(), rotol(), 0.2 * R).solve_power_fixed_pitch(170e3, RPM)
    assert fp.converged and rel(fp.power, 170e3) < 1e-6
    assert fp.front["rpm"] == pytest.approx(fp.rear["rpm"])

    g = RotorPair(rotol(), rotol(18.0), 0.2 * R, gear_ratio=1.2).solve_rpm(RPM)
    assert g.converged and g.rear["rpm"] == pytest.approx(1.2 * RPM)


def test_iteration_converges_geometrically(isolated):
    res = RotorPair(rotol(), rotol(), 0.2 * isolated.s.rad).solve_rpm(RPM)
    h = np.array(res.history)
    assert res.converged and h[-1] < 1e-6
    assert np.all(h[3:] < h[2:-1])                  # monotone after the first passes


def test_front_unloads_and_pair_field_is_the_sum(isolated):
    """The rear draws extra flow through the front, which unloads it; the pair's field is
    the sum of the two rotors' fields (swirl in the front's direction)."""
    R = isolated.s.rad
    pair = RotorPair(rotol(), rotol(), 0.2 * R)
    res = pair.solve_rpm(RPM)
    assert res.front["thrust"] < isolated.performance()["thrust"]
    r, z = 0.6 * R, 0.5 * R
    u, w, v = pair.field(r, z)
    a = pair.front.induced_field().velocity(r, z)
    b = pair.rear.induced_field().velocity(r, z - 0.2 * R)
    assert u == pytest.approx(a[0] + b[0]) and v == pytest.approx(a[2] - b[2])


def test_mismatched_flight_conditions_rejected():
    a, b = rotol(), rotol()
    b.velocity(60.0)
    with pytest.raises(ValueError):
        RotorPair(a, b, 0.1)
    with pytest.raises(ValueError):
        RotorPair(rotol(), rotol(), -0.1)


def test_nonconvergence_is_reported():
    pair = RotorPair(rotol(), rotol(), 0.2, max_iter=2)
    res = pair.solve_rpm(RPM)
    assert not res.converged and res.iterations == 2
    with pytest.raises(PairConvergenceError):
        RotorPair(rotol(), rotol(), 0.2).solve_rpm(1.0)     # a hopeless operating point
