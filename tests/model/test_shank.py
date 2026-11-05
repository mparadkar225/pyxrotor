"""The round-shank drag model (pyxrotor.shank): its per-station forces are XROTOR's own
viscous thrust and power formulas, the placeholder section really has no lift or drag in
XROTOR, the blend weights follow GETCLCDCM's interpolation, and the constant-power solve
absorbs the requested power."""
import numpy as np
import pytest

from pyxrotor import shank
from pyxrotor.aero import clcdcm
from test_rotor import rotol_both

pytestmark = pytest.mark.model


@pytest.fixture
def solved():
    _, r = rotol_both(22.0)
    assert r.solve_rpm(1100.0)
    return r


def test_station_forces_are_tpqs(solved):
    """With XROTOR's own CD at every station, the sums are TPQ's TVIS and PVIS, bit for bit."""
    s = solved.s
    sh = shank.RoundShank(0.3, 0.4)
    dtv, dpv, mach = sh.station_drag(solved, cd=s.cd[:s.ii])
    assert np.add.accumulate([0.0, *dtv])[-1] == s.tvis
    assert np.add.accumulate([0.0, *dpv])[-1] == s.pvis
    assert np.all((mach > 0.1) & (mach < 0.6))


def test_placeholder_has_no_lift_or_drag():
    sec = shank.round_section(0.25)
    for alf in (-0.5, 0.0, 0.3, 1.0):
        for m in (0.2, 0.6, 0.95):
            a = clcdcm(alf, m / 0.5, 1.0e6, 0.5, 1.0, sec)
            assert abs(a.clift) <= 0.01 * abs(alf) * 10.1 and a.cdrag == 0.0


def test_cylinder_curve():
    """TR-852 force tests: 0.51-0.58 below Mach 0.5, peak 1.09 at 0.70, held outside."""
    assert shank.cylinder_cd(0.1) == shank.cylinder_cd(0.30) == 0.58
    assert shank.cylinder_cd(0.70) == 1.09 and shank.cylinder_cd(0.95) == 0.86
    assert 0.5 < shank.cylinder_cd(0.45) < 0.53


def _rotol_with_shank(sh):
    _, r = rotol_both(22.0)
    secs = shank_sections(r, sh)
    r.set_sections(secs)
    return r


def shank_sections(r, sh):
    from pyxrotor.rotor import sections_of
    secs = [x for x in sections_of(r.s) if x.xisect >= sh.r_airfoil]
    if not secs or secs[0].xisect > sh.r_airfoil:
        first = sections_of(r.s)[-1] if not secs else secs[0]
        secs.insert(0, type(first)(**{**first.__dict__, "xisect": sh.r_airfoil}))
    return [sh.round_section(float(r.s.xi[0]) - 0.01), sh.round_section(sh.r_round), *secs]


def test_weights_follow_the_section_blend():
    sh = shank.RoundShank(0.30, 0.40)
    r = _rotol_with_shank(sh)
    w = sh.weights(r)
    xi = r.s.xi[:r.s.ii]
    assert np.all(w[xi <= 0.30] == 1.0) and np.all(w[xi >= 0.40] == 0.0)
    mid = (xi > 0.30) & (xi < 0.40)
    assert np.allclose(w[mid], (0.40 - xi[mid]) / 0.10, rtol=0, atol=1e-12)


def test_round_stations_carry_no_lift_and_the_power_is_absorbed():
    sh = shank.RoundShank(0.30, 0.40)
    r = _rotol_with_shank(sh)
    assert r.solve_rpm(1100.0)
    s = r.s
    round_ = s.xi[:s.ii] <= 0.30
    assert np.all(np.abs(s.cl[:s.ii][round_]) < 0.01) and np.all(s.cd[:s.ii][round_] == 0.0)
    d = sh.drag(r)
    assert d["thrust"] < 0.0 < d["power"]
    target = 1.2 * r.performance()["power"]
    sol = sh.solve_power(r, target, 1100.0)
    assert sol.usable and sol.exact and abs(sol.power / target - 1.0) < 1e-4
    p = sh.performance(r)
    assert abs(p["power"] / target - 1.0) < 1e-4
    assert p["thrust"] == pytest.approx(p["thrust_blades"] + sh.drag(r)["thrust"])
    assert abs(sol.efficiency - p["efficiency"]) < 1e-9


def test_oval_shank_is_drag_on_its_thickness(solved):
    """thickness_ratio h/b scales the cylinder drag on the chord by h/b."""
    sh1 = shank.RoundShank(0.30, 0.40)
    r = _rotol_with_shank(sh1)
    assert r.solve_rpm(1100.0)
    sh2 = shank.RoundShank(0.30, 0.40, thickness_ratio=([0.0, 1.0], [0.5, 0.5]))
    assert sh2.drag(r)["power"] == pytest.approx(0.5 * sh1.drag(r)["power"], rel=1e-12)
