"""The swept-blade model in the compiled core (pyxrotor.sweep), by limits and independent
checks:
  * no sweep, or a Sweep with zero offsets, is XROTOR's straight blade (bit for bit with
    sweep None, which the differential tests check everywhere; within rounding with zero
    offsets);
  * the effective sweep CSCALC computes equals the angle between the 3-D quarter-chord curve,
    built in Cartesian coordinates and differentiated numerically, and the local relative
    velocity, on a twisted blade with sweep and lean;
  * GETCLCDCM on a swept station is the section evaluated in the normal plane: CL = CLn cos^2
    at tan(an) = tan(a)/cos and Mach M cos; CD's friction part unchanged and the rest cos^3;
    the derivatives match finite differences;
  * the incompressible lift slope of a swept station is the straight one times cos L at small
    angles (simple sweep theory);
  * a converged swept solution is self-consistent: each station's stored CL is GETCLCDCM at
    its final angle, W and effective sweep;
  * DESI refuses a swept blade.
"""
import math

import numpy as np
import pytest

from pyxrotor import Rotor, aero, sections
from pyxrotor.sweep import Sweep
from pyxrotor.wake import cscalc
from pyxrotor import _speedups

PI = math.pi


def rotol_like(nblds=3):
    r = Rotor(fix_all=True)
    r.atmosphere(3.0)
    r.set_sections([sections.naca16(t, x, 0.4) for x, t in
                    ((0.2, 0.25), (0.5, 0.10), (0.75, 0.07), (1.0, 0.05))])
    xi = np.linspace(0.2, 1.0, 12)
    r.arbi(nblds=nblds, vel=120.0, tip_radius=1.7, hub_radius=0.2 * 1.7, r_over_R=xi,
           c_over_R=0.12 + 0.05 * np.sin(PI * xi), beta_deg=60.0 - 40.0 * (xi - 0.2))
    return r


def test_zero_offsets_are_straight():
    a, b = rotol_like(), rotol_like()
    b.set_sweep(Sweep([0.2, 1.0], [0.0, 0.0]))
    assert a.solve_rpm(1800.0) and b.solve_rpm(1800.0)
    pa, pb = a.performance(), b.performance()
    for k in ("thrust", "power", "efficiency"):
        assert pb[k] == pytest.approx(pa[k], rel=1e-12, abs=1e-12)
    assert np.all(np.asarray(b.s.sweep_cos[:b.s.ii]) == 1.0)


def _cartesian_sweep(r, i, h=1e-5):
    """Effective sweep at station i from the 3-D quarter-chord curve: rotating frame, blade
    along +x at theta = 0, rotation toward +theta, z downstream. The point at radius xi R is
    displaced in its section plane by s along the chord (toward the trailing edge, i.e. -theta
    and +z) and n normal to it; the tangential offset is an arc."""
    s = r.s
    sw = s.sweep

    def point(x):
        b = _beta_spline(s, x)
        sv, _, nv, _ = sw.at(x)
        y = sv * math.cos(b) - nv * math.sin(b)        # against rotation, arc length
        z = sv * math.sin(b) + nv * math.cos(b)
        th = -y / x
        return np.array([x * math.cos(th), x * math.sin(th), z])
    x = float(s.xi[i])
    t = (point(x + h) - point(x - h)) / (2 * h)
    v = cscalc(s, i)
    sv, _, nv, _ = sw.at(x)
    b = _beta_spline(s, x)
    th = -(sv * math.cos(b) - nv * math.sin(b)) / x
    e_th = np.array([-math.sin(th), math.cos(th), 0.0])
    wvec = -v.ci * e_th + np.array([0.0, 0.0, v.si])  # relative flow: against rotation
    sinl = float(np.dot(t, wvec) / (np.linalg.norm(t) * np.linalg.norm(wvec)))
    return math.sqrt(1.0 - sinl * sinl)


def _beta_spline(s, x):
    from pyxrotor.numerics import seval, spline
    xi = np.array(s.xi[: s.ii])
    be = np.array(s.beta[: s.ii])
    return seval(x, be, spline(be, xi, s.ii), xi, s.ii)


def test_effective_sweep_geometry():
    r = rotol_like()
    xi = np.linspace(0.2, 1.0, 9)
    r.set_sweep(Sweep(xi, s=0.3 * (xi - 0.5) ** 2 - 0.02, n=0.05 * (xi - 0.2)))
    assert r.solve_rpm(1800.0)
    for i in range(1, r.s.ii - 1):
        want = _cartesian_sweep(r, i)
        assert float(r.s.sweep_cos[i]) == pytest.approx(want, abs=2e-6), i


def test_from_angle():
    sw = Sweep.from_angle([0.2, 0.5, 1.0], [-10.0, 0.0, 40.0])
    assert sw.at(0.5)[0] == pytest.approx(0.0, abs=1e-12)          # crosses the axis at 0
    assert sw.at(1.0)[0] > 0 and sw.at(0.2)[0] > 0      # forward-swept root: s' < 0 inboard


def _station(cos_l, mach=0.6):
    """A one-station context: the State's section blend at station 0, with sweep_cos set."""
    r = rotol_like()
    s = r.s
    s.vel, s.vso = 100.0, 100.0 / mach
    s.sweep = Sweep([0.0, 1.0], [0.0, 0.0])
    s.sweep_cos[:] = cos_l
    return s


@pytest.mark.parametrize("cos_l", [1.0, 0.85, 0.6])
@pytest.mark.parametrize("alf", [-0.05, 0.03, 0.12, 0.4])
def test_normal_plane_section(cos_l, alf):
    s = _station(cos_l)
    i, w, rey = 5, 1.0, 2.0e6
    got = _speedups.getclcdcm(s, i, alf, w, rey)
    s.sweep = None
    alfn = math.atan2(math.sin(alf), math.cos(alf) * cos_l)
    n = _speedups.getclcdcm(s, i, alfn, w * cos_l, rey)
    # friction part of each blended section (CDMIN (Re/REREF)^REXP), blended like CD
    k = int(s.iaero[i])
    frac = (s.xi[i] - s.xiaero[k]) / (s.xiaero[k + 1] - s.xiaero[k])
    cdf = [s.aerodata[aero.AERODATA_ROW["cdmin"] - 1, j]
           * (rey / s.aerodata[aero.AERODATA_ROW["reref"] - 1, j])
           ** s.aerodata[aero.AERODATA_ROW["rexp"] - 1, j] for j in (k, k + 1)]
    cdf = (1 - frac) * cdf[0] + frac * cdf[1]
    assert got[0] == pytest.approx(n[0] * cos_l ** 2, rel=1e-13, abs=1e-15)      # CL
    assert got[7] == pytest.approx(cdf + (n[7] - cdf) * cos_l ** 3, rel=1e-12)    # CD
    # CL derivatives against finite differences (CLCDCM's CL derivatives are exact)
    s.sweep = Sweep([0.0, 1.0], [0.0, 0.0])
    h = 1e-7
    for j, (da, dw) in ((1, (h, 0)), (2, (0, h))):
        hi = _speedups.getclcdcm(s, i, alf + da, w + dw, rey)[0]
        lo = _speedups.getclcdcm(s, i, alf - da, w - dw, rey)[0]
        assert got[j] == pytest.approx((hi - lo) / (2 * h), rel=2e-5, abs=1e-8), j
    # CD derivatives by the chain rule from the normal-plane ones (CLCDCM's own CD_ALF is
    # approximate: its drag-rise term uses |dCL/da|)
    dan = cos_l / (cos_l ** 2 * math.cos(alf) ** 2 + math.sin(alf) ** 2)
    assert got[8] == pytest.approx(n[8] * dan * cos_l ** 3, rel=1e-12, abs=1e-18)
    assert got[9] == pytest.approx(n[9] * cos_l * cos_l ** 3, rel=1e-12, abs=1e-18)


def test_lift_slope_cos():
    """Simple sweep: dCL/da at small angles, incompressible, is the straight slope x cos L."""
    for cos_l in (0.9, 0.7, 0.5):
        s = _station(cos_l, mach=0.01)
        a0 = s.aerodata[aero.AERODATA_ROW["a0"] - 1, int(s.iaero[5]) if s.iaero[5] >= 0 else 0]
        swept = _speedups.getclcdcm(s, 5, 0.0, 1.0, 0.0)[1]
        s.sweep = None
        straight = _speedups.getclcdcm(s, 5, 0.0, 1.0, 0.0)[1]
        assert swept / straight == pytest.approx(cos_l, rel=2e-3), cos_l


def test_converged_swept_solution_is_consistent():
    r = rotol_like()
    xi = np.linspace(0.2, 1.0, 9)
    r.set_sweep(Sweep.from_angle(xi, np.clip(60.0 * (xi - 0.45), -10.0, 40.0)))
    assert r.solve_rpm(1800.0)
    s = r.s
    for i in range(s.ii):
        v = cscalc(s, i)
        alfa = float(s.beta[i]) - v.phi
        cl = _speedups.getclcdcm(s, i, alfa, v.w, float(s.re[i]))[0]
        assert cl == pytest.approx(float(s.cl[i]), rel=1e-9, abs=1e-12), i
    lam = r.effective_sweep_deg()
    # this blade twists 40 deg: the offsets along the chord turn with it, so the effective
    # tip sweep (25 deg) is well below the 40 deg planform angle
    assert 20.0 < lam[-1] < 40.0 and lam[0] < 12.0


def test_sweep_needs_more_blade_angle():
    """At the same power a swept blade needs more blade angle (its lift slope is lower)."""
    out = []
    for lam in (0.0, 35.0):
        r = rotol_like()
        if lam:
            xi = np.linspace(0.2, 1.0, 9)
            r.set_sweep(Sweep.from_angle(xi, np.where(xi > 0.5, (xi - 0.5) / 0.5 * lam, 0.0)))
        assert r.solve_power(600e3, fix="rpm", rpm=1800.0)
        out.append(float(np.interp(0.75, r.s.xi[: r.s.ii], r.s.beta[: r.s.ii])))
    assert out[1] > out[0] + math.radians(0.3)


def test_design_refuses_sweep():
    r = rotol_like()
    r.set_sweep(Sweep([0.2, 1.0], [0.0, 0.1]))
    with pytest.raises(NotImplementedError):
        r.design(nblds=3, tip_radius=1.7, hub_radius=0.34, vel=120.0, rpm=1800.0,
                 power=600e3, cl=0.5)
