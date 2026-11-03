"""The induced-velocity field evaluator (pyxrotor.field), a new model rather than a port.

Validation chain, each step independent of the code it checks:
  1. elliptic integrals against SciPy;
  2. the semi-infinite vortex cylinder's closed form against brute-force Biot-Savart (rings
     of straight segments, integrated downstream; no elliptic integrals);
  3. the rotor field's exact limits: XROTOR's own disk-averaged velocities on the disk,
     twice them far downstream, nothing far upstream, swirl only behind the disk, and zero
     divergence;
  4. the whole rotor field against an azimuth-averaged discrete helical vortex wake (bound
     and trailing filaments, Biot-Savart), whose sign conventions are fixed by physics.
"""
import numpy as np
import pytest
import scipy.special as sp

from pyxrotor.field import InducedField, cylinder_velocity, ellipe, ellipk, ellippi
from pyxrotor.slipstream import profile
from rotors import solved_state

pytestmark = pytest.mark.model


# ------------------------------------------------------------------------- 1. elliptic
def test_elliptic_integrals_match_scipy():
    rng = np.random.default_rng(1)
    m = np.concatenate([np.linspace(0, 0.999999, 3000), 1 - np.logspace(-14, -1, 300)])
    np.testing.assert_allclose(ellipk(m), sp.ellipk(m), rtol=5e-15)
    np.testing.assert_allclose(ellipe(m), sp.ellipe(m), rtol=5e-14)
    n = np.concatenate([rng.uniform(-5, 0.999999, 5000), 1 - np.logspace(-12, -1, 50)])
    mm = np.concatenate([rng.uniform(0, 0.99999, 5000), np.full(50, 0.5)])
    ref = sp.elliprf(0, 1 - mm, 1) + n / 3 * sp.elliprj(0, 1 - mm, 1, 1 - n)
    np.testing.assert_allclose(ellippi(n, mm), ref, rtol=1e-12)


# ------------------------------------------------------------------- Biot-Savart tools
def segment_velocity(a_pts, b_pts, gam, p):
    """Velocity at points p (3, N) of straight vortex segments a->b (3, M) with
    circulation gam (M): exact Biot-Savart, summed over segments. Returns (3, N)."""
    out = np.zeros((3, p.shape[1]))
    for k in range(p.shape[1]):
        a = a_pts - p[:, k:k + 1]
        b = b_pts - p[:, k:k + 1]
        am, bm = np.linalg.norm(a, axis=0), np.linalg.norm(b, axis=0)
        axb = np.cross(a.T, b.T).T
        den = am * bm * (am * bm + np.einsum("ij,ij->j", a, b))
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(den > 1e-30, gam * (am + bm) / den, 0.0) / (4.0 * np.pi)
        out[:, k] = (axb * f).sum(axis=1)
    return out


def ring(a, z, n):
    t = np.linspace(0.0, 2.0 * np.pi, n + 1)
    pts = np.array([a * np.cos(t), a * np.sin(t), np.full_like(t, z)])
    return pts[:, :-1], pts[:, 1:]


def brute_cylinder(gamma, a, r, z, nring=720, nz=400):
    """Semi-infinite cylinder as rings at z' in [0, inf) (z' = a tan(pi s / 2), Gauss-
    Legendre in s). Ring orientation: +theta, which accelerates the flow inside along +z."""
    s, w = np.polynomial.legendre.leggauss(nz)
    s = 0.5 * (s + 1.0)
    w = 0.5 * w
    zp = a * np.tan(0.5 * np.pi * s)
    dz = a * 0.5 * np.pi / np.cos(0.5 * np.pi * s) ** 2
    p = np.array([[r], [0.0], [z]])
    u = np.zeros(3)
    for zk, wk, dk in zip(zp, w, dz):
        A, B = ring(a, zk, nring)
        u += wk * dk * segment_velocity(A, B, np.full(nring, gamma), p)[:, 0]
    return u[2], u[0]                       # axial, radial (field point at theta = 0)


# ------------------------------------------------------------------ 2. vortex cylinder
@pytest.mark.parametrize("r,z", [(0.3, 0.0), (0.3, 0.4), (0.3, -0.4), (0.8, 1.5), (0.8, -1.5),
                                 (1.3, 0.2), (1.3, -0.2), (2.0, 3.0), (0.95, 0.3), (1.05, -0.3),
                                 (0.5, -5.0), (0.5, 8.0)])
def test_cylinder_closed_form_vs_biot_savart(r, z):
    uz, ur = cylinder_velocity(1.0, 1.0, r, z)
    bz, br = brute_cylinder(1.0, 1.0, r, z)
    assert abs(uz - bz) < 2e-5 and abs(ur - br) < 2e-5, (uz, bz, ur, br)


def test_cylinder_limits():
    r = np.array([0.0, 0.2, 0.5, 0.9, 1.1, 2.0])
    uz0, ur0 = cylinder_velocity(1.0, 1.0, r, 0.0)
    np.testing.assert_allclose(uz0, np.where(r < 1, 0.5, 0.0), atol=1e-14)       # disk: half
    uzd, urd = cylinder_velocity(1.0, 1.0, r, 1e6)
    np.testing.assert_allclose(uzd, np.where(r < 1, 1.0, 0.0), atol=1e-9)        # far wake
    uzu, uru = cylinder_velocity(1.0, 1.0, r, -1e6)
    np.testing.assert_allclose(uzu, 0.0, atol=1e-9)                              # far upstream
    # radial velocity even in z, inward for a contracting (accelerating) cylinder
    z = np.array([0.1, 0.5, 2.0])
    _, a = cylinder_velocity(1.0, 1.0, 0.7, z)
    _, b = cylinder_velocity(1.0, 1.0, 0.7, -z)
    np.testing.assert_allclose(a, b, rtol=1e-12)
    assert (a < 0).all()
    # on the axis: the solenoid formula, continuous with points just off it
    zz = np.linspace(-3, 3, 13)
    np.testing.assert_allclose(cylinder_velocity(1.0, 1.0, 0.0, zz)[0],
                               cylinder_velocity(1.0, 1.0, 1e-7, zz)[0], atol=1e-9)


def test_cylinder_is_divergence_free():
    h = 1e-5
    for r, z in [(0.4, 0.3), (0.6, -0.7), (1.4, 0.5), (2.5, -1.0)]:
        dudz = (cylinder_velocity(1.0, 1.0, r, z + h)[0] - cylinder_velocity(1.0, 1.0, r, z - h)[0]) / (2 * h)
        rur = lambda rr: rr * cylinder_velocity(1.0, 1.0, rr, z)[1]   # noqa: E731
        div = dudz + (rur(r + h) - rur(r - h)) / (2 * h) / r
        assert abs(div) < 1e-6, (r, z, div)


# ------------------------------------------------------------------------ 3. rotor limits
@pytest.fixture(scope="module")
def rotor():
    s = solved_state()
    assert s.conv
    return s, InducedField.from_state(s)


def test_disk_values_are_xrotors_averages(rotor):
    """On the disk: half of what VPUT writes (its far-wake 2*Va, 2*Vt), station by station."""
    s, f = rotor
    r, u, v = f.disk_average_profile()
    rp, up, vp = profile(s)
    np.testing.assert_allclose(r, rp, rtol=1e-14)
    np.testing.assert_allclose(u, 0.5 * up, rtol=1e-8)         # XROTOR's PI = 3.141592654
    np.testing.assert_allclose(v, 0.5 * vp, rtol=1e-8)


def test_far_field_limits(rotor):
    s, f = rotor
    r, u0, v0 = f.disk_average_profile()
    u, w, v = f.velocity(r, 1e5 * s.rad)
    np.testing.assert_allclose(u, 2 * u0, rtol=1e-6)                # far wake: twice
    np.testing.assert_allclose(v, 2 * v0, rtol=1e-14)
    assert np.abs(w).max() < 1e-6 * np.abs(u).max()
    u, w, v = f.velocity(r, -1e5 * s.rad)
    assert np.abs(u).max() < 1e-6 * np.abs(u0).max()                 # far upstream: nothing
    assert (v == 0).all()                                            # no swirl ahead, anywhere
    _, _, v = f.velocity(r, -1e-9)
    assert (v == 0).all()
    u, _, v = f.velocity(1.5 * s.rad, np.array([-0.3, 0.3, 3.0]) * s.rad)
    assert (v == 0).all()                                            # outside the wake


def test_upstream_decay_matches_actuator_disk_on_axis_shape(rotor):
    """Upstream the axial velocity falls monotonically with distance, faster near the loaded
    outer span (close to the trailing sheets) than the on-axis actuator-disk formula
    1 - d/sqrt(d^2 + R^2) suggests: at 0.2 R about 0.6 at r = 0.7 R, against 0.80 on axis."""
    s, f = rotor
    zs = np.array([0.0, -0.05, -0.1, -0.2, -0.3, -0.5, -1.0, -2.0]) * s.rad
    for rr in (0.3, 0.5, 0.7, 0.9):
        u = f.velocity(np.full_like(zs, rr * s.rad), zs)[0]
        assert np.all(np.diff(u) < 0)
        assert 0.4 < u[3] / u[0] < 0.9


def test_rotor_field_divergence_free(rotor):
    s, f = rotor
    h = 1e-6 * s.rad
    for rr, zz in [(0.5, 0.2), (0.75, -0.3), (1.2, 0.4), (0.33, 1.0)]:
        r, z = rr * s.rad, zz * s.rad
        dudz = (f.velocity(r, z + h)[0] - f.velocity(r, z - h)[0]) / (2 * h)
        dr = (f.velocity(r + h, z)[1] * (r + h) - f.velocity(r - h, z)[1] * (r - h)) / (2 * h) / r
        assert abs(dudz + dr) < 1e-5 * s.vel / s.rad, (rr, zz, dudz + dr)


# ------------------------------------------------------- 4. discrete helical vortex wake
def discrete_wake(s, length=40.0, seg_per_turn=72):
    """Bound and trailing vortex filaments of XROTOR's rigid helical wake, in m and m^2/s.
    Blades rotate toward +theta about +z (downstream); the wake of blade b at distance z
    left the blade when it was z/h radians behind."""
    B, R, V = int(s.nblds), float(s.rad), float(s.vel)
    ii = s.ii
    xv = np.array(s.xv[: ii + 1]) * R
    gam = np.array(s.gam[:ii]) * V * R
    h = float(s.adw) * R
    shed = np.diff(np.concatenate([[0.0], gam, [0.0]]))          # at xv[j], downstream
    nphi = int(length * R / (2 * np.pi * h) * seg_per_turn) + 1
    phi = np.linspace(0.0, length * R / h, nphi)
    A, Bp, G = [], [], []
    for b in range(B):
        tb = 2 * np.pi * b / B
        # bound vortex: from tip toward hub on the blade (thrust with positive GAM)
        for i in range(ii):
            A.append([[xv[i + 1] * np.cos(tb)], [xv[i + 1] * np.sin(tb)], [0.0]])
            Bp.append([[xv[i] * np.cos(tb)], [xv[i] * np.sin(tb)], [0.0]])
            G.append([gam[i]])
        for j, a in enumerate(xv):
            pts = np.array([a * np.cos(tb - phi), a * np.sin(tb - phi), h * phi])
            A.append(pts[:, :-1])
            Bp.append(pts[:, 1:])
            G.append(np.full(nphi - 1, shed[j]))
    return (np.concatenate(A, axis=1), np.concatenate(Bp, axis=1), np.concatenate(G))


def averaged(s, segs, r, z, naz=24):
    """Azimuth average of the discrete field at (r, z): (axial, radial, swirl in +theta)."""
    A, Bp, G = segs
    B = int(s.nblds)
    psi = (np.arange(naz) + 0.5) * (2 * np.pi / B) / naz
    p = np.array([r * np.cos(psi), r * np.sin(psi), np.full(naz, z)])
    u = segment_velocity(A, Bp, G, p)
    ur = u[0] * np.cos(psi) + u[1] * np.sin(psi)
    ut = -u[0] * np.sin(psi) + u[1] * np.cos(psi)
    return u[2].mean(), ur.mean(), ut.mean()


@pytest.fixture(scope="module")
def discrete(rotor):
    s, _ = rotor
    return discrete_wake(s)


def test_discrete_wake_physics(rotor, discrete):
    """Independent of the cylinder model: a thrusting rotor accelerates the flow downstream
    and drags it round in its direction of rotation behind the disk, not ahead of it."""
    s, _ = rotor
    u, w, v = averaged(s, discrete, 0.7 * s.rad, 2.0 * s.rad)
    assert u > 0 and v > 0 and w < 0
    u, w, v = averaged(s, discrete, 0.7 * s.rad, -0.5 * s.rad)
    assert u > 0 and abs(v) < 0.02 * u


@pytest.mark.parametrize("rr,zz", [(0.7, 0.3), (0.7, -0.3), (0.5, 1.0), (0.9, -0.15),
                                   (0.4, -1.0), (1.3, 0.5), (1.3, -0.5), (0.85, 2.0)])
def test_field_matches_discrete_vortex_wake(rotor, discrete, rr, zz):
    """The averaged closed-form field against the averaged discrete wake, off the disk.
    Tolerance: the wake is truncated at 40 R (about 3e-4 of the far-wake velocity) and the
    helices are polygons."""
    s, f = rotor
    r, z = rr * s.rad, zz * s.rad
    ref = averaged(s, discrete, r, z, naz=192)     # filaments pass close to the samples
    got = [float(x) for x in f.velocity(r, z)]
    scale = float(np.abs(f.velocity(np.array(s.xi[: s.ii]) * s.rad, 1e4 * s.rad)[0]).max())
    for name, g, d in zip(("axial", "radial", "swirl"), got, ref):
        assert abs(g - d) < 3e-3 * scale, (name, rr, zz, g, d, scale)
