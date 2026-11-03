"""Induced-velocity field of a solved rotor anywhere in space (not part of XROTOR).

For rotor-rotor interaction (contra-rotating and tandem rotors) each rotor must know the
velocity the other induces at its blades, upstream or downstream of the other's disk.
XROTOR only evaluates induced velocities on its own disk. This module extends its model
off the disk, circumferentially averaged (the steady part a rotor in another's flow sees).

Model: XROTOR's rigid helical wake. B blades carry circulation GAM(i) on station i (between
vortex radii XV(i) and XV(i+1)); trailing vortices of strength GAM(j) - GAM(j-1) leave each
radius XV(j) on helices of pitch ADW*R per radian, with no contraction. Averaged over
azimuth, this vorticity is exactly
  * semi-infinite axial vortex sheets plus the radial bound-vortex sheet on the disk, whose
    swirl is, by Stokes' theorem, exactly v = B*Gamma(r)/(2 pi r) behind the disk, half of it
    on the disk and zero ahead of it (Gamma piecewise constant per station, 0 outside the
    blade); and no axial or radial velocity;
  * semi-infinite vortex cylinders of ring vorticity gamma_j = -B*dGamma_j/(2 pi h) at each
    XV(j) (h = ADW*R), whose axial and radial velocities have closed forms in complete
    elliptic integrals (Branlard & Gaunaa 2015, "Cylindrical vortex wake model: right
    cylinder", Wind Energy 18; gamma is the far-wake axial velocity inside the cylinder).

On the disk this reproduces the circumferentially averaged velocities XROTOR's VPUT uses
(Va = B*Gamma/(4 pi h), Vt = B*Gamma/(4 pi r)), and far downstream twice those. The radial
velocity, which XROTOR's model ignores, is returned as well.

Conventions: z along the flight direction's opposite (positive downstream, disk at z = 0),
r from the axis; axial velocity positive downstream (the direction the rotor accelerates
flow), radial positive outward, swirl positive in the rotor's direction of rotation. All
are perturbations (no free stream, no imposed or nacelle velocities), in m/s.

Not modeled: slipstream contraction (and so the radial displacement of the wake), blade-
passing unsteadiness, ducts. Field points exactly on a vortex sheet (r equal to an XV
radius, or on the disk with r = a) are singular; the jump there is resolved one-sided.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PI = np.pi


# ----------------------------------------------------------------------------- elliptic
def _carlson_rf(x, y, z):
    """Carlson's R_F(x, y, z) by duplication (x, y, z >= 0, at most one zero)."""
    x, y, z = (np.array(v, dtype=float, copy=True) for v in np.broadcast_arrays(x, y, z))
    for _ in range(60):
        lam = np.sqrt(x * y) + np.sqrt(y * z) + np.sqrt(z * x)
        x, y, z = 0.25 * (x + lam), 0.25 * (y + lam), 0.25 * (z + lam)
        mu = (x + y + z) / 3.0
        if np.all(np.maximum.reduce([abs(x - mu), abs(y - mu), abs(z - mu)]) <= 1e-9 * abs(mu)):
            break
    mu = (x + y + z) / 3.0
    dx, dy, dz = 1.0 - x / mu, 1.0 - y / mu, 1.0 - z / mu
    e2 = dx * dy - dz * dz
    e3 = dx * dy * dz
    return (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0) / np.sqrt(mu)


def _carlson_rj(x, y, z, p):
    """Carlson's R_J(x, y, z, p) (x, y, z >= 0, at most one zero, p > 0), by the duplication
    algorithm of B. C. Carlson, Numer. Algorithms 10 (1995), eqs. 2.36-2.41; R_D = R_J(x,y,z,z)."""
    x, y, z, p = (np.array(v, dtype=float, copy=True) for v in np.broadcast_arrays(x, y, z, p))
    a0 = (x + y + z + 2.0 * p) / 5.0
    a = a0.copy()
    delta = (p - x) * (p - y) * (p - z)
    x0, y0, z0 = x.copy(), y.copy(), z.copy()
    total = np.zeros_like(x)
    m = 0
    for m in range(40):
        sx, sy, sz, sp = np.sqrt(x), np.sqrt(y), np.sqrt(z), np.sqrt(p)
        lam = sx * sy + sy * sz + sz * sx
        d = (sp + sx) * (sp + sy) * (sp + sz)
        e = 4.0 ** (-3 * m) * delta / (d * d)
        total = total + 4.0 ** (-m) * _carlson_rc(1.0, 1.0 + e) / d
        x, y, z, p = 0.25 * (x + lam), 0.25 * (y + lam), 0.25 * (z + lam), 0.25 * (p + lam)
        a = 0.25 * (a + lam)
        spread = np.maximum.reduce([abs(a0 - x0), abs(a0 - y0), abs(a0 - z0), abs(a0 - (2.5 * a0 - 0.5 * (x0 + y0 + z0)))])
        if np.all(4.0 ** (-(m + 1)) * spread <= 1e-4 * abs(a)):
            break
    fm = 4.0 ** (m + 1)
    X = (a0 - x0) / (fm * a)
    Y = (a0 - y0) / (fm * a)
    Z = (a0 - z0) / (fm * a)
    P = -(X + Y + Z) / 2.0
    e2 = X * Y + X * Z + Y * Z - 3.0 * P * P
    e3 = X * Y * Z + 2.0 * e2 * P + 4.0 * P * P * P
    e4 = (2.0 * X * Y * Z + e2 * P + 3.0 * P * P * P) * P
    e5 = X * Y * Z * P * P
    series = (1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0 - 3.0 * e4 / 22.0
              - 9.0 * e2 * e3 / 52.0 + 3.0 * e5 / 26.0)
    return series / (fm * a * np.sqrt(a)) + 6.0 * total


def _carlson_rc(x, y):
    """Carlson's R_C(x, y) = R_F(x, y, y), y > 0, by duplication. (The closed forms with
    arccos/arccosh lose all precision as x -> y, which is where R_J's sum drives it.)"""
    x, y = (np.array(v, dtype=float, copy=True) for v in np.broadcast_arrays(x, y))
    for _ in range(60):
        lam = 2.0 * np.sqrt(x) * np.sqrt(y) + y
        x, y = 0.25 * (x + lam), 0.25 * (y + lam)
        mu = (x + 2.0 * y) / 3.0
        s = (y - x) / (3.0 * mu)
        if np.all(abs(s) <= 1e-4):
            break
    mu = (x + 2.0 * y) / 3.0
    s = (y - x) / (3.0 * mu)
    return (1.0 + s * s * (0.3 + s * (1.0 / 7.0 + s * (0.375 + s * 9.0 / 22.0)))) / np.sqrt(mu)


def ellipk(m):
    """K(m), complete elliptic integral of the first kind (parameter m = k^2 < 1)."""
    return _carlson_rf(0.0, 1.0 - np.asarray(m, float), 1.0)


def ellipe(m):
    """E(m), complete elliptic integral of the second kind."""
    m = np.asarray(m, float)
    return _carlson_rf(0.0, 1.0 - m, 1.0) - m / 3.0 * _carlson_rj(0.0, 1.0 - m, 1.0, 1.0)


def ellippi(n, m):
    """Pi(n, m), complete elliptic integral of the third kind (n < 1, m < 1)."""
    n, m = np.asarray(n, float), np.asarray(m, float)
    return _carlson_rf(0.0, 1.0 - m, 1.0) + n / 3.0 * _carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n)


# ---------------------------------------------------------------------- vortex cylinder
def cylinder_velocity(gamma, a, r, z):
    """Axial and radial velocity at (r, z) of a semi-infinite vortex cylinder of radius a,
    from z = 0 to +inf, with ring vorticity gamma per unit length (so the axial velocity far
    inside it is gamma). Broadcasts over all arguments. On r = a, the inside value is used."""
    gamma, a, r, z = np.broadcast_arrays(*(np.asarray(v, float) for v in (gamma, a, r, z)))
    r = np.where(r == a, a * (1.0 - 1e-12), r)          # on the sheet: inside limit
    rr = np.maximum(r, 1e-300)
    k2 = 4.0 * rr * a / ((a + rr) ** 2 + z * z)
    k02 = 4.0 * rr * a / (a + rr) ** 2
    k = np.sqrt(k2)
    step = np.where(r < a, 1.0, 0.0)                    # (R - r + |R - r|) / (2 |R - r|)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = z * k / (2.0 * PI * np.sqrt(rr * a)) * (
            ellipk(k2) + (a - rr) / (a + rr) * ellippi(k02, k2))
        uz = 0.5 * gamma * (step + term)
        ur = -gamma / (2.0 * PI) * np.sqrt(a / rr) * ((2.0 - k2) / k * ellipk(k2) - 2.0 / k * ellipe(k2))
    ur = np.where(r <= 0.0, 0.0, ur)                    # on the axis by symmetry
    uz = np.where(r <= 0.0, 0.5 * gamma * (1.0 + z / np.sqrt(a * a + z * z)), uz)
    return uz, ur


# --------------------------------------------------------------------------- rotor field
@dataclass(frozen=True)
class InducedField:
    """Circumferentially averaged induced-velocity field of a rotor (see module docstring).
    Build it with InducedField.from_state(state) from a solved XROTOR State."""
    nblds: int
    rad: float              # tip radius, m
    vel: float              # flight speed, m/s (velocity scale of GAM)
    pitch: float            # wake helix pitch per radian, m (ADW * RAD)
    xv: np.ndarray          # vortex radii XV(0..II), / RAD
    gam: np.ndarray         # circulation per station GAM(0..II-1), / (VEL * RAD)

    @classmethod
    def from_state(cls, s) -> "InducedField":
        if getattr(s, "duct", False):
            raise NotImplementedError("InducedField: ducted rotors are not modeled")
        ii = s.ii
        return cls(nblds=int(s.nblds), rad=float(s.rad), vel=float(s.vel),
                   pitch=float(s.adw) * float(s.rad), xv=np.array(s.xv[: ii + 1], dtype=float),
                   gam=np.array(s.gam[:ii], dtype=float))

    # circulation (m^2/s) and ring vorticity of the trailing cylinders (m/s)
    def _gamma_dim(self) -> np.ndarray:
        return self.gam * (self.vel * self.rad)

    def _cylinders(self) -> tuple[np.ndarray, np.ndarray]:
        g = np.concatenate([[0.0], self._gamma_dim(), [0.0]])
        dgam = np.diff(g)                               # trailing strength at XV(j)
        return self.xv * self.rad, -self.nblds * dgam / (2.0 * PI * self.pitch)

    def circulation(self, r) -> np.ndarray:
        """Total bound circulation B*Gamma(r) (m^2/s), piecewise constant per station."""
        r = np.asarray(r, float) / self.rad
        j = np.searchsorted(self.xv, r, side="right") - 1
        inside = (j >= 0) & (j < len(self.gam))
        g = np.where(inside, self.gam[np.clip(j, 0, len(self.gam) - 1)], 0.0)
        return self.nblds * g * (self.vel * self.rad)

    def velocity(self, r, z):
        """(axial, radial, swirl) induced velocity in m/s at radius r (m) and axial position z
        (m, positive downstream of the disk). Broadcasts over r and z."""
        r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
        a, gam = self._cylinders()
        uz, ur = cylinder_velocity(gam[:, None], a[:, None], r.ravel()[None, :], z.ravel()[None, :])
        u = np.add.reduce(uz, axis=0).reshape(r.shape)
        w = np.add.reduce(ur, axis=0).reshape(r.shape)
        step = np.where(z > 0.0, 1.0, np.where(z < 0.0, 0.0, 0.5))
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(r > 0.0, step * self.circulation(r) / (2.0 * PI * r), 0.0)
        return u, w, v

    def disk_average_profile(self):
        """At the station radii: the on-disk (axial, swirl) averages and the radii (m), the
        quantities XROTOR's VPUT uses (it writes twice these, the far-wake values)."""
        xi = 0.5 * (self.xv[:-1] + self.xv[1:])
        r = xi * self.rad
        u, _, v = self.velocity(r, 0.0)
        return r, u, v
