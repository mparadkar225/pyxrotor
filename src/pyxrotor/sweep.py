"""Swept (scimitar) blades: the geometry of a curved quarter-chord line, and the effective
sweep that the compiled core's section model uses (not part of XROTOR).

XROTOR's blade is a straight lifting line along the radius, in the disk plane. Here each
station's quarter-chord point is displaced, in its own section plane, by

    s(r)  along the chord, toward the trailing edge (positive = swept back), and
    n(r)  normal to the chord, toward the pressure side's normal n = (-sin b, cos b)
          in (tangential, axial) coordinates (positive = raked downstream at b = 0),

both in units of the tip radius R, as functions of r/R. The offsets are fixed to the blade
sections, so a pitch change rotates them with the blade (about the radial pitch-change axis,
through s = n = 0). CH and BETA keep XROTOR's meaning: the chord and blade angle of the
cylindrical (streamwise) section. The AeroSection parameters describe the airfoil normal to
the quarter-chord line, as simple sweep theory treats it.

Effective sweep (CSCALC, every iteration): the angle between the line and the plane normal
to the local relative velocity W = (0, CI, SI) in (radial, tangential, axial) components,

    T = (1, y' - y/xi, z'),  y = s cos b - n sin b,  z = s sin b + n cos b,
    sin L = (T . W) / (|T| |W|).

y is the tangential offset (an arc at radius xi, hence the y/xi term of d(r theta)/dr) and z
the axial one; primes are d/d(r/R), with b' from the blade's twist. GETCLCDCM then evaluates
each section at the normal angle of attack tan(an) = tan(a) / cos L and the normal Mach
number M cos L, and returns CL = CLn cos^2 L and CD = CDf + (CDn - CDf) cos^3 L, where
CDf = CDMIN (Re/REREF)^REXP is the friction part (along W, unchanged per unit radius: the
longer swept strip and the slower normal flow cancel) and the rest is pressure drag in the
normal plane (cos^2 from the normal dynamic pressure and normal chord per unit radius, cos
from its component along W). The derivatives with respect to the angle and W follow; the
dependence of L on the velocities is left out of the Newton Jacobian (the residual is
exact, so the solution is too).

Induction (Rotor.set_sweep(..., induction=True), the default; `vrtxco` below): the vortex
formulation with XROTOR's discrete helical wake leaving from the swept and raked quarter-chord
points, and the curved bound vortex's own induction, smoothed within a chord-sized core
(Sweep.core) so that it converges as the stations are refined.

Not modeled: the other blades' bound vortices (XROTOR leaves them out too), lifting-surface
(chordwise) effects, the loss of sweep effectiveness near the root and tip beyond what the
curved line's induction gives, rake's dihedral-like effect on the section, and blade
deflection.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import math

import numpy as np

from .numerics import deval, seval, spline


@dataclass
class Sweep:
    """Quarter-chord offsets (units of R) at r/R = xi: s along the chord (aft positive), n
    normal to it. Splined (XROTOR's SPLINE) for values and slopes at the stations."""
    xi: np.ndarray
    s: np.ndarray
    n: np.ndarray = field(default=None)
    core: float = 0.25       # bound-vortex core radius / local chord (vrtxco)

    def __post_init__(self):
        self.xi = np.asarray(self.xi, dtype=float)
        self.s = np.asarray(self.s, dtype=float)
        self.n = np.zeros_like(self.s) if self.n is None else np.asarray(self.n, dtype=float)
        if not (self.xi.shape == self.s.shape == self.n.shape) or self.xi.size < 2:
            raise ValueError("xi, s and n need the same length, at least 2")
        if np.any(np.diff(self.xi) <= 0):
            raise ValueError("xi must increase")
        self._ss = spline(self.s, self.xi, self.xi.size)
        self._ns = spline(self.n, self.xi, self.xi.size)

    @classmethod
    def from_angle(cls, xi, sweep_deg, zero_at=None, n=None):
        """A blade swept in its section planes by the angle sweep_deg(xi) (the planform
        angle of the quarter-chord line from the radial direction, aft positive): s' =
        tan(sweep). zero_at: r/R where the line crosses the pitch-change axis (s = 0); by
        default where the sweep changes sign, or the first station if it doesn't."""
        xi = np.asarray(xi, dtype=float)
        t = np.tan(np.radians(np.asarray(sweep_deg, dtype=float)))
        s = np.concatenate([[0.0], np.cumsum(0.5 * (t[1:] + t[:-1]) * np.diff(xi))])
        if zero_at is None:
            k = np.nonzero(np.diff(np.sign(t)))[0]
            zero_at = xi[0] if k.size == 0 else float(
                xi[k[0]] - t[k[0]] * (xi[k[0] + 1] - xi[k[0]]) / (t[k[0] + 1] - t[k[0]]))
        s = s - float(np.interp(zero_at, xi, s))
        return cls(xi, s, n)

    def at(self, x):
        """(s, s', n, n') at r/R = x."""
        k = self.xi.size
        return (seval(x, self.s, self._ss, self.xi, k), deval(x, self.s, self._ss, self.xi, k),
                seval(x, self.n, self._ns, self.xi, k), deval(x, self.n, self._ns, self.xi, k))


def station_arrays(s) -> np.ndarray:
    """(5, IX) C-contiguous: s, s', n, n', BETA' at the State's stations (for the core)."""
    ix = s.xi.size
    out = np.zeros((5, ix))
    ii = s.ii
    sw = s.sweep
    if ii >= 2:
        xi = np.array(s.xi[:ii], dtype=float)
        beta = np.array(s.beta[:ii], dtype=float)
        bs = spline(beta, xi, ii)
        for i in range(ii):
            out[:4, i] = sw.at(float(xi[i]))
            out[4, i] = deval(float(xi[i]), beta, bs, xi, ii)
    return out


def effective_sweep(s) -> np.ndarray:
    """Effective sweep angle (rad) at the stations, from the last CSCALC (arccos of the
    State's sweep_cos)."""
    return np.arccos(np.clip(np.array(s.sweep_cos[: s.ii]), -1.0, 1.0))


# ------------------------------------------------------------- induction of a swept blade
def blade_points(s, r):
    """Quarter-chord points of blade 1 at radii r (array, units of R) in VRTXCO's frame:
    x axial (downstream, with XROTOR's RAKE), blade 1 along +y, theta toward +z against the
    rotation (the way the wake helices turn). Returns (x, theta)."""
    from .numerics import seval, spline
    ii = s.ii
    xi = np.array(s.xi[:ii], dtype=float)
    be = np.array(s.beta[:ii], dtype=float)
    bs = spline(be, xi, ii)
    tanrak = np.tan(float(s.rake))
    xs, ths = [], []
    for rv in np.atleast_1d(r):
        rv = float(rv)
        b = seval(rv, be, bs, xi, ii)
        sv, _, nv, _ = s.sweep.at(rv)
        y = sv * np.cos(b) - nv * np.sin(b)
        z = sv * np.sin(b) + nv * np.cos(b)
        xs.append(rv * tanrak + z)
        ths.append(y / rv)
    return np.array(xs), np.array(ths)


def vrtxco(s) -> None:
    """VRTXCO for a swept blade (not XROTOR): XROTOR's discrete helical vortex wake, with each
    trailing leg leaving from its swept and raked quarter-chord point, plus the induction of
    the blade's own bound vortex, now curved (straight segments between the leg points
    XV(J), panel J carrying GAM(J)). Velocities are evaluated at each station's point on its
    panel's bound segment (the curve point at XI(I), projected onto that segment), in the
    local radial and tangential directions there.

    A curved lifting line's self-induction diverges logarithmically as the panels shrink;
    a real blade's bound vorticity is spread over its chord, so the bound vortex's induction
    is smoothed within a core of Sweep.core times the station's chord (default 0.25: the
    quarter-chord line to the 3/4-chord point is c/2, so a core of c/4 is the middle of
    the range a chordwise distribution spans). The results depend on it; see
    tests/model/test_sweep.py and validation/nasa_propfan.

    Kept from VRTXCO: its helix spacing (induced.theta_spacing) and pitch (ADW), the sign
    conventions, and the omission of the *other* blades' bound vortices. A straight blade
    (any rigid rotation or axial shift of one) gives XROTOR's coefficients within rounding
    (tests/model/test_sweep.py)."""
    from .induced import PI_EXACT, theta_spacing, vorsegvel
    ii, nblds, adw = s.ii, s.nblds, float(s.adw)
    if s.duct:
        raise NotImplementedError("swept-blade induction has no ducted model")
    s.vind_gam[:, :ii, :ii] = 0.0
    s.vind_adw[:, :ii] = 0.0
    thetspc = theta_spacing(adw)
    xitip = float(s.xv[ii])
    xv = np.array(s.xv[: ii + 1], dtype=float)
    x0, th0 = blade_points(s, xv)
    dtbld = 2.0 * PI_EXACT / float(nblds)

    # trailing legs: (3, NSEG) start/end points per leg J, all blades, in VRTXCO's order
    r1s, r2s, a1s, a2s = [], [], [], []
    for j in range(ii + 1):
        px, py, pz, pa = [], [], [], []
        for nb in range(nblds):
            ang = th0[j] + thetspc + nb * dtbld
            px.append(x0[j] + thetspc * xitip * adw)
            py.append(xv[j] * np.cos(ang))
            pz.append(xv[j] * np.sin(ang))
            pa.append(thetspc * xitip)
        r1s.append(np.array([np.concatenate([p[:-1] for p in q]) for q in (px, py, pz)]))
        r2s.append(np.array([np.concatenate([p[1:] for p in q]) for q in (px, py, pz)]))
        a1s.append(np.concatenate([p[:-1] for p in pa]))
        a2s.append(np.concatenate([p[1:] for p in pa]))
    nseg = r1s[0].shape[1]
    r1 = np.concatenate(r1s, axis=1)
    r2 = np.concatenate(r2s, axis=1)
    a1 = np.concatenate(a1s)
    a2 = np.concatenate(a2s)
    # bound vortex of blade 1: points at XV(J)
    pb = np.array([x0, xv * np.cos(th0), xv * np.sin(th0)])
    # control points
    xc, thc = blade_points(s, np.array(s.xi[:ii], dtype=float))
    q = np.array([xc, np.array(s.xi[:ii]) * np.cos(thc), np.array(s.xi[:ii]) * np.sin(thc)])

    gam = np.array(s.gam[:ii], dtype=float)
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        for i in range(ii):
            d = pb[:, i + 1] - pb[:, i]
            tpar = float(np.clip(np.dot(q[:, i] - pb[:, i], d) / np.dot(d, d), 0.0, 1.0))
            c = (pb[:, i] + tpar * d)[:, None]
            th = math.atan2(float(c[2, 0]), float(c[1, 0]))
            e_r = np.array([0.0, math.cos(th), math.sin(th)])
            e_t = np.array([0.0, -math.sin(th), math.cos(th)])

            def local(v):                       # (3, ...) -> XROTOR's (axial, radial, -tang)
                return np.array([v[0], np.tensordot(e_r, v, 1), -np.tensordot(e_t, v, 1)])
            uvw, ua1, ub1 = vorsegvel(r1 - c, r2 - c, full=False)
            leg = local(uvw.reshape(3, ii + 1, nseg).sum(axis=2))
            dadw = local((ua1 * a1 + ub1 * a2).reshape(3, ii + 1, nseg).sum(axis=2))
            ba, bb = pb[:, :-1] - c, pb[:, 1:] - c
            bu, _, _ = vorsegvel(ba, bb, full=False)
            bu[:, i] = 0.0                      # a point on its own segment: no induction
            # core: the bound vorticity is spread over the chord, so its induction is smoothed
            # within delta = core x chord of the segment's line, by h^2/(h^2 + delta^2)
            seg = bb - ba
            h2 = np.sum(np.cross(ba.T, bb.T) ** 2, axis=1) / np.sum(seg * seg, axis=0)
            delta = s.sweep.core * float(s.ch[i])
            bnd = local(np.nan_to_num(bu) * (h2 / (h2 + delta * delta)))
            for j in range(ii + 1):
                if j < ii:
                    s.vind_gam[:, i, j] += -leg[:, j]
                    s.vind_adw[:, i] -= dadw[:, j] * gam[j]
                if j > 0:
                    s.vind_gam[:, i, j - 1] += leg[:, j]
                    s.vind_adw[:, i] += dadw[:, j] * gam[j - 1]
            s.vind_gam[:, i, :ii] += bnd
