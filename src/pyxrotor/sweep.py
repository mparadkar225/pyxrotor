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

Not modeled (stage 3): the induction of a curved bound vortex and of a wake leaving from
swept and raked positions, the loss of sweep effectiveness near the root and tip, rake's
dihedral-like effect on the section, and blade deflection.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .numerics import deval, seval, spline


@dataclass
class Sweep:
    """Quarter-chord offsets (units of R) at r/R = xi: s along the chord (aft positive), n
    normal to it. Splined (XROTOR's SPLINE) for values and slopes at the stations."""
    xi: np.ndarray
    s: np.ndarray
    n: np.ndarray = field(default=None)

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
