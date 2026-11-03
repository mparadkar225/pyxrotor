"""Contra-rotating (and co-rotating tandem) rotor pairs (not part of XROTOR).

Each rotor is a full XROTOR solve (a Rotor). It sees the other's circumferentially
averaged induced velocity (field.InducedField) at its own blade stations, imposed as
XROTOR's external inflow (UADD, VADD, used by UVADD). The front rotor sees the rear's
field from upstream, which is axial only, since there is no average swirl ahead of a disk.
The rear rotor sees the front's slipstream (axial and swirl) at the spacing, with the
swirl's sign flipped when it turns the other way. XROTOR's free-wake condition uses the
imposed inflow too (through its "equivalent prop"), so each rotor's own wake pitch responds
to the other rotor.

The two solves are iterated, with under-relaxation of the imposed inflow, until the
imposed velocities stop changing. Operating conditions (RotorPair.solve_*):
    solve_rpm(rpm)                         both blade angles fixed, rpm given
    solve_power(power, rpm, torque_split)  rpm and shaft power given, pitch of each rotor
                                           set to give the torque split (0.5: equal torque)
    solve_thrust(thrust, rpm, torque_split) the same with total thrust given
    solve_power_fixed_pitch(power, rpm)    both blade angles fixed, rpm found for the power
The rear runs at gear_ratio times the front's rpm.

Coordinates: z along the axis, positive downstream, front disk at z = 0, rear at the
spacing; swirl positive in the *front* rotor's direction of rotation.

Model limits (see field.py): no slipstream contraction, each wake a rigid helix of its
own pitch, no radial velocity at the blades, no viscous wakes or blade-passing unsteadiness.
Both rotors must have the same flight speed and fluid properties.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field as dc_field

import numpy as np

from .constants import PI
from .field import InducedField
from .slipstream import set_profile


class PairConvergenceError(RuntimeError):
    """The coupled iteration (or a rotor's own solve inside it) did not converge."""


@dataclass
class PairResult:
    converged: bool
    iterations: int                # coupled iterations in total (all outer evaluations)
    front: dict                    # Rotor.performance() of each rotor
    rear: dict
    thrust: float                  # N, total
    power: float                   # W, total shaft power
    efficiency: float              # T V / P
    torque_split: float            # front share of the total torque magnitude
    net_torque: float              # N m on the airframe axis (front direction positive)
    swirl_power_fraction: float    # far-wake swirl kinetic-energy flux / shaft power
    residual_swirl: dict = dc_field(repr=False)   # far wake: r (m), v (m/s, front direction)
    #: per coupled iteration, the largest change of the imposed velocities / VEL; for the
    #: thrust and fixed-pitch modes, the relative target error per outer evaluation
    history: list = dc_field(default_factory=list, repr=False)


class RotorPair:
    """Two coaxial rotors: front at z = 0, rear at z = spacing (m) behind it.

    front, rear: Rotor objects already set up (sections, blade, atmosphere, nacelle UBODY if
    any). counter_rotating: the rear turns the other way (else a co-rotating tandem).
    The pair writes each rotor's imposed slipstream (NADD, RADD, UADD, VADD), replacing any
    profile set before.
    """

    def __init__(self, front, rear, spacing: float, counter_rotating: bool = True,
                 gear_ratio: float = 1.0, relaxation: float = 0.7, tol: float = 1e-6,
                 max_iter: int = 60):
        self.front, self.rear = front, rear
        self.spacing = float(spacing)
        if self.spacing < 0.0:
            raise ValueError("spacing must be >= 0 (rear behind front)")
        self.sign = -1.0 if counter_rotating else 1.0
        self.gear_ratio = float(gear_ratio)
        self.relaxation = float(relaxation)
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self._check_common()

    # ------------------------------------------------------------------- set-up
    def _check_common(self):
        a, b = self.front.s, self.rear.s
        for name in ("vel", "rho", "vso", "rmu"):
            va, vb = getattr(a, name), getattr(b, name)
            if not math.isclose(va, vb, rel_tol=1e-12, abs_tol=0.0):
                raise ValueError(f"front and rear differ in {name.upper()} ({va} vs {vb})")

    def velocity(self, vel: float) -> None:
        """OPER VELO on both rotors."""
        self.front.velocity(vel)
        self.rear.velocity(vel)

    def atmosphere(self, alt_km: float) -> None:
        self.front.atmosphere(alt_km)
        self.rear.atmosphere(alt_km)

    def nacelle(self, z, area, wake_area: float, front_z: float = 0.0,
                fix_source_positions: bool = False) -> None:
        """A shared nacelle/spinner body (NACE) for both rotors: body stations z (m, same
        axis as the front rotor's position front_z) and areas (m^2; negative = radius).
        Each rotor gets the body's UBODY at its own position."""
        for rotor, zr in ((self.front, front_z), (self.rear, front_z + self.spacing)):
            rotor.nacelle(z=z, area=area, rotor_z=zr, wake_area=wake_area,
                          fix_source_positions=fix_source_positions)

    # --------------------------------------------------------------- interaction
    @staticmethod
    def _stations(rotor) -> np.ndarray:
        s = rotor.s
        return np.array(s.xi[: s.ii], dtype=float) * s.rad

    def _imposed(self, which: str, other_field: InducedField | None):
        """(r, u, v) to impose on `which` ('front' or 'rear'), at its station radii (m/s),
        from the other rotor's field; zeros if the other has no solution yet."""
        rotor = self.front if which == "front" else self.rear
        r = self._stations(rotor)
        if other_field is None:
            return r, np.zeros_like(r), np.zeros_like(r)
        z = -self.spacing if which == "front" else self.spacing
        u, _, v = other_field.velocity(r, np.full_like(r, z))
        # swirl in this rotor's own rotation direction (UVADD's WT convention)
        return r, u, self.sign * v

    @staticmethod
    def _impose(rotor, r, u, v):
        set_profile(rotor.s, r, u, v, 1.0, 1.0)
        rotor.s.conv = False

    # ------------------------------------------------------------------ iteration
    def _iterate(self, solve_front, solve_rear) -> tuple[bool, int, list]:
        """Fixed-point iteration: front, then rear, each solved with the other's latest
        field, the imposed inflow under-relaxed. Converged when the largest change of the
        imposed velocities (before relaxation), over both rotors, is below tol * VEL.
        Starts from the rotors' current solutions (the other's field is zero if a rotor has
        none). Returns (converged, iterations, history of that change)."""
        w = self.relaxation
        vel = self.front.s.vel
        imposed = {}
        for which, rotor in (("front", self.front), ("rear", self.rear)):
            r = self._stations(rotor)
            imposed[which] = (r, np.zeros_like(r), np.zeros_like(r))
        fields = {"front": self._field(self.front), "rear": self._field(self.rear)}
        hist = []
        for it in range(1, self.max_iter + 1):
            change = 0.0
            for which, solve in (("front", solve_front), ("rear", solve_rear)):
                rotor = self.front if which == "front" else self.rear
                other = fields["rear" if which == "front" else "front"]
                r, u, v = self._imposed(which, other)
                _, u0, v0 = imposed[which]
                change = max(change, float(np.max(np.abs(u - u0) + np.abs(v - v0))) / vel)
                if it > 1:                                 # first pass: take it whole
                    u, v = u0 + w * (u - u0), v0 + w * (v - v0)
                imposed[which] = (r, u, v)
                self._impose(rotor, r, u, v)
                if not solve(rotor):
                    raise PairConvergenceError(f"{which} rotor did not converge (iteration {it})")
                fields[which] = self._field(rotor)
            hist.append(change)
            if it > 1 and change < self.tol:
                return True, it, hist
        return False, self.max_iter, hist

    @staticmethod
    def _field(rotor) -> InducedField | None:
        return InducedField.from_state(rotor.s) if rotor.s.conv else None

    def _rpm_rear(self, rpm: float) -> float:
        return self.gear_ratio * float(rpm)

    # ----------------------------------------------------------- operating modes
    def solve_rpm(self, rpm: float) -> PairResult:
        """Both blade angles fixed; the front at rpm, the rear at gear_ratio * rpm."""
        rf, rr = float(rpm), self._rpm_rear(rpm)
        ok, n, h = self._iterate(lambda r: r.solve_rpm(rf), lambda r: r.solve_rpm(rr))
        return self._result(ok, n, h)

    def _torques(self, total_q: float, split: float):
        return split * total_q, (1.0 - split) * total_q

    def _solve_torques(self, qf: float, qr: float, rpm: float):
        rf, rr = float(rpm), self._rpm_rear(rpm)
        return self._iterate(lambda r: r.solve_torque(qf, fix="rpm", rpm=rf),
                             lambda r: r.solve_torque(qr, fix="rpm", rpm=rr))

    def solve_power(self, power: float, rpm: float, torque_split: float = 0.5) -> PairResult:
        """Shaft power (W) at rpm: each rotor's blade angle is set (XROTOR TORQ, fixed rpm)
        to absorb its share of the torque, front share = torque_split."""
        if not 0.0 < torque_split < 1.0:
            raise ValueError("torque_split must be between 0 and 1")
        wf = float(rpm) * PI / 30.0
        wr = self._rpm_rear(rpm) * PI / 30.0
        total_q = float(power) / (torque_split * wf + (1.0 - torque_split) * wr)
        ok, n, h = self._solve_torques(*self._torques(total_q, torque_split), rpm)
        return self._result(ok, n, h)

    def solve_thrust(self, thrust: float, rpm: float, torque_split: float = 0.5,
                     tol: float = 1e-6, max_iter: int = 30) -> PairResult:
        """Total thrust (N) at rpm with the given torque split: secant iteration on the
        total torque around the coupled solve."""
        thrust = float(thrust)
        count, hist = [0], []

        def t_of(q):
            ok, n, h = self._solve_torques(*self._torques(q, torque_split), rpm)
            count[0] += n
            if not ok:
                raise PairConvergenceError("coupled iteration did not converge")
            err = self._totals()[0] - thrust
            hist.append(abs(err) / abs(thrust))
            return err

        q0 = self._guess_torque(thrust, rpm)
        q1 = 1.1 * q0
        f0, f1 = t_of(q0), t_of(q1)
        for _ in range(max_iter):
            if abs(f1) <= tol * abs(thrust):
                return self._result(True, count[0], hist)
            q0, q1 = q1, q1 - f1 * (q1 - q0) / (f1 - f0)
            f0, f1 = f1, t_of(q1)
        return self._result(False, count[0], hist)

    def solve_power_fixed_pitch(self, power: float, rpm_guess: float, tol: float = 1e-6,
                                max_iter: int = 30) -> PairResult:
        """Both blade angles fixed: the rpm (front; rear = gear_ratio times it) at which the
        pair absorbs the shaft power (W). Secant iteration on rpm."""
        power = float(power)
        count, hist = [0], []

        def p_of(rpm):
            res = self.solve_rpm(rpm)
            count[0] += res.iterations
            if not res.converged:
                raise PairConvergenceError("coupled iteration did not converge")
            hist.append(abs(res.power - power) / abs(power))
            return res.power - power

        r0 = float(rpm_guess)
        r1 = 1.05 * r0
        f0, f1 = p_of(r0), p_of(r1)
        for _ in range(max_iter):
            if abs(f1) <= tol * abs(power):
                return self._result(True, count[0], hist)
            r0, r1 = r1, r1 - f1 * (r1 - r0) / (f1 - f0)
            f0, f1 = f1, p_of(r1)
        return self._result(False, count[0], hist)

    def _guess_torque(self, thrust, rpm):
        """Starting total torque for solve_thrust: momentum-theory power for the thrust."""
        s = self.front.s
        area = PI * s.rad * s.rad
        vi = 0.5 * (-s.vel + math.sqrt(s.vel * s.vel + 2.0 * abs(thrust) / (s.rho * area)))
        p = abs(thrust) * (s.vel + vi) / 0.85
        return p / (float(rpm) * PI / 30.0)

    # -------------------------------------------------------------------- results
    def _totals(self):
        pf, pr = self.front.performance(), self.rear.performance()
        return pf["thrust"] + pr["thrust"], pf["power"] + pr["power"], pf, pr

    def field(self, r, z):
        """Combined induced velocity of the pair at (r, z): (axial, radial, swirl) m/s,
        z from the front disk, swirl positive in the front rotor's direction."""
        ff, fr = self._field(self.front), self._field(self.rear)
        if ff is None or fr is None:
            raise ValueError("both rotors need a converged solution")
        r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
        a = ff.velocity(r, z)
        b = fr.velocity(r, z - self.spacing)
        return a[0] + b[0], a[1] + b[1], a[2] + self.sign * b[2]

    def _result(self, ok, n, hist) -> PairResult:
        t, p, pf, pr = self._totals()
        s = self.front.s
        qf, qr = pf["torque"], pr["torque"]
        net = qf + self.sign * qr
        # far-wake swirl of the pair at the larger rotor's radii (no contraction modeled)
        rmax = max(self.front.s.rad, self.rear.s.rad)
        r = np.linspace(0.0, rmax, 401)[1:]
        zfar = 1.0e3 * rmax
        u, _, v = self.field(r, zfar)
        dr = r[1] - r[0]
        swirl_flux = float(np.add.reduce(0.5 * s.rho * v * v * (s.vel + u) * 2.0 * PI * r * dr))
        return PairResult(
            converged=bool(ok), iterations=int(n), front=pf, rear=pr,
            thrust=t, power=p, efficiency=(t * s.vel / p) if p else float("nan"),
            torque_split=abs(qf) / (abs(qf) + abs(qr)) if (qf or qr) else float("nan"),
            net_torque=net, swirl_power_fraction=swirl_flux / p if p else float("nan"),
            residual_swirl={"r": r, "v": v}, history=list(hist))
