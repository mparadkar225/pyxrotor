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
imposed velocities stop changing.

interaction="weights" replaces the field with XROTOR's own coupling, as CROTOR automates
it (Philip Carter's counter-rotation extension, crotor.f v0.85 in CROTOR 755-ES1.3, GPL v2
or later; a port, checked bit for bit against a transcription in the reference shim): each
rotor's VPUT slipstream (far-downstream velocities at the disk radii) is imposed on the
other with fixed weights, the default being CROTOR's close-spacing set (front: axial 0.5,
swirl 0; rear: axial 0.5, swirl -1 for counter-rotation). The spacing does not enter. The
rotors are solved alternately, front first, starting from the profiles they already carry
(CROTOR's LFWD/LAFT take them from the restart files), with no relaxation, until each
rotor's thrust changes by at most thrust_tol from its previous solve, and then each is
solved once more (CROTOR's RPM command). These weights are the field model's zero-spacing
limit (tests/model/test_pair.py). Operating conditions (RotorPair.solve_*):
    solve_rpm(rpm)                         both blade angles fixed, rpm given
    solve_power(power, rpm, torque_split)  rpm and shaft power given, pitch of each rotor
                                           set to give the torque split (0.5: equal torque)
    solve_thrust(thrust, rpm, torque_split) the same with total thrust given
    solve_power_fixed_pitch(power, rpm)    both blade angles fixed, rpm found for the power
    solve_power_split(power, rpm, power_split, fix, locked)
                                           shaft power split by power (CROTOR's POWE): at
                                           fixed rpm, at fixed pitch with each rotor's own
                                           rpm, or at fixed pitch geared together
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
from .fmath import div
from .numerics import spline
from .oper import aper
from .slipstream import set_profile
from .wake import uvadd


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
    profile set before (in the "weights" interaction, starting from it).

    interaction: "field" (the induced field at the spacing; relaxation, tol, max_iter
    apply) or "weights" (CROTOR's VPUT/VGET coupling; weights, thrust_tol, max_solves
    apply). weights = (front axial, front swirl, rear axial, rear swirl), CROTOR's UWTF,
    VWTF, UWTA, VWTA; default (0.5, 0.0, 0.5, -1.0), with the rear swirl weight +1.0 for a
    co-rotating pair. thrust_tol and max_solves are CROTOR's CCON (0.001) and CRIT (20,
    counted in single-rotor solves).
    """

    def __init__(self, front, rear, spacing: float, counter_rotating: bool = True,
                 gear_ratio: float = 1.0, relaxation: float = 0.7, tol: float = 1e-6,
                 max_iter: int = 60, interaction: str = "field", weights=None,
                 thrust_tol: float = 1e-3, max_solves: int = 20):
        self.front, self.rear = front, rear
        self.spacing = float(spacing)
        if self.spacing < 0.0:
            raise ValueError("spacing must be >= 0 (rear behind front)")
        self.sign = -1.0 if counter_rotating else 1.0
        if interaction not in ("field", "weights"):
            raise ValueError("interaction must be 'field' or 'weights'")
        self.interaction = interaction
        if weights is None:
            weights = (0.5, 0.0, 0.5, self.sign)
        self.weights = tuple(float(w) for w in weights)
        if len(self.weights) != 4:
            raise ValueError("weights = (front axial, front swirl, rear axial, rear swirl)")
        self.thrust_tol = float(thrust_tol)
        self.max_solves = int(max_solves)
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
        if self.interaction == "weights":
            return self._iterate_weights(solve_front, solve_rear)
        return self._iterate_field(solve_front, solve_rear)

    @staticmethod
    def _putvel(rotor, onto, axial_weight: float, swirl_weight: float) -> None:
        """Port of CROTOR's SUBROUTINE PUTVEL (crotor.f): the solved rotor's slipstream
        (SAVVEL's far-downstream induced velocities at its disk radii), weighted and splined,
        becomes the other rotor's imposed profile (CROTOR stores it, and LOADCR loads it
        before that rotor's next solve). Same values as slipstream.profile except for a
        ducted rotor, where PUTVEL forms VA*VADUCT_VA rather than VA + VD.

        CROTOR: PUTVEL sets NADD = II before its UVADD calls, so the solved rotor's own
        inflow is read over its first II profile entries, whatever its profile length (stale
        entries after VCLR, which zeroes only the counts). WA and WT cancel from the result
        algebraically, so this only changes rounding; reproduced, and left in the solved
        rotor's State as in the COMMON."""
        s, o = rotor.s, onto.s
        ii = s.ii
        s.nadd = ii                                    # CROTOR: NADD=II
        blds = float(s.nblds)
        uduct, vaduct_va = 0.0, 1.0
        if s.duct:
            uduct = s.urduct - 1.0
            vaduct_va = 2.0 * s.urduct
        r = np.empty(ii)
        u = np.empty(ii)
        v = np.empty(ii)
        for i in range(ii):
            xi = float(s.xi[i])
            vt = div(blds * float(s.gam[i]), (4.0 * PI) * xi)
            va = div(vt * xi, s.adw)
            va = va * vaduct_va
            utot = 1.0 + uduct + float(s.ubody[i])
            wa, wt = uvadd(s, xi)
            ci = div(xi, s.adv) - wt - vt
            si = utot + wa + va
            r[i] = xi * s.rad                                  # RDIM
            u[i] = (2.0 * (si - (utot + wa))) * s.vel          # UDIM
            v[i] = (-2.0 * (ci - (div(xi, s.adv) - wt))) * s.vel   # VDIM
        o.nadd = ii
        o.radd[:ii] = r
        o.uadd[:ii] = u * axial_weight
        o.vadd[:ii] = v * swirl_weight
        o.uaddr[:ii] = spline(o.uadd[:ii], o.radd[:ii], ii)
        o.vaddr[:ii] = spline(o.vadd[:ii], o.radd[:ii], ii)

    @staticmethod
    def _shared_tail(last, rotor) -> None:
        """CROTOR: both rotors share XROTOR's COMMON, and LOADCR copies only the first NADD
        entries of a rotor's stored profile into it, so the entries beyond come from the
        rotor in the COMMON before (the last one solved, or loaded). PUTVEL reads them
        (NADD=II). Reproduced by copying them into this rotor's State. The same goes for
        ADV, which the fixed-pitch power solve (APER(3,2)) starts from: it is the previous
        rotor's."""
        if last is rotor:
            return
        a, b = last.s, rotor.s
        n = b.nadd
        for k in ("radd", "uadd", "vadd", "uaddr", "vaddr"):
            getattr(b, k)[n:] = getattr(a, k)[n:]
        b.adv = a.adv                               # CROTOR: ADV is COMMON

    @staticmethod
    def _restore_pitch(rotor) -> None:
        """CROTOR label 590 (POWE): a failed solve's blade-angle change is undone."""
        s = rotor.s
        s.beta[:s.ii] = s.beta[:s.ii] - s.dbeta
        s.beta0[:s.ii] = s.beta0[:s.ii] - s.dbeta

    def _iterate_weights(self, solve_front, solve_rear, powe=None) -> tuple[bool, int, list]:
        """Port of CROTOR's CR iteration (SUBROUTINE CROTOR, labels 410-520). Solves front
        and rear alternately (ICRITR counts the solves), each with the profile the other
        last left it (PUTVEL). A rotor is converged (LCONF/LCONA) once its thrust changes
        by at most thrust_tol relative to its previous solve; then the system converges
        (LCON) after one more solve of each, front first, ending on the rear. At max_solves
        without that, it stops the same way (LNCON). Returns (LCON, ICRITR, the relative
        thrust change after each solve). CROTOR's convergence menu (label 590) is a
        PairConvergenceError here, after POWE's restore of the blade angles.

        powe: None for solves given as functions (the RPM command, or this package's
        torque-split modes), else the POWE command (label 430): a dict with "power" (W,
        POWINF, POWINA), "rpm" (RPMINF, RPMINA), "pitch" (LPITCH: blade angles fixed) and
        "locked" (LSRPM: rpms synchronized at the gear ratio).
          pitch False: each rotor at its rpm, blade angle set for its power, APER(3,1);
              a converged solve with a blade-angle change over 12 deg and negative thrust
              is a failure.
          pitch True, locked False: blade angles fixed, each rotor's rpm found for its
              power, APER(3,2), starting from the ADV in the COMMON.
          pitch True, locked True: the same until either rotor has converged; then (LSWIT)
              both run at fixed pitch and rpm, the front at (RPMFWD + RPMAFT) / (1 + gear
              ratio) and the rear at gear ratio times that, APER(4,2), the convergence flags
              cleared once. Total power is then only near the target (CROTOR's "(Near)
              specified power")."""
        uwtf, vwtf, uwta, vwta = self.weights
        ttot_prev = {"front": 0.0, "rear": 0.0}     # TTOTF, TTOTA
        flag = {"front": False, "rear": False}      # LCONF, LCONA
        lcon = lncon = False
        lswit = False                               # LSWIT
        rpm_found = {"front": 0.0, "rear": 0.0}     # RPMFWD, RPMAFT
        rpm_locked = {"front": 0.0, "rear": 0.0}    # RPMTF, RPMTA
        icritr = 0
        forward = False                             # LFORWARD
        hist = []
        last = self.rear                            # in the COMMON after LFWD, LAFT
        while True:
            icritr += 1
            forward = not forward
            which = "front" if forward else "rear"
            rotor, other = (self.front, self.rear) if forward else (self.rear, self.front)
            self._shared_tail(last, rotor)
            last = rotor
            if powe is None:
                solve = solve_front if forward else solve_rear
                if not solve(rotor):
                    raise PairConvergenceError(f"{which} rotor did not converge (solve {icritr})")
            else:
                s = rotor.s
                s.pspec = powe["power"][0 if forward else 1]      # PSPEC
                rpm = powe["rpm"][0 if forward else 1]
                s.conv = False
                if not powe["pitch"]:
                    s.adv = div(s.vel, s.rad * rpm * PI / 30.0)
                    aper(s, 3, 1, s.loprini)
                elif not powe["locked"]:
                    aper(s, 3, 2, s.loprini)
                    rpm_found[which] = div(s.vel, s.rad * s.adv * PI / 30.0)
                else:
                    if (flag["rear"] or flag["front"]) and not lswit:
                        rpm_locked["front"] = div(rpm_found["front"] + rpm_found["rear"],
                                                  1.0 + self.gear_ratio)
                        rpm_locked["rear"] = rpm_locked["front"] * self.gear_ratio
                        lswit = True
                        flag["front"] = flag["rear"] = lcon = False
                    if lswit:
                        s.adv = div(s.vel, s.rad * rpm_locked[which] * PI / 30.0)
                        s.conv = False
                        aper(s, 4, 2, s.loprini)
                    else:
                        aper(s, 3, 2, s.loprini)
                        rpm_found[which] = div(s.vel, s.rad * s.adv * PI / 30.0)
                bad = not s.conv
                if s.conv and not powe["pitch"]:
                    bad = abs(s.dbeta) > PI / 15.0 and s.ttot < 0.0
                if bad:
                    self._restore_pitch(rotor)
                    why = "did not converge" if not s.conv else \
                        "negative thrust and excessive blade angle change"
                    raise PairConvergenceError(f"{which} rotor {why} (solve {icritr})")
            if forward:
                self._putvel(rotor, other, uwta, vwta)
            else:
                self._putvel(rotor, other, uwtf, vwtf)
            ttot = rotor.s.ttot
            tdiff = abs(div(ttot - ttot_prev[which], ttot))
            ttot_prev[which] = ttot
            hist.append(tdiff)
            if tdiff <= self.thrust_tol:
                flag[which] = True
            if flag["front"] and flag["rear"]:
                if powe is not None and powe["pitch"] and powe["locked"] and not lswit:
                    continue                        # on to the synchronized phase
                # label 500
                if not forward and not lcon:
                    continue
                lcon = True
            elif icritr >= self.max_solves:         # label 510
                if not forward and not lncon:
                    continue
                lncon = True
            else:
                continue
            if forward:                             # label 520: the aft rotor last
                continue
            return lcon, icritr, hist

    def _iterate_field(self, solve_front, solve_rear) -> tuple[bool, int, list]:
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

    def solve_power_split(self, power: float, rpm: float | None = None,
                          power_split: float = 0.5, fix: str = "rpm",
                          locked: bool = True) -> PairResult:
        """Shaft power (W), split between the rotors by power (front share power_split),
        as CROTOR's POWE command does (equal halves by default):
          fix="rpm": the front at rpm, the rear at gear_ratio * rpm; each blade angle set
              to absorb its share (XROTOR POWE, fixed rpm).
          fix="pitch", locked=False: blade angles fixed, each rotor turns at the rpm that
              absorbs its share (rotors not geared together; rpm unused).
          fix="pitch", locked=True: rotors geared together at gear_ratio (CROTOR's SYNC):
              CROTOR's procedure, which first finds each rotor's rpm for its share, then
              runs both at their gear-ratio average and converges at fixed pitch, so the
              total power is only near the target ("weights" interaction only; for the
              exact fixed-pitch power use solve_power_fixed_pitch).
        In the "weights" interaction this is CROTOR's POWE (FIX, SYNC) bit for bit."""
        if not 0.0 < power_split < 1.0:
            raise ValueError("power_split must be between 0 and 1")
        fix = fix.lower()
        if fix not in ("rpm", "pitch"):
            raise ValueError("fix must be 'rpm' or 'pitch'")
        if fix == "rpm" and rpm is None:
            raise ValueError("fix='rpm' needs the rpm")
        pf, pa = power_split * float(power), (1.0 - power_split) * float(power)  # POWINF, POWINA
        rf = float(rpm) if rpm is not None else 0.0
        powe = {"power": (pf, pa), "rpm": (rf, self._rpm_rear(rf)),
                "pitch": fix == "pitch", "locked": bool(locked)}
        if self.interaction == "weights":
            ok, n, h = self._iterate_weights(None, None, powe)
            return self._result(ok, n, h)
        if powe["pitch"] and powe["locked"]:
            raise NotImplementedError("locked rpm at fixed pitch: use solve_power_fixed_pitch "
                                      "(or interaction='weights' for CROTOR's procedure)")

        def solver(p, r):
            def solve(rotor):
                s = rotor.s
                s.pspec = p
                s.conv = False
                if powe["pitch"]:
                    aper(s, 3, 2, s.loprini)
                else:
                    s.adv = div(s.vel, s.rad * r * PI / 30.0)
                    aper(s, 3, 1, s.loprini)
                return s.conv
            return solve
        ok, n, h = self._iterate(solver(pf, powe["rpm"][0]), solver(pa, powe["rpm"][1]))
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
        if self.front.s.duct or self.rear.s.duct:      # InducedField has no duct model
            v = np.full_like(r, np.nan)
            swirl_flux = float("nan")
        else:
            u, _, v = self.field(r, zfar)
            dr = r[1] - r[0]
            swirl_flux = float(np.add.reduce(0.5 * s.rho * v * v * (s.vel + u) * 2.0 * PI * r * dr))
        return PairResult(
            converged=bool(ok), iterations=int(n), front=pf, rear=pr,
            thrust=t, power=p, efficiency=(t * s.vel / p) if p else float("nan"),
            torque_split=abs(qf) / (abs(qf) + abs(qr)) if (qf or qr) else float("nan"),
            net_torque=net, swirl_power_fraction=swirl_flux / p if p else float("nan"),
            residual_swirl={"r": r, "v": v}, history=list(hist))
