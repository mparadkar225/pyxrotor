"""User-facing XROTOR session: start-up, blade input, and the OPER menu commands.

    from pyxrotor import Rotor
    r = Rotor()                                  # XROTOR start-up state (INIT, SETDEF)
    r.atmosphere(0.0)                            # top-level ATMO, km
    r.set_sections([...])                        # AERO menu: parametric sections
    r.arbi(nblds=2, vel=50.0, tip_radius=1.64, hub_radius=0.33,
           r_over_R=..., c_over_R=..., beta_deg=...)
    r.velocity(30.0)                             # OPER VELO
    r.solve_rpm(1000.0)                          # OPER RPM
    r.performance()                              # thrust, power, efficiency, CT, CP, ...

Each method is a port of the corresponding XROTOR routine or OPER command handler,
with its terminal prompts turned into arguments; the State (r.s) holds everything.

Ported here: INIT, SETDEF, ATMO, ARBI, the OPER commands ATMO, VELO, ANGL, ADVA, RPM,
THRU, TORQ, POWE, ITER, INIT (toggle), and the DESI commands INPU, EDIT (B, RT, RH, RW, V,
A, R, T, P, CC, CL), FORM, ITER, N (design, redesign, ...), the OPER sweeps ASEQ, RSEQ,
BSEQ, VSEQ and case commands ADDC, CLRC, CASE, CPUT, CGET, and top-level SAVE and LOAD.
Not yet: DESI CR/CW/CX (CL from file, to file, by cursor), sequences (ASEQ/RSEQ/BSEQ/VSEQ, executed in
xcasepl.f) are ported here; still missing: PVAR file reading, CPUT/CGET, the xrotor.def file (GETDEF/WDEF).
"""
from __future__ import annotations

from dataclasses import asdict
from typing import Sequence

import numpy as np

from . import fmath as fm
from .aero import AERODATA_ROW, AeroSection, setiaero, sort_sections
from .constants import IX, PI
from .fmath import div
from .geometry import setx, tinvrt
from .numerics import seval, spline
from .design import BetzLimitError, desgen, setcld
from .oper import aper
from .state import State
from .wake import xwinit

# ---- ATMO tables (U.S. Standard Atmosphere)
_ALT = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
        20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
        30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
        40.0, 45.0, 60.0, 75.0]
_VSO = [340.0, 336.0, 332.0, 329.0, 325.0, 320.0, 316.0, 312.0, 308.0, 304.0,
        299.0, 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, 295.0, 295.0,
        295.0, 295.8, 296.4, 297.1, 297.8, 298.5, 299.1, 299.8, 300.5, 301.1,
        301.8, 302.5, 303.1, 305.0, 306.8, 308.7, 310.5, 312.3, 314.0, 316.0,
        318.0, 355.0, 372.0, 325.0]
_RHO = [1.226, 1.112, 1.007, 0.909, 0.820, 0.737, 0.660, 0.589, 0.526, 0.467,
        0.413, 0.364, 0.311, 0.265, 0.227, 0.194, 0.163, 0.141, 0.121, 0.103,
        .0880, .0749, .0637, .0543, .0463, .0395, .0338, .0288, .0246, .0210,
        .0180, .0154, .0132, .0113, .0096, .0082, .0070, .0060, .0052, .0044,
        0.004, 0.002, 3.9E-4, 8.0E-5]
_RMU = [1.780, 1.749, 1.717, 1.684, 1.652, 1.619, 1.586, 1.552, 1.517, 1.482,
        1.447, 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, 1.418, 1.418,
        1.418, 1.427, 1.433, 1.438, 1.444, 1.449, 1.454, 1.460, 1.465, 1.471,
        1.476, 1.481, 1.487, 1.502, 1.512, 1.532, 1.546, 1.561, 1.580, 1.600,
        1.700, 1.912, 2.047, 1.667]


def atmo(alspec: float) -> tuple[float, float, float]:
    """Port of SUBROUTINE ATMO: (speed of sound m/s, density kg/m^3, viscosity kg/m-s) at
    altitude alspec km, linearly interpolated; alspec = -1 gives water at 15 C."""
    alspec = float(alspec)
    if alspec == -1.0:
        return 1500.0, 1000.0, 1.15e-3
    for i in range(1, len(_ALT)):
        if alspec > _ALT[i]:
            continue
        dalt = _ALT[i] - _ALT[i - 1]
        alfrac = (alspec - _ALT[i - 1]) / dalt
        vso = _VSO[i - 1] + (_VSO[i] - _VSO[i - 1]) * alfrac
        rho = _RHO[i - 1] + (_RHO[i] - _RHO[i - 1]) * alfrac
        rmu = (_RMU[i - 1] + (_RMU[i] - _RMU[i - 1]) * alfrac) * 1.0e-5
        return vso, rho, rmu
    if alspec > _ALT[-1]:                          # 'You're in low earth orbit.  Good luck.'
        return _VSO[-1], _RHO[-1], _RMU[-1] * 1.0e-5
    return 0.0, 0.0, 0.0                           # NaN altitude: outputs left unset


def put_aero(s: State, n: int, sec: AeroSection) -> None:
    """Port of SUBROUTINE PUTAERO (n is 0-based)."""
    for name, row in AERODATA_ROW.items():
        s.aerodata[row - 1, n] = getattr(sec, name)
    s.xiaero[n] = sec.xisect


def startup() -> State:
    """XROTOR's state after program start: PROGRAM XROTOR's set-up, INIT and SETDEF.
    (GETDEF's optional xrotor.def file is not read.)"""
    s = State()
    s.version = 7.55
    s.luread, s.luwrit, s.lutemp, s.lusave = 5, 6, 3, 4
    # ---- INIT
    s.greek = False
    s.idev, s.idevrp, s.ipslu = 1, 2, 0
    s.scrnfr, s.size = 0.70, 10.0
    s.xpage, s.ypage, s.xmarg, s.ymarg = 11.0, 8.5, 0.0, 0.0
    s.legend = False
    s.ncolor = 64
    s.lslope = False
    s.urduct = 1.0
    # ---- SETDEF
    s.rake = 0.0
    s.vel = 1.0
    s.alt = 0.0
    s.vso, s.rho, s.rmu = atmo(s.alt)
    s.naero = 1
    put_aero(s, 0, AeroSection())
    s.iaero[:] = 0                                # IAERO(I) = 1 (0-based 0)
    s.xpitch = 0.3
    s.ii = 30
    s.incr = 2
    s.ixspac = 2
    s.vrtx, s.fast, s.free, s.duct, s.terse = False, False, True, False, False
    s.lland, s.lgrid, s.lvnorm = True, True, True
    s.par, s.csize = 0.6, 0.014
    # ---- rest of INIT
    s.xinf = 3.0
    s.nn = 32
    s.iinf = s.ii + s.ii // 2
    s.conv = False
    s.lplot = False
    s.wind = False
    s.lstruc = False
    s.name = ""
    s.savfil = ""
    s.gee = 9.81
    s.adwfctr = 1.0
    s.xitip = 1.0
    s.ubody[:s.ii] = 0.0
    s.nadd = 0
    s.ncase = 0
    s.kcase = 0
    s.niterd = 40
    s.nitera = 40
    s.ldesini = False
    s.loprini = True
    s.lpwrvar = False
    s.npwrvar = 0
    s.lrotor = False
    s.iaero[:] = -1                               # IAERO(I) = 0: unset
    s.xyoff[:] = 0.0
    s.xyfac[:] = 0.0
    return s


def arbi(s: State, nblds: int, vel: float, tip_radius: float, hub_radius: float,
         r_over_R: Sequence[float], c_over_R: Sequence[float], beta_deg: Sequence[float]) -> None:
    """Port of SUBROUTINE ARBI: arbitrary blade geometry from stations of r/R, chord/R and
    blade angle (deg), splined onto XROTOR's radial grid. A negative tip_radius means r and
    c are given in meters, as in the Fortran."""
    x = np.asarray(r_over_R, dtype=float).copy()
    c = np.asarray(c_over_R, dtype=float).copy()
    b = np.asarray(beta_deg, dtype=float).copy()
    nst = x.size
    if nst < 3:
        raise ValueError("ARBI: must have 3 or more stations")
    s.greek = False
    s.conv = False
    s.nblds = int(nblds)
    s.vel = float(vel)
    rad = abs(float(tip_radius))
    if rad <= hub_radius:
        raise ValueError("ARBI: tip radius must exceed hub radius")
    s.rad = rad
    s.xi0 = div(hub_radius, rad)
    s.xw0 = div(hub_radius, rad)
    grid = setx(s.ii, s.xi0, s.xitip, s.ixspac)
    s.dt = grid.dt
    s.t[:s.ii], s.xi[:s.ii], s.dxi[:s.ii], s.xv[:s.ii + 1] = grid.t, grid.xi, grid.dxi, grid.xv
    xwinit(s)

    w1, w2, w3 = x, c, b
    if tip_radius < 0.0:
        w1 = w1 / rad
        w2 = w2 / rad
    for n in range(nst):
        if w1[n] < s.xi0 or w1[n] > s.xitip:
            raise ValueError(f"ARBI: radii must be between {s.xi0} and 1.0")
        if n > 0 and w1[n] <= w1[n - 1]:
            raise ValueError("ARBI: radii must monotonically increase")
    w3 = w3 * PI / 180.0
    w1 = np.array([tinvrt(v, s.xi0, s.xitip, s.ixspac) for v in w1])
    w4 = spline(w2, w1, nst)
    w5 = spline(w3, w1, nst)
    for i in range(s.ii):
        s.ch[i] = abs(seval(s.t[i], w2, w4, w1, nst))
        s.beta[i] = seval(s.t[i], w3, w5, w1, nst)
        s.beta0[i] = s.beta[i]
        s.ubody[i] = 0.0
    if s.name == "":
        s.name = "Arbitrary blade"

    # ---- reasonable advance ratio to start the iterative routines
    is_ = s.ii // 2                               # Fortran IS = II/2 + 1
    a0 = 0.0
    ang = float(s.beta[is_]) - a0
    s.adv = div(float(s.xi[is_]) * fm.sin(ang), fm.cos(ang))
    s.adv = fm.fmax(0.1, s.adv)
    s.adw = s.adv
    s.iaero[:s.ii] = setiaero(s.xi[:s.ii], sections_of(s))


def sections_of(s: State) -> list[AeroSection]:
    out = []
    for n in range(s.naero):
        kw = {name: float(s.aerodata[row - 1, n]) for name, row in AERODATA_ROW.items()}
        out.append(AeroSection(xisect=float(s.xiaero[n]), **kw))
    return out


class Rotor:
    """An XROTOR session. r.s is the State (all COMMON variables).

    refresh_design_sections: False (default) reproduces 7.55, where a design keeps the
    aero-section indices (IAERO) of the previous radial grid, so with several sections the
    outer stations extrapolate past the last one (see CLAUDE.md). True assigns sections on
    the design's own grid (SETIAERO after SETX in DESGEN).
    """

    def __init__(self, refresh_design_sections: bool = False):
        self.s = startup()
        self.refresh_design_sections = bool(refresh_design_sections)

    # ------------------------------------------------------------- top level / AERO
    def atmosphere(self, alt_km: float) -> None:
        """ATMO r: fluid properties from the standard atmosphere (alt -1: water)."""
        s = self.s
        s.alt = float(alt_km)
        s.vso, s.rho, s.rmu = atmo(s.alt)

    def set_sections(self, sections: Sequence[AeroSection]) -> None:
        """AERO menu: replace the section list (sorted by r/R, as SORTAR does)."""
        s = self.s
        secs = sort_sections(list(sections))
        if not 1 <= len(secs) <= s.xiaero.size:
            raise ValueError(f"need 1 to {s.xiaero.size} sections")
        for n, sec in enumerate(secs):
            put_aero(s, n, sec)
        s.naero = len(secs)
        if not s.conv:
            s.iaero[:s.ii] = setiaero(s.xi[:s.ii], secs)

    def arbi(self, nblds, vel, tip_radius, hub_radius, r_over_R, c_over_R, beta_deg) -> None:
        """ARBI: arbitrary rotor geometry."""
        arbi(self.s, nblds, vel, tip_radius, hub_radius, r_over_R, c_over_R, beta_deg)

    # ------------------------------------------------------------------- OPER menu
    def velocity(self, vel: float) -> None:
        """VELO r: change flight speed, keeping dimensional thrust, torque and power."""
        s = self.s
        velold = s.vel
        s.vel = float(vel)
        r2, r3 = s.rad * s.rad, s.rad * s.rad * s.rad          # RAD**2, RAD**3
        vo2, v2 = velold * velold, s.vel * s.vel               # VEL**2
        thr = s.ttot * (s.rho * vo2 * r2)
        s.ttot = div(thr, s.rho * v2 * r2)
        trq = s.qtot * (s.rho * vo2 * r3)
        s.qtot = div(trq, s.rho * v2 * r3)
        pwr = s.ptot * (s.rho * (vo2 * velold) * r2)
        s.ptot = div(pwr, s.rho * (v2 * s.vel) * r2)
        s.conv = False

    def change_pitch(self, delta_deg: float) -> None:
        """ANGL r: change the collective blade angle by delta_deg."""
        s = self.s
        d = float(delta_deg) * PI / 180.0
        s.beta[:s.ii] = s.beta[:s.ii] + d
        s.beta0[:s.ii] = s.beta0[:s.ii] + d
        s.conv = False

    def solve_advance_ratio(self, adv: float) -> bool:
        """ADVA r: prescribe advance ratio (V/(Omega R)) and solve."""
        s = self.s
        s.adv = float(adv)
        s.conv = False
        aper(s, 4, 2, s.loprini)
        return s.conv

    def solve_rpm(self, rpm: float) -> bool:
        """RPM r: prescribe rpm and solve."""
        s = self.s
        s.adv = div(s.vel, s.rad * float(rpm) * PI / 30.0)
        s.conv = False
        aper(s, 4, 2, s.loprini)
        return s.conv

    def _solve_target(self, ispec: int, fix: str, rpm: float | None) -> bool:
        s = self.s
        fix = fix.upper()[:1]
        if fix not in ("P", "R"):
            raise ValueError("fix must be 'pitch' or 'rpm'")
        s.conv = False
        if fix == "P":
            aper(s, ispec, 2, s.loprini)
        else:
            if rpm is None:
                rpm = div(s.vel, s.rad * s.adv * PI / 30.0)
            s.adv = div(s.vel, s.rad * float(rpm) * PI / 30.0)
            aper(s, ispec, 1, s.loprini)
            if not s.conv:
                # ---- convergence failed: restore clobbered blade angles
                s.beta[:s.ii] = s.beta[:s.ii] - s.dbeta
                s.beta0[:s.ii] = s.beta0[:s.ii] - s.dbeta
        return s.conv

    def solve_thrust(self, thrust: float, fix: str = "pitch", rpm: float | None = None) -> bool:
        """THRU r: solve for thrust (N) at fixed pitch (rpm free) or fixed rpm (pitch free)."""
        self.s.tspec = float(thrust)
        return self._solve_target(1, fix, rpm)

    def solve_torque(self, torque: float, fix: str = "pitch", rpm: float | None = None) -> bool:
        """TORQ r: solve for torque (N m)."""
        self.s.qspec = float(torque)
        return self._solve_target(2, fix, rpm)

    def solve_power(self, power: float, fix: str = "pitch", rpm: float | None = None) -> bool:
        """POWE r: solve for power (W)."""
        self.s.pspec = float(power)
        return self._solve_target(3, fix, rpm)

    def set_iterations(self, n: int) -> None:
        """ITER i: maximum Newton iterations."""
        self.s.nitera = int(n)

    def toggle_initialization(self) -> bool:
        """INIT: toggle re-initialization of each analysis case (LOPRINI)."""
        self.s.loprini = not self.s.loprini
        return self.s.loprini

    # ------------------------------------------------------------------- DESI menu
    def design(self, nblds: int, tip_radius: float, hub_radius: float, vel: float, *,
               wake_hub_radius: float = 0.0, adv: float | None = None, rpm: float | None = None,
               thrust: float | None = None, power: float | None = None, cl: float = 0.0,
               betz_ok: bool = False) -> bool:
        """DESI INPU (then <Return>): design a minimum-induced-loss rotor.

        Give adv or rpm, and thrust (N) or power (W). cl is the constant design CL; 0 means
        XROTOR's default, 0.8 with the sign of the thrust or power. Returns CONV.
        Raises BetzLimitError where XROTOR would ask whether to continue past the Betz
        windmill limits and betz_ok is False; ValueError where DESGEN rejects the inputs.
        """
        s = self.s
        s.nblds = int(nblds)
        s.raddes, s.r0des, s.rwdes = float(tip_radius), float(hub_radius), float(wake_hub_radius)
        s.veldes = float(vel)
        if (adv is None) == (rpm is None):
            raise ValueError("give exactly one of adv and rpm")
        s.advdes = 0.0 if adv is None else float(adv)
        s.rpmdes = 0.0 if rpm is None else float(rpm)
        if s.advdes == 0.0 and rpm is None:
            raise ValueError("adv = 0 prescribes rpm: give rpm instead")
        if (thrust is None) == (power is None):
            raise ValueError("give exactly one of thrust and power")
        s.tddes = 0.0 if thrust is None else float(thrust)
        s.pddes = 0.0 if power is None else float(power)
        if s.tddes == 0.0 and s.pddes == 0.0:
            raise ValueError("thrust or power must be nonzero")      # INPU asks again
        s.dest = s.tddes != 0.0
        s.desp = not s.dest
        s.cldes0 = float(cl)
        if s.cldes0 == 0.0:
            cl1 = fm.sign(0.8, s.tddes if s.dest else s.pddes)
            setcld(s, cl1, cl1)
        else:
            setcld(s, s.cldes0, s.cldes0)
        if s.name == "":
            s.name = "Designed blade"
        s.lrotor = False
        return self._desgen(betz_ok)

    def redesign(self, *, betz_ok: bool = False, **changes) -> bool:
        """DESI EDIT: change design parameters, then design again (<Return>).

        Keywords, applied in the order given, as the EDIT commands are:
            nblds (B), tip_radius (RT), hub_radius (RH), wake_hub_radius (RW), vel (V),
            adv (A; clears rpm), rpm (R; clears adv), thrust (T; design for thrust),
            power (P; design for power), cl (CC: constant CL),
            cl_root_tip=(root, tip) (CL: linear in r/R on the *current* grid, as in 7.55).
        """
        s = self.s
        for key, value in changes.items():
            if key == "nblds":
                s.nblds = int(value)
            elif key == "tip_radius":
                s.raddes = float(value)
            elif key == "hub_radius":
                s.r0des = float(value)
            elif key == "wake_hub_radius":
                s.rwdes = float(value)
            elif key == "vel":
                s.veldes = float(value)
            elif key == "adv":
                s.advdes, s.rpmdes = float(value), 0.0
            elif key == "rpm":
                s.rpmdes, s.advdes = float(value), 0.0
            elif key == "thrust":
                s.tddes, s.dest, s.desp, s.pddes = float(value), True, False, 0.0
            elif key == "power":
                s.pddes, s.dest, s.desp, s.tddes = float(value), False, True, 0.0
            elif key == "cl":
                s.cldes0 = float(value)
                setcld(s, s.cldes0, s.cldes0)
            elif key == "cl_root_tip":
                clroot, cltip = (float(v) for v in value)
                setcld(s, clroot, cltip)
                s.cldes0 = 0.5 * (clroot + cltip)
            else:
                raise TypeError(f"unknown design parameter {key!r}")
        return self._desgen(betz_ok)

    def _desgen(self, betz_ok: bool) -> bool:
        s = self.s
        s.ldesini = True                                  # 'HHY added initialization'
        sec1 = sections_of(s)[0]
        if not desgen(s, betz_ok=betz_ok, refresh_sections=self.refresh_design_sections):
            if any(c > sec1.clmax or c < sec1.clmin for c in s.cldes[:s.ii]):
                raise ValueError(f"design CL is past section 1's stall limits "
                                 f"({sec1.clmin}, {sec1.clmax})")
            raise BetzLimitError("disk loading above the ultimate Betz limits "
                                 "(pass betz_ok=True to design anyway)")
        return bool(s.conv)

    def design_formulation(self, model: str | None = None, *, toggle_free_wake: bool = False) -> None:
        """DESI FORM (also OPER FORM): 'grad' (graded momentum), 'pot' (Goldstein potential),
        'vrtx' (discrete vortex wake); toggle_free_wake is WAKE. As in 7.55, 'vrtx' leaves
        FAST alone, and graded momentum still wins while FAST is set."""
        s = self.s
        if model is not None:
            m = model.upper()
            if m == "GRAD":
                s.vrtx, s.fast = False, True
            elif m == "POT":
                s.vrtx, s.fast = False, False
            elif m == "VRTX":
                s.vrtx = True
            else:
                raise ValueError("model must be 'grad', 'pot' or 'vrtx'")
        if toggle_free_wake:
            s.free = not s.free

    def set_design_iterations(self, n: int) -> None:
        """DESI ITER (prompted form): maximum design iterations NITERD. (In 7.55, ITER with
        an inline argument, 'ITER 20', sets the OPER limit NITERA instead.)"""
        self.s.niterd = int(n)

    def set_stations(self, ii: int) -> None:
        """DESI N: change the number of radial stations, re-splining chord, blade angle,
        UBODY and the design CL onto the new grid if there is a rotor (LROTOR)."""
        s = self.s
        ii = int(ii)
        if not 1 <= ii <= IX:
            raise ValueError(f"number of stations must be 1 to {IX}")
        if s.lrotor:
            iisav = s.ii
            w1 = np.array(s.xi[:iisav])
            w2, w4, w6, w8 = (np.array(a[:iisav]) for a in (s.ch, s.beta, s.ubody, s.cldes))
            w3, w5, w7, w9 = (spline(w, w1, iisav) for w in (w2, w4, w6, w8))
        s.ii = ii
        s.iinf = ii + ii // 2
        grid = setx(s.ii, s.xi0, s.xitip, s.ixspac)
        s.dt = grid.dt
        s.t[:ii], s.xi[:ii], s.dxi[:ii], s.xv[:ii + 1] = grid.t, grid.xi, grid.dxi, grid.xv
        if s.lrotor:
            for i in range(ii):
                x = float(s.xi[i])
                s.ch[i] = seval(x, w2, w3, w1, iisav)
                s.beta[i] = seval(x, w4, w5, w1, iisav)
                s.ubody[i] = seval(x, w6, w7, w1, iisav)
                s.cldes[i] = seval(x, w8, w9, w1, iisav)
                s.beta0[i] = s.beta[i]

    # ------------------------------------------------------------- OPER case sweeps
    def _sweep(self, itype: int, first, last, step, **opts) -> list[dict]:
        from . import cases
        added = cases.setcas(self.s, itype, first, last, step, **opts)
        cases.run_cases(self.s)                        # CASPLT: solves *all* stored cases
        return self.cases()[-added:] if added else []

    def sweep_advance_ratio(self, first: float, last: float, step: float, *,
                            engine_line: bool = False) -> list[dict]:
        """ASEQ: a sequence of advance ratios at the current blade angle.

        As in 7.55, new cases are appended to the stored ones and every stored case is
        re-solved (clear_cases() first for a fresh table); afterwards the original operating
        point is restored and solved again. Returns the new cases (see cases()).
        engine_line: hold power on the engine power-rpm line (if one is set, LPWRVAR)."""
        return self._sweep(1, first, last, step, engine_line=engine_line)

    def sweep_rpm(self, first: float, last: float, step: float, *, hold: str = "pitch",
                  engine_line: bool = False) -> list[dict]:
        """RSEQ: a sequence of rpms, holding the blade 'pitch', or the current 'thrust',
        'torque' or 'power' spec (TSPEC, QSPEC, PSPEC) by changing pitch."""
        fix = {"pitch": "A", "thrust": "T", "torque": "Q", "power": "P"}[hold]
        return self._sweep(2, first, last, step, engine_line=engine_line, fix=fix)

    def sweep_blade_angle(self, first_deg: float, last_deg: float, step_deg: float, *,
                          engine_line: bool = False) -> list[dict]:
        """BSEQ: a sequence of tip blade angles (deg) at the current advance ratio.
        (7.55 leaves the step, in degrees, in the solver's DBET.)"""
        return self._sweep(3, first_deg, last_deg, step_deg, engine_line=engine_line)

    def sweep_velocity(self, first: float, last: float, step: float, *,
                       constant_speed: bool = False, rpm: float | None = None,
                       power: float | None = None, engine_line: bool = False) -> list[dict]:
        """VSEQ: a sequence of flight speeds at fixed pitch (FP: advance ratio scales with
        speed), or constant speed (CS) at rpm (default the current) and power (default
        PSPEC, or the current power), varying pitch."""
        return self._sweep(4, first, last, step, constant_speed=constant_speed, rpm=rpm,
                           power=power, engine_line=engine_line)

    def add_case(self) -> bool:
        """ADDC: store the current operating point (unsolved) in the case table."""
        from .cases import add_case
        return add_case(self.s)

    def clear_cases(self) -> None:
        """CLRC: clear the case table."""
        from .cases import clear_cases
        clear_cases(self.s)

    def rerun_case(self, n: int) -> bool:
        """CASE n (1-based): solve stored case n at its advance ratio and tip angle."""
        from .cases import rerun_case
        return rerun_case(self.s, n)

    def cases(self) -> list[dict]:
        """The case table as SHOCAS lists it (power etc. None where not solved)."""
        from .cases import NOT_SOLVED, PI_EXACT
        s = self.s
        out = []
        for n in range(s.ncase):
            p = [float(v) for v in s.caspar[:, n]]
            solved = p[8] != NOT_SOLVED
            out.append({
                "n": n + 1, "adv": p[1], "vel": p[2], "beta_tip_deg": (p[3] * 180.0) / PI_EXACT,
                "rpm": (div(p[2], s.rad * p[1]) * 30.0) / PI_EXACT,
                "alt": p[4], "rho": p[5], "rmu": p[6], "vso": p[7],
                "power": p[8] if solved else None, "thrust": p[9] if solved else None,
                "torque": p[10] if solved else None, "efficiency": p[11] if solved else None,
            })
        return out

    def write_cases(self, path) -> None:
        """CPUT f: the case table in XROTOR's listing format."""
        from .cases import case_lines
        s = self.s
        with open(path, "w") as f:
            f.write("\n".join(case_lines(s.caspar, s.ncase, s.rad, s.name)) + "\n")

    def read_cases(self, path, full_lines: bool = False) -> int:
        """CGET f: replace the case table from a CPUT file; returns the number of cases.
        As in 7.55, results (power, thrust, torque, efficiency) are lost, because lines are
        read into 128 characters; full_lines=True reads them too."""
        from .cases import cases_from_lines
        with open(path) as f:
            cases_from_lines(self.s, f.read().splitlines(), full_lines)
        return self.s.ncase

    # ------------------------------------------------------------- restart files
    def save(self, path, overwrite: bool = True) -> None:
        """SAVE f: write the rotor, aero sections and operating point (7.55 format)."""
        from .restart import save
        save(self.s, path, overwrite)

    def load(self, path, stations: int | None = None, keep_free_wake: bool = False) -> bool:
        """LOAD f: read a rotor saved by XROTOR (>= 6.9 format) or pyxrotor and solve its
        operating point. stations: the station count to use if the file's differs (7.55
        prompts; None keeps the current II). keep_free_wake=True avoids 7.55's reading of
        the file's LVDUCT flag as FREE. Returns CONV."""
        from .restart import load
        return load(self.s, path, stations, keep_free_wake)

    # --------------------------------------------------------------------- results
    def performance(self) -> dict:
        """Dimensional and coefficient results, defined as XROTOR's OUTPUT prints them."""
        s = self.s
        rv2r2 = s.rho * s.vel ** 2 * s.rad ** 2
        tdim = s.ttot * rv2r2
        qdim = s.qtot * rv2r2 * s.rad
        pdim = s.ptot * s.rho * s.vel ** 3 * s.rad ** 2
        rpm = div(s.vel, s.rad * s.adv * PI / 30.0)
        en, dia = rpm / 60.0, 2.0 * s.rad
        return {
            "converged": bool(s.conv),
            "thrust": tdim, "torque": qdim, "power": pdim,
            "efficiency": div(s.ttot, s.ptot), "induced_efficiency": div(s.twak, s.pwak),
            "rpm": rpm, "adv": s.adv, "J": s.adv * PI,
            "CT": div(tdim, s.rho * en ** 2 * dia ** 4),
            "CP": div(pdim, s.rho * en ** 3 * dia ** 5),
            "beta_change_deg": s.dbeta * 180.0 / PI,
        }

    def stations(self) -> dict:
        """Radial distributions at the II stations."""
        s = self.s
        ii = s.ii
        return {k: np.array(getattr(s, k)[:ii]) for k in ("xi", "ch", "beta", "cl", "cd", "gam", "re", "effp")}
