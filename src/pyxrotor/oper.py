"""Operating-point solver (xoper.f).

Ported routines:
    TPQ     tpq()      thrust, torque, power and their sensitivities (real and equivalent prop)
    APINIT  apinit()   initial circulation and wake advance ratio (graded momentum)
    APITER  apiter()   Newton iteration for an operating point
    APER    aper()     APINIT (optional) + APITER

ISPEC (target): 1 thrust TSPEC, 2 torque QSPEC, 3 power PSPEC, 4 fixed advance ratio,
                5 power from the engine power-rpm line.
ICON (constraint): 1 advance ratio (rpm) fixed, 2 blade pitch fixed.

All state lives in the State; outputs are the COMMON variables XROTOR sets (TTOT, PTOT,
QTOT, GAM, CL, CD, ..., CONV). Faithfully reproduced 7.55 quirks are marked "7.55:".

TPQ, APINIT and APITER are compiled: their code (with CSCALC, UVADD, SEVAL, GETCLCDCM,
GETALF and CLCDCM as C helpers) is in _speedups.pyx (Cython); the functions here call it.
"""
from __future__ import annotations

from . import _speedups
from .aero import AERODATA_ROW, AeroSection
from .constants import IQ
from .induced import gradmo, helico, vcalc, vrtxco
from .numerics import deval, filter_, seval
from .state import State
from .wake import setxw


# --------------------------------------------------------------------------- helpers
def sections(s: State) -> list[AeroSection]:
    """The aero sections stored in XIAERO/AERODATA/NAERO."""
    out = []
    for n in range(s.naero):
        kw = {name: float(s.aerodata[row - 1, n]) for name, row in AERODATA_ROW.items()}
        out.append(AeroSection(xisect=float(s.xiaero[n]), **kw))
    return out


def induced_model(s: State) -> None:
    """The APITER/DESI selection: GRADMO if FAST, else HELICO, or VRTXCO if VRTX. Sets IWTYP."""
    if s.fast:
        gradmo(s)
        s.iwtyp = 1
    elif not s.vrtx:
        helico(s)
        s.iwtyp = 2
    else:
        vrtxco(s)
        s.iwtyp = 3


# ------------------------------------------------------------------------------- TPQ
def tpq(s: State, itype: int) -> None:
    """Port of SUBROUTINE TPQ(ITYPE): thrust, power and torque, and their sensitivities to
    GAM(J), ADV, ADW and the blade-angle change.

    ITYPE 1: analysis (fixed blade angles). 2: design, fixed CL, chord from circulation.
    3: design, fixed chord, angles from CL.

    Not ported: momentum-theory checks (TMOM, PMOM, VAAavg, VATavg, TDIM, PDIM) that TPQ
    computes into unused locals, and the unused equivalent-prop axial velocity SW.
    """
    _speedups.tpq(s, itype)


# ---------------------------------------------------------------------------- APINIT
def apinit(s: State) -> None:
    """Port of SUBROUTINE APINIT: circulation from the blade-element lift with no induced
    velocity, then up to 10 graded-momentum passes converging the wake advance ratio."""
    _speedups.apinit(s, gradmo, vcalc, filter_)


# ---------------------------------------------------------------------------- APITER
def apiter(s: State, ispec: int, icon: int) -> None:
    """Port of SUBROUTINE APITER(ISPEC,ICON): Newton iteration on circulation, ADV, ADW and
    blade-angle change, up to NITERA iterations, to RMS(dGAM) <= 1e-7. Sets CONV on success.

    Not ported: the 0.75R "ADVFACT" correction, which 7.55 computes and then overwrites
    with 1.0 (the computation has no side effects).
    """
    _speedups.apiter(s, ispec, icon, IQ, induced_model, vcalc, setxw, tpq, filter_,
                      seval, deval)


def aper(s: State, ispec: int, icon: int, linit: bool) -> None:
    """Port of SUBROUTINE APER(ISPEC,ICON,LINIT): optional initialization, then APITER.
    On failure CONV stays False (XROTOR prints 'Iteration limit exceeded')."""
    if linit:
        apinit(s)
    apiter(s, ispec, icon)
