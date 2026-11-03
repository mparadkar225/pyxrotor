"""pyxrotor: a Python port of Mark Drela and Harold Youngren's XROTOR 7.55.

    from pyxrotor import Rotor, AeroSection
    r = Rotor()                       # XROTOR start-up state
    r.atmosphere(0.0)
    r.set_sections([AeroSection(...)])
    r.arbi(nblds, vel, tip_radius, hub_radius, r_over_R, c_over_R, beta_deg)
    r.velocity(40.0)
    r.solve_rpm(1200.0)               # or solve_advance_ratio / solve_thrust / ...
    r.performance()
    r.design(nblds, tip_radius, hub_radius, vel, rpm=..., power=..., cl=...)   # DESI

Modules (each a port of the named XROTOR routines, bit-identical to 7.55):
    state      State: all 27 COMMON blocks
    rotor      INIT, SETDEF, ATMO, ARBI, OPER and DESI commands (Rotor)
    design     SETCLD, DESINI, DEITER, DESMIL, DESGEN
    cases      SETCAS, CASPLT (solves), SHOCAS, GETCAS, ADDC/CLRC/CASE
    restart    SAVE, LOAD
    nacelle    NACELL (NACE), SRCVEL
    fortio     gfortran formatted output and list-directed input, GETFLT
    oper       TPQ, APINIT, APITER, APER
    induced    VCALC, GRADMO, HELICO, VRTXCO, VORSEGVEL
    wake       UVADD, CSCALC, XWINIT, SETXW
    aero       CLCDCM, GETCLCDCM, GETALF, SETIAERO, SORTAR
    numerics   GAUSS, TRISOL, SPLINE, SEVAL, DEVAL, SEVLIN, FILTER
    geometry   SETX, TINVRT
    fmath      Fortran IEEE semantics for Python scalar math

See CLAUDE.md for the porting rules and README.md for status.
"""

from .aero import AeroSection
from .rotor import Rotor
from .state import State

__version__ = "0.1.0"
__all__ = ["AeroSection", "Rotor", "State"]
