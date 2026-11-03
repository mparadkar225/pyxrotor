"""Nacelle (spinner/cowling) perturbation velocities: SUBROUTINE NACELL (xrotor.f, the
top-level NACE command) and SUBROUTINE SRCVEL (srclin.f).

The body is an axisymmetric area distribution A(z) behind (or around) the rotor. NACELL
models it with point sources on the axis, one per interval of the input stations, and
evaluates their axial and radial velocities at each blade station (UBODY, VBODY, as
fractions of the flight speed). The far-wake area sets the wake hub radius XW0.

File format (lengths in m, areas in m^2; a negative area is a radius):
    zprop  awake          axial position of the rotor, far-wake area
    z(1)   a(1)
    ...                   up to 200 stations
A file name starting with CONST (any case) instead sets a constant UBODY.

7.55 behavior reproduced (see CLAUDE.md):
  * The sources come from area changes ("HHY FIX", replacing SRCLIN's solve), but shifted by
    one interval: the source at the midpoint of interval i has the area change of interval
    i-1, the first interval gets none, and the last interval's change is dropped.
    (fix: fix_source_positions=True).
  * CONST sets XW0 = 0, discarding the hub wake radius.
  * A bad number in the file sets UBODY = 0 and returns, leaving VBODY and XW0 unchanged.
Not ported: SRCLIN (dead code in 7.55) and the plot (NACPLT).
"""
from __future__ import annotations

import math
import os

import numpy as np

from . import fmath as fm
from .constants import PI
from .fmath import div
from .fortio import ListReader
from .state import State
from .wake import xwinit

PI_EXACT = 4.0 * math.atan(1.0)          # SRCVEL does not include XROTOR.INC
IZX = 200


def srcvel(zb, src, z: float, r: float) -> tuple[float, float]:
    """Port of SUBROUTINE SRCVEL(NB,ZB,SRC,Z,R,U,V): axial and radial velocity at (z, r)
    from point sources of strength src[i] at the midpoints of the zb intervals."""
    z, r = float(z), float(r)
    usum = vsum = 0.0
    for i in range(len(zb) - 1):
        zs = (float(zb[i + 1]) + float(zb[i])) * 0.5
        dz = z - zs
        rs = fm.sqrt(dz * dz + r * r)
        r3 = rs * rs * rs
        usum = usum + div(float(src[i]) * dz, r3)
        vsum = vsum + div(float(src[i]) * r, r3)
    return div(usum, 4.0 * PI_EXACT), div(vsum, 4.0 * PI_EXACT)


def body_velocities(s: State, zprop: float, awake: float, z, area,
                    fix_source_positions: bool = False) -> None:
    """NACELL from the non-file data on: dimensional body stations (radius if area < 0),
    rotor position zprop and far-wake area awake (radius if < 0). Sets UBODY, VBODY, XW0
    (then XWINIT) and CONV = .FALSE."""
    if s.rad == 0.0:
        raise ValueError("NACELL: rotor geometry must be defined first")
    z = [float(v) for v in z]
    a = [float(v) for v in area]
    if not 1 <= len(z) == len(a) <= IZX:
        raise ValueError(f"need 1 to {IZX} body stations")
    awake = float(awake)
    if awake < 0.0:
        awake = PI * (awake * awake)
    a = [PI * (v * v) if v < 0.0 else v for v in a]
    _finish(s, float(zprop), awake, z, a, fix_source_positions)


def _finish(s, zprop, awake, z, a, fix):
    rad = s.rad
    zprop = div(zprop, rad)
    awake = div(awake, rad * rad)
    nz = len(z)
    zb = [div(v, rad) - zprop for v in z]
    ab = [div(v, rad * rad) for v in a]
    src = [0.0] * nz
    for i in range(1, nz):
        src[i] = ab[i] - ab[i - 1]                    # 7.55: placed one interval late
    if fix:
        for i in range(nz - 1):
            src[i] = ab[i + 1] - ab[i]
    for i in range(s.ii):
        s.ubody[i], s.vbody[i] = srcvel(zb, src, 0.0, float(s.xi[i]))
    _wake(s, awake)


def _wake(s, awake):
    s.xw0 = fm.sqrt(div(awake, PI))                  # label 20
    xwinit(s)
    s.conv = False


def constant_ubody(s: State, ubody: float) -> None:
    """NACE CONST: UBODY = ubody, VBODY = 0 at every station; 7.55 also sets XW0 = 0."""
    if s.rad == 0.0:
        raise ValueError("NACELL: rotor geometry must be defined first")
    s.greek = False
    s.ubody[:s.ii] = float(ubody)
    s.vbody[:s.ii] = 0.0
    _wake(s, 0.0)                                    # 7.55: AWAKE = 0, so XW0 = 0


def nacelle(s: State, path: str | os.PathLike, ubody: float = 0.0,
            fix_source_positions: bool = False) -> bool:
    """NACE f: read a body area file (or CONST..., with constant UBODY = ubody).
    Returns False where 7.55 reports a read error (it then sets UBODY = 0)."""
    s.greek = False
    if s.rad == 0.0:
        raise ValueError("NACELL: rotor geometry must be defined first")
    path = os.fspath(path)
    s.fname = path[:80]
    if path[:5].upper() == "CONST":                  # FNTEMP(1:5).EQ.'CONST' (the name as given)
        constant_ubody(s, ubody)
        return True
    try:
        with open(path) as f:
            rd = ListReader(f.read())
    except OSError:
        s.ubody[:s.ii] = 0.0                          # label 80: OPEN error
        return False
    try:
        zprop, awake = rd.read("rr")                  # 7.55: no END= (an empty file stops)
    except ValueError:
        s.ubody[:s.ii] = 0.0
        return False
    if zprop is None or awake is None:
        raise ValueError("NACELL: null value in the first record")
    if awake < 0.0:
        awake = PI * (awake * awake)
    z, a = [], []
    for _ in range(IZX):
        try:
            zk, ak = rd.read("rr")
        except EOFError:
            break
        except ValueError:
            s.ubody[:s.ii] = 0.0                      # label 90: READ error
            return False
        if zk is None or ak is None:
            raise ValueError("NACELL: null values (7.55 would keep stale values)")
        if ak < 0.0:
            ak = PI * (ak * ak)
        z.append(zk)
        a.append(ak)
    _finish(s, zprop, awake, z, a, fix_source_positions)
    return True
