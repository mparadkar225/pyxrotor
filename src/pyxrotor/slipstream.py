"""External slipstream profiles: SUBROUTINE SAVVEL (VPUT) and SUBROUTINE GETVEL (VGET),
xrotor.f. This is how XROTOR couples rotors: VPUT writes a solved rotor's slipstream, VGET
imposes such a profile on another rotor's inflow (RADD, UADD, VADD, used by UVADD).

VPUT writes one record per blade station, list-directed:
    r (m)    Vaxial (m/s)    Vrot (m/s)
the circumferentially averaged induced velocities *far downstream* (twice the disk
values, from the helical-wake relation Va = Vt*r/ADW), at the *disk* radii. The rotor's
own imposed inflow and nacelle velocities are not included. Vrot is positive in the
direction of rotation, so a counter-rotating rotor downstream imports it with weight -1.

VGET reads at least 2 such records (at most IX), multiplies the axial and tangential
velocities by the given weights, and splines them in r (UADDR, VADDR). UVADD imposes
nothing on stations outside the profile's radius range.

7.55 behavior reproduced (see CLAUDE.md):
  * The radii are the disk radii, not the contracted far-wake radii the velocities belong to.
  * LOAD reads RADD/UADD/VADD from a restart file but never splines them, so UVADD uses
    stale derivatives (zeros in a fresh session). Fix: State.fixes.slipstream_splines,
    applied in restart.load.
"""
from __future__ import annotations

import os

import numpy as np

from .constants import IX, PI
from .fmath import div
from .fortio import ListReader, write_list
from .numerics import spline
from .state import State
from .wake import uvadd


def profile(s: State) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The (r, Vaxial, Vrot) arrays SAVVEL writes, from the current solution."""
    if s.rad == 0.0 or not s.conv:
        raise ValueError("SAVVEL: no converged operating point ('Must define operating point first.')")
    blds = float(s.nblds)
    ii = s.ii
    r = np.empty(ii)
    u = np.empty(ii)
    v = np.empty(ii)
    for i in range(ii):
        xi = float(s.xi[i])
        gam = float(s.gam[i])
        vt = div(blds * gam, (4.0 * PI) * xi)       # circumferentially averaged
        va = div(vt * xi, s.adw)
        uduct, vaduct_va = 0.0, 1.0
        if s.duct:
            uduct = s.urduct - 1.0
            vaduct_va = 2.0 * s.urduct
        vd = va * (vaduct_va - 1.0)
        utot = 1.0 + uduct + float(s.ubody[i])
        wa, wt = uvadd(s, xi)
        ci = div(xi, s.adv) - wt - vt
        si = utot + wa + va + vd
        r[i] = xi * s.rad
        u[i] = (2.0 * (si - (utot + wa))) * s.vel
        v[i] = (-2.0 * (ci - (div(xi, s.adv) - wt))) * s.vel
    return r, u, v


def savvel(s: State, path: str | os.PathLike, overwrite: bool = True) -> None:
    """VPUT f. XROTOR asks before overwriting (default yes); overwrite=False raises."""
    path = os.fspath(path)
    s.greek = False
    r, u, v = profile(s)
    if os.path.exists(path) and not overwrite:
        raise FileExistsError(path)
    s.fname = path[:80]
    with open(path, "w") as f:
        for k in range(len(r)):
            f.write(write_list(r[k], u[k], v[k]) + "\n")


def set_profile(s: State, r, u, v, axial_weight: float, swirl_weight: float) -> None:
    """GETVEL from arrays: impose (r m, Vaxial m/s, Vrot m/s), weighted, and spline it."""
    r = np.asarray(r, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    n = r.size
    if not 2 <= n <= IX or u.size != n or v.size != n:
        raise ValueError(f"need 2 to {IX} points with r, u, v of equal length")
    s.nadd = n
    s.radd[:n], s.uadd[:n], s.vadd[:n] = r, u, v
    uwt, vwt = float(axial_weight), float(swirl_weight)
    s.uadd[:n] = uwt * s.uadd[:n]
    s.vadd[:n] = vwt * s.vadd[:n]
    respline(s)


def respline(s: State) -> None:
    """CALL SPLINE(UADD,UADDR,RADD,NADD) and the same for VADD."""
    n = s.nadd
    if n >= 2:
        s.uaddr[:n] = spline(s.uadd[:n], s.radd[:n], n)
        s.vaddr[:n] = spline(s.vadd[:n], s.radd[:n], n)


def getvel(s: State, path: str | os.PathLike, axial_weight: float, swirl_weight: float) -> bool:
    """VGET f: read a VPUT file and impose it with the given weights (XROTOR's prompts:
    axial 0 to 1, tangential 0 or +-1). Returns False where XROTOR prints 'New wake
    velocities not read' (a read error, or fewer than 2 records); NADD is then unchanged."""
    path = os.fspath(path)
    s.greek = False
    s.fname = path[:80]
    try:
        with open(path) as f:
            rd = ListReader(f.read())
    except OSError:
        return False
    rows = []
    for _ in range(IX):
        try:
            vals = rd.read("rrr")
        except EOFError:
            break
        except ValueError:
            return False                            # 'File read error'
        if any(x is None for x in vals):
            raise ValueError("VGET: null values (7.55 would keep stale values)")
        rows.append(vals)
    if len(rows) < 2:
        return False
    a = np.array(rows, dtype=float)
    set_profile(s, a[:, 0], a[:, 1], a[:, 2], axial_weight, swirl_weight)
    return True
