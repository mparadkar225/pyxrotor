"""XROTOR restart files: SUBROUTINE SAVE and SUBROUTINE LOAD (xio.f), the 7.55 format
("XROTOR VERSION: 7.55" header, '!' comment lines, aero sections included).

    save_lines(s)  -> the file's records, byte for byte as 7.55 writes them
    save(s, path)
    load(s, path, stations=None)

Not ported: SAVEOLD/LOADOLD (the pre-6.9 format; load raises NotImplementedError for it).

7.55 behavior reproduced (see CLAUDE.md):
  * SAVE writes the flags LVDUCT, DUCT, WIND but LOAD reads that record as FREE, DUCT, WIND,
    so after LOAD the free-wake flag is "ADW differs from ADV*URDUCT by 5e-5 or more".
    (fix: load(..., keep_free_wake=True) keeps the session's FREE instead).
  * SAVE writes the blade angles BETA0, not BETA.
  * RDLINE treats a line reading exactly END or ERR as the end of the file.
  * LOAD reads XI(1..IIX) from the file into COMMON before SETX regenerates XI(1..II), so
    with IIX > II stations II+1..IIX keep the file's radii.
Malformed files: 7.55 either reports 'incompatible format' or stops with a runtime error
(READs without END=/ERR=); here both raise LoadError, and the State is not guaranteed to
match what 7.55 leaves behind on those paths.
"""
from __future__ import annotations

import os

from .aero import setiaero
from .constants import IX, PI
from .design import _set_grid
from .fmath import div
from .fortio import fmt_f, fmt_g, fmt_i, fmt_l, read_list
from .numerics import seval, spline
from .oper import aper, sections
from .rotor import put_aero
from .aero import AeroSection
from .state import State
from .wake import xwinit


class LoadError(ValueError):
    """LOAD's 'File ... has incompatible format' (or a Fortran runtime READ error)."""


def _g(*vals) -> str:                                   # FORMAT(5(1X,G12.5))
    return "".join(" " + fmt_g(v, 12, 5) for v in vals)


def save_lines(s: State) -> list[str]:
    """The records SUBROUTINE SAVE writes. Also sets GREEK = .FALSE. and, if the case name
    is blank, NAME = 'Saved blade', as SAVE does."""
    s.greek = False
    lvduct = abs(s.adw - s.adv * s.urduct) >= 5.0e-5
    if s.name.strip() == "":
        s.name = "Saved blade"
    out = ["XROTOR VERSION: " + fmt_f(s.version, 5, 2), s.name.ljust(32)[:32]]
    out += ["!         Rho          Vso          Rmu           Alt", _g(s.rho, s.vso, s.rmu, s.alt)]
    out += ["!         Rad          Vel          Adv          Rake", _g(s.rad, s.vel, s.adv, s.rake)]
    out += ["!         XI0          XIW", _g(s.xi0, s.xw0)]
    out += ["!  Naero", " " + fmt_i(s.naero, 5)]
    for sec in sections(s):
        a0deg = (sec.a0 * 180.0) / PI
        out += ["!   Xisection", _g(sec.xisect)]
        out += ["!       A0deg        dCLdA        CLmax         CLmin",
                _g(a0deg, sec.dclda, sec.clmax, sec.clmin)]
        out += ["!  dCLdAstall     dCLstall      Cmconst         Mcrit",
                _g(sec.dclda_stall, sec.dcl_stall, sec.cmcon, sec.mcrit)]
        out += ["!       CDmin      CLCDmin     dCDdCL^2", _g(sec.cdmin, sec.cldmin, sec.dcdcl2)]
        out += ["!       REref        REexp", _g(sec.reref, sec.rexp)]
    out += ["!LVDuct  LDuct   LWind",
            "".join(" " + fmt_l(v, 2) + "     " for v in (lvduct, s.duct, s.wind)).rstrip()]
    out += ["!   II Nblds", " " + fmt_i(s.ii, 5) + " " + fmt_i(s.nblds, 5),
            "!         r/R          C/R     Beta0deg         Ubody"]
    for i in range(s.ii):
        out.append(_g(s.xi[i], s.ch[i], (float(s.beta0[i]) * 180.0) / PI, s.ubody[i]))
    out += ["!      URDuct", _g(s.urduct)]
    if s.nadd > 1:
        out += ["!Nadd", " " + fmt_i(s.nadd, 5), "!        Radd         Uadd         Vadd"]
        for i in range(s.nadd):
            out.append(_g(s.radd[i], s.uadd[i], s.vadd[i]))
    return out


def save(s: State, path: str | os.PathLike, overwrite: bool = True) -> None:
    """SAVE f. XROTOR asks before overwriting (default yes); overwrite=False raises
    FileExistsError instead."""
    path = os.fspath(path)
    if os.path.exists(path) and not overwrite:
        raise FileExistsError(path)
    lines = save_lines(s)
    s.fname = path[:80]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


class _Reader:
    """RDLINE over the file's records: skips '!'/'#' comments and blank lines, 128-char LINE."""

    def __init__(self, text: str):
        self.recs = text.split("\n")
        if self.recs and self.recs[-1] == "":
            self.recs.pop()
        self.k = 0

    def line(self) -> str:
        while self.k < len(self.recs):
            rec = self.recs[self.k].rstrip("\r")[:128]
            self.k += 1
            if rec[:1] in ("!", "#"):
                continue
            if rec.strip(" ") == "":                    # Fortran blanks only, not tabs
                continue
            return rec
        return "END"


def _is_end(line: str) -> bool:
    return line.rstrip() in ("END", "ERR")


def _read(line: str, kinds: str) -> list:
    try:
        vals = read_list(line, kinds)
    except (ValueError, EOFError) as e:
        raise LoadError(f"cannot read {line!r}: {e}") from e
    if any(v is None for v in vals):
        raise LoadError(f"null values in {line!r} (7.55 would keep stale values)")
    return vals


def load(s: State, path: str | os.PathLike, stations: int | None = None,
         keep_free_wake: bool = False) -> bool:
    """LOAD f: read a rotor and its aero sections, re-spline the blade onto II stations
    (stations: the answer to 7.55's prompt when the file's count differs; None keeps II),
    solve the operating point at the file's advance ratio (APER 4,2) and set the DESI
    parameters from it. Returns CONV."""
    path = os.fspath(path)
    s.greek = False
    s.fname = path[:80]
    with open(path) as f:
        rd = _Reader(f.read())
    line = rd.line()
    if _is_end(line):
        raise LoadError(f"{path}: empty file")
    if line[:6] != "XROTOR":
        raise NotImplementedError("old-format (pre-6.9) XROTOR file: LOADOLD is not ported")
    try:
        float(line.ljust(22)[16:22].replace("D", "E"))      # FILEVERS (7.55 stops if bad)
    except ValueError as e:
        raise LoadError(f"bad version field in {line!r}") from e

    s.name = rd.line()[:32].rstrip()
    s.rho, s.vso, s.rmu, s.alt = _read(rd.line(), "rrrr")
    s.rad, s.vel, s.adv, s.rake = _read(rd.line(), "rrrr")
    s.xi0, s.xw0 = _read(rd.line(), "rr")
    (naero,) = _read(rd.line(), "i")
    if not 0 <= naero <= s.xiaero.size:
        raise LoadError(f"NAERO = {naero} out of range")
    s.naero = naero
    for n in range(naero):
        (xisect,) = _read(rd.line(), "r")
        a0deg, dclda, clmax, clmin = _read(rd.line(), "rrrr")
        dclda_stall, dcl_stall, cmcon, mcrit = _read(rd.line(), "rrrr")
        cdmin, cldmin, dcdcl2 = _read(rd.line(), "rrr")
        reref, rexp = _read(rd.line(), "rr")
        put_aero(s, n, AeroSection(xisect=xisect, a0=(a0deg * PI) / 180.0, clmax=clmax,
                                   clmin=clmin, dclda=dclda, dclda_stall=dclda_stall,
                                   dcl_stall=dcl_stall, cdmin=cdmin, cldmin=cldmin,
                                   dcdcl2=dcdcl2, cmcon=cmcon, mcrit=mcrit, reref=reref, rexp=rexp))

    free, s.duct, s.wind = _read(rd.line(), "lll")         # 7.55: SAVE wrote LVDUCT here
    if not keep_free_wake:
        s.free = free
    line = rd.line()
    if _is_end(line):
        raise LoadError("missing blade stations")
    iix, s.nblds = _read(line, "ii")
    if not 1 <= iix <= IX:
        raise LoadError(f"{iix} stations in file (1 to {IX} allowed)")
    for i in range(iix):
        x, c, bdeg, ub = _read(rd.line(), "rrrr")
        s.xi[i], s.ch[i], s.ubody[i] = x, c, ub
        s.beta[i] = (bdeg * PI) / 180.0
        s.beta0[i] = s.beta[i]

    # ---- optional duct velocity and slipstream profiles
    s.urduct = 1.0
    line = rd.line()
    if not _is_end(line):
        try:
            (s.urduct,) = read_list(line, "r", current=[s.urduct])
        except EOFError:
            pass
        except ValueError as e:
            raise LoadError(str(e)) from e
    s.nadd = 0
    line = rd.line()
    if not _is_end(line):
        try:
            (nadd,) = read_list(line, "i", current=[0])
        except EOFError:
            nadd = None
        except ValueError as e:
            raise LoadError(str(e)) from e
        if nadd is not None:
            s.nadd = min(nadd, IX)
            for i in range(s.nadd):
                line = rd.line()
                ok = not _is_end(line)
                if ok:
                    try:
                        s.radd[i], s.uadd[i], s.vadd[i] = read_list(line, "rrr")
                    except (ValueError, EOFError):
                        ok = False
                if not ok:                                 # label 20: I = i+1
                    if i + 1 > 2:
                        s.nadd = i
                    break
    s.conv = False

    if iix != s.ii and stations is not None:
        s.ii = int(stations)
        if not 1 <= s.ii <= IX:
            raise ValueError(f"stations must be 1 to {IX}")

    # ---- spline the file's blade onto the generated stations
    w1 = s.xi[:iix].copy()
    w2, w4, w6 = s.ch[:iix].copy(), s.beta[:iix].copy(), s.ubody[:iix].copy()
    w3, w5, w7 = spline(w2, w1, iix), spline(w4, w1, iix), spline(w6, w1, iix)
    _set_grid(s)
    xwinit(s)
    ii = s.ii
    for i in range(ii):
        x = float(s.xi[i])
        s.ch[i] = seval(x, w2, w3, w1, iix)
        s.beta[i] = seval(x, w4, w5, w1, iix)
        s.ubody[i] = seval(x, w6, w7, w1, iix)
        s.beta0[i] = s.beta[i]
    s.iinf = ii + ii // 2

    s.iaero[:ii] = setiaero(s.xi[:ii], sections(s))
    aper(s, 4, 2, True)

    # ---- design quantities for an MIL design with the same parameters
    v, r = s.vel, s.rad
    s.raddes, s.veldes, s.advdes = r, v, 0.0
    s.rpmdes = (div(v, r * s.adv) * 30.0) / PI
    s.r0des, s.rwdes = s.xi0 * r, s.xw0 * r
    s.tddes = ((s.ttot * s.rho) * (v * v)) * (r * r)
    s.pddes = ((s.ptot * s.rho) * (v * v * v)) * (r * r)
    s.dest, s.desp = False, True
    s.cldes[:ii] = s.cl[:ii]
    s.cldes0 = 0.0
    s.lrotor = True
    return bool(s.conv)
