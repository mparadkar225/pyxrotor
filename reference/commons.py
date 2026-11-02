"""Read the COMMON-block memory layout out of a Fortran 77 INCLUDE file.

XROTOR keeps all of its state in COMMON blocks declared in XROTOR.INC. Parsing that file
gives the offset, type, and shape of every variable, which lets the Python harness read
and write the compiled reference's state directly (see ref755.py).

Only the declaration forms XROTOR.INC actually uses are supported: IMPLICIT, PARAMETER,
DIMENSION, typed declarations (REAL, INTEGER, LOGICAL, CHARACTER*n) and COMMON.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# Byte sizes under the reference build flags (-fdefault-real-8, default INTEGER and LOGICAL).
ITEMSIZE = {"real": 8, "integer": 4, "logical": 4}
DTYPE = {"real": "<f8", "integer": "<i4", "logical": "<i4"}


@dataclass
class Var:
    name: str
    kind: str                      # real | integer | logical | character
    bounds: tuple = ()             # ((lo, hi), ...) per dimension; () for scalars
    itemsize: int = 8
    offset: int = 0                # bytes from the start of the block

    @property
    def shape(self) -> tuple:
        return tuple(hi - lo + 1 for lo, hi in self.bounds)

    @property
    def count(self) -> int:
        n = 1
        for s in self.shape:
            n *= s
        return n

    @property
    def nbytes(self) -> int:
        return self.count * self.itemsize


@dataclass
class Block:
    name: str
    vars: list = field(default_factory=list)
    size: int = 0

    @property
    def symbol(self) -> str:
        """gfortran's linker symbol for the block."""
        return self.name.lower() + "_"


def _statements(text: str) -> list[str]:
    """Join fixed-form continuation lines and drop comments."""
    out: list[str] = []
    for raw in text.splitlines():
        if not raw.strip() or raw[0] in "Cc*!":
            continue
        line = raw[:72]
        if "!" in line:
            line = line[: line.index("!")]
        if len(line) > 5 and line[5] not in " 0":
            out[-1] += line[6:]
        else:
            out.append(line[6:])
    return [s.strip() for s in out if s.strip()]


def _split_top(s: str) -> list[str]:
    """Split on commas that are not inside parentheses."""
    parts, depth, cur = [], 0, ""
    for ch in s:
        depth += ch == "("
        depth -= ch == ")"
        if ch == "," and depth == 0:
            parts.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        parts.append(cur.strip())
    return parts


def _eval(expr: str, params: dict) -> float | int:
    expr = expr.strip().upper()
    if "." not in expr:
        expr = expr.replace("/", "//")       # Fortran integer division
    return eval(expr, {"__builtins__": {}}, dict(params))


def _entity(text: str) -> tuple[str, str | None]:
    """'NAME(dims)' -> ('NAME', 'dims') ; 'NAME' -> ('NAME', None)."""
    m = re.fullmatch(r"([A-Z_][A-Z0-9_]*)\s*(?:\((.*)\))?", text.strip().upper())
    if not m:
        raise ValueError(f"cannot parse declaration entity {text!r}")
    return m.group(1), m.group(2)


def parse_include(path: str | Path) -> dict[str, Block]:
    params: dict = {}
    implicit = {c: ("integer" if "I" <= c <= "N" else "real") for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
    types: dict[str, tuple[str, int]] = {}
    dims: dict[str, str] = {}
    blocks: dict[str, Block] = {}

    for st in _statements(Path(path).read_text()):
        up = st.upper()
        if up.startswith("IMPLICIT"):
            m = re.fullmatch(r"IMPLICIT\s+(REAL|INTEGER|LOGICAL)\s*\((.*)\)", up)
            for item in m.group(2).split(","):
                a, _, b = item.strip().partition("-")
                for code in range(ord(a), ord(b or a) + 1):
                    implicit[chr(code)] = m.group(1).lower()
        elif up.startswith("PARAMETER"):
            for item in _split_top(up[up.index("(") + 1: up.rindex(")")]):
                k, v = item.split("=")
                params[k.strip()] = _eval(v, params)
        elif up.startswith("DIMENSION"):
            for ent in _split_top(up[len("DIMENSION"):]):
                n, d = _entity(ent)
                dims[n] = d
        elif re.match(r"(REAL|INTEGER|LOGICAL|CHARACTER)\b", up):
            m = re.match(r"(REAL|INTEGER|LOGICAL|CHARACTER)(\*\d+)?\s+(.*)", up)
            kind = m.group(1).lower()
            size = int(m.group(2)[1:]) if m.group(2) else ITEMSIZE.get(kind, 1)
            for ent in _split_top(m.group(3)):
                n, d = _entity(ent)
                types[n] = (kind, size)
                if d:
                    dims[n] = d
        elif up.startswith("COMMON"):
            for bname, body in re.findall(r"/\s*([A-Z0-9_]+)\s*/([^/]*)", up):
                blk = Block(bname)
                off = 0
                for ent in _split_top(body):
                    n, d = _entity(ent)
                    d = d or dims.get(n)
                    kind, size = types.get(n, (implicit[n[0]], None))
                    size = size or ITEMSIZE[kind]
                    bounds = ()
                    if d:
                        bl = []
                        for dd in _split_top(d):
                            lo, _, hi = dd.partition(":")
                            bl.append((_eval(lo, params), _eval(hi, params)) if hi else (1, _eval(lo, params)))
                        bounds = tuple(bl)
                    align = 1 if kind == "character" else size
                    off += (-off) % align
                    v = Var(n, kind, bounds, size, off)
                    blk.vars.append(v)
                    off += v.nbytes
                blk.size = off
                blocks[bname] = blk
    return blocks


def parameters(path: str | Path) -> dict:
    """The PARAMETER constants of an include file (IX, NAX, PI, ...)."""
    params: dict = {}
    for st in _statements(Path(path).read_text()):
        up = st.upper()
        if up.startswith("PARAMETER"):
            for item in _split_top(up[up.index("(") + 1: up.rindex(")")]):
                k, v = item.split("=")
                params[k.strip()] = _eval(v, params)
    return params
