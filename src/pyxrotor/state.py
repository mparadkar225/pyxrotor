"""XROTOR's complete state as one object.

XROTOR 7.55 keeps everything in 27 COMMON blocks (XROTOR.INC). State carries exactly those
variables, so ported routines read and write the same quantities as the Fortran and a test
can compare the whole state after any routine.

Conventions
    * Attribute names are the Fortran names in lower case: XI -> s.xi, VW_GAM -> s.vw_gam.
    * Arrays keep Fortran's shapes and are indexed from 0: Fortran A(I,J) is a[I-1, J-1].
      For an array declared A(0:N), a[k] is A(k).
    * Variables that hold 1-based indices into other arrays (INDEX_VARS) store 0-based
      values, with Fortran's "none" (0) stored as -1.
    * LOGICAL -> bool, INTEGER -> int, REAL -> float, CHARACTER -> str.
"""
from __future__ import annotations

import copy

import numpy as np

from ._layout import BLOCKS

# Variables whose values are 1-based indices in Fortran.
INDEX_VARS = frozenset({"IAERO"})

_DTYPE = {"real": np.float64, "integer": np.int64, "logical": np.bool_}
_SCALAR = {"real": 0.0, "integer": 0, "logical": False}


def _shape(bounds):
    return tuple(hi - lo + 1 for lo, hi in bounds)


class State:
    """All XROTOR COMMON variables, zero-initialized (Fortran COMMON starts zeroed)."""

    #: Fortran variable name -> COMMON block name
    BLOCK_OF = {name: blk for blk, vs in BLOCKS.items() for name, _, _ in vs}
    #: Fortran variable name -> (kind, fortran bounds)
    SPEC = {name: (kind, bounds) for vs in BLOCKS.values() for name, kind, bounds in vs}

    def __init__(self):
        for name, (kind, bounds) in self.SPEC.items():
            if kind.startswith("character"):
                value = ""
            elif bounds:
                value = np.zeros(_shape(bounds), dtype=_DTYPE[kind])
            else:
                value = _SCALAR[kind]
            setattr(self, name.lower(), value)
        if "IAERO" in self.SPEC:
            self.iaero[:] = -1                      # Fortran 0 = unset

    def copy(self) -> "State":
        return copy.deepcopy(self)

    @classmethod
    def names(cls, blocks=None) -> list[str]:
        """Fortran names of all variables, or of those in the given blocks."""
        return [n for n, b in cls.BLOCK_OF.items() if blocks is None or b in blocks]

    def __repr__(self):
        return f"<State II={self.ii} NBLDS={self.nblds} ADV={self.adv:.5g} CONV={self.conv}>"
