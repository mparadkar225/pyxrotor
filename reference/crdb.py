"""The compiled CR-176883 data banks (reference/cr176883, REAL*8 build) for differential tests.

    lib = CRDB()
    lib.call("CLCD16", xm, alpha, ...)   # every argument by reference; returns them converted

Python float -> REAL*8, int -> INTEGER*4, bool -> LOGICAL. The routines keep no state
between calls except DATA tables, so one shared library instance is enough.
"""
from __future__ import annotations

import ctypes
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIB = HERE / "build" / "libcr176883.so"


def available() -> bool:
    return LIB.exists()


class CRDB:
    def __init__(self):
        self.lib = ctypes.CDLL(str(LIB))

    def call(self, name, *args):
        refs = []
        for a in args:
            if isinstance(a, bool):
                refs.append(ctypes.c_int32(1 if a else 0))
            elif isinstance(a, int):
                refs.append(ctypes.c_int32(a))
            else:
                refs.append(ctypes.c_double(float(a)))
        getattr(self.lib, name.lower() + "_")(*[ctypes.byref(r) for r in refs])
        return [r.value for r in refs]
