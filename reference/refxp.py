"""Python access to the compiled xrotor-python reference (DARcorporation, 2019).

The build patches its C interface from float32 to float64 and compiles it with the same
flags as the 7.55 reference, so the two can be compared at full double precision.

    XRotor, Case = load()
    xr = XRotor()                  # each instance loads a private copy of the library
    xr.case = Case.from_dict({...})
    xr.operate(rpm=2000)

xrotor-python has no DESI, MODI, or parametric section aero; sections are tabulated polars.
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

os.environ.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")

HERE = Path(__file__).resolve().parent
PKG_DIR = HERE / "build" / "xrotor-python"


def available() -> bool:
    return (PKG_DIR / "xrotor" / "libxrotor.so").exists()


def load():
    """Import the patched package. Returns (XRotor, Case)."""
    if not available():
        raise FileNotFoundError(f"{PKG_DIR} not built; run reference/build.sh")
    if str(PKG_DIR) not in sys.path:
        sys.path.insert(0, str(PKG_DIR))
    mod = importlib.import_module("xrotor")
    if not Path(mod.__file__).resolve().is_relative_to(PKG_DIR):
        raise ImportError(f"'xrotor' resolved to {mod.__file__}, not the patched reference build")
    return mod.XRotor, mod.Case


def test_case() -> dict:
    """The NASA Report 640 propeller case bundled with xrotor-python."""
    load()
    return importlib.import_module("xrotor.test").case
