"""Python access to the compiled XROTOR 7.55 reference library.

    ref = Ref755()                  # a fresh, independent copy of XROTOR's state
    ref["II"] = 40                  # read/write any COMMON variable by name
    ref.call("SETX")                # call any Fortran subroutine (arguments by reference)
    x = ref.get("XI")[:40]          # copies out; ref["XI"] is a live, writable view

Each Ref755 loads its own private copy of libxrotor755.so, so instances never share
COMMON-block state. Fortran STOP statements still end the whole Python process, so
reference calls that might STOP belong in a subprocess.
"""
from __future__ import annotations

import contextlib
import ctypes
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

from .commons import DTYPE, parse_include, parameters

# Unbuffered Fortran terminal output, so quiet() can silence it reliably. Must be set
# before libgfortran initializes, i.e. before the first reference library is loaded.
os.environ.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")

HERE = Path(__file__).resolve().parent
LIB = HERE / "build" / "libxrotor755.so"
INCLUDE = HERE / "xrotor-7.55" / "src" / "XROTOR.INC"


def available() -> bool:
    return LIB.exists() and INCLUDE.exists()


@contextlib.contextmanager
def quiet():
    """Send the Fortran library's terminal output (file descriptor 1) to /dev/null."""
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        yield
    finally:
        os.dup2(saved, 1)
        os.close(saved)
        os.close(devnull)


def _arg(x):
    """Convert a Python value into something Fortran can take by reference."""
    if isinstance(x, np.ndarray):
        if x.dtype == np.float64 or x.dtype == np.int32:
            return np.asfortranarray(x)
        raise TypeError(f"array arguments must be float64 or int32, not {x.dtype}")
    if isinstance(x, (bool, np.bool_)):
        return np.array(1 if x else 0, dtype=np.int32)      # gfortran .TRUE. == 1
    if isinstance(x, (int, np.integer)):
        return np.array(x, dtype=np.int32)
    if isinstance(x, (float, np.floating)):
        return np.array(x, dtype=np.float64)
    raise TypeError(f"unsupported argument type {type(x).__name__}")


class Ref755:
    """One isolated instance of XROTOR 7.55, initialized as the program does at start-up."""

    def __init__(self, lib: str | Path = LIB, start: bool = True, silent: bool = True):
        if not Path(lib).exists():
            raise FileNotFoundError(f"{lib} not built; run reference/build.sh")
        tmpdir = tempfile.mkdtemp(prefix="ref755-")
        private = Path(tmpdir) / "libxrotor755.so"
        shutil.copy2(lib, private)
        self._lib = ctypes.CDLL(str(private))
        shutil.rmtree(tmpdir)                 # the mapping outlives the file

        self.silent = silent
        self.params = parameters(INCLUDE)
        self.blocks = parse_include(INCLUDE)
        self._views: dict[str, np.ndarray] = {}
        self._where: dict[str, str] = {}
        for blk in self.blocks.values():
            base = ctypes.addressof(ctypes.c_char.in_dll(self._lib, blk.symbol))
            buf = (ctypes.c_char * blk.size).from_address(base)
            for v in blk.vars:
                dtype = f"S{v.itemsize}" if v.kind == "character" else DTYPE[v.kind]
                view = np.ndarray(v.shape, dtype=dtype, buffer=buf, offset=v.offset, order="F")
                self._views[v.name] = view
                self._where[v.name] = blk.name
        if start:
            self.call("REFSTART")

    # ------------------------------------------------------------------ state
    def __getitem__(self, name: str) -> np.ndarray:
        """Live view of a COMMON variable. Scalars are 0-d arrays; use .get() for values."""
        return self._views[name.upper()]

    def __setitem__(self, name: str, value) -> None:
        view = self._views[name.upper()]
        if view.dtype.kind == "S":
            view[...] = str(value).ljust(view.dtype.itemsize).encode()
        else:
            view[...] = value

    def get(self, name: str):
        """A copy of a COMMON variable: a Python scalar, string, or NumPy array."""
        view = self._views[name.upper()]
        if view.dtype.kind == "S":
            return view[()].decode().rstrip()
        return view.item() if view.ndim == 0 else view.copy()

    def snapshot(self, blocks=None) -> dict:
        """Copies of every variable in the given COMMON blocks (default: all)."""
        names = [n for n, b in self._where.items() if blocks is None or b in blocks]
        return {n: self.get(n) for n in names}

    def block_of(self, name: str) -> str:
        return self._where[name.upper()]

    # ------------------------------------------------------------------ calls
    def _run(self, name: str, args, restype=None):
        fn = getattr(self._lib, name.lower() + "_")
        fn.restype = restype
        conv = [_arg(a) for a in args]
        ptrs = [ctypes.c_void_p(c.ctypes.data) for c in conv]
        with (quiet() if self.silent else contextlib.nullcontext()):
            result = fn(*ptrs)
        return result, conv

    def call(self, name: str, *args) -> list[np.ndarray]:
        """Call a SUBROUTINE. Returns the arguments as arrays, so outputs can be read back."""
        return self._run(name, args)[1]

    def func(self, name: str, *args, kind: str = "real"):
        """Call a FUNCTION returning REAL (default), INTEGER, or LOGICAL."""
        restype = {"real": ctypes.c_double, "integer": ctypes.c_int32, "logical": ctypes.c_int32}[kind]
        value = self._run(name, args, restype)[0]
        return bool(value) if kind == "logical" else value
