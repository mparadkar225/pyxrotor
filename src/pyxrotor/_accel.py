"""Switch between the compiled hot loops (pyxrotor._speedups, built from _speedups.pyx by
setup.py) and the pure-Python implementations, which are the reference and always present.

Both give bit-identical results. The compiled ones are used when the extension was built,
unless the environment variable PYXROTOR_PURE=1 is set, or code calls use_compiled(False):

    from pyxrotor import _accel
    with _accel.pure():                 # run something on the pure-Python code
        ...
"""
from __future__ import annotations

import contextlib
import os

try:
    from . import _speedups as ext
except ImportError:                     # not built: no compiler, PYXROTOR_PURE=1 at install, ...
    ext = None

available = ext is not None
enabled = available and os.environ.get("PYXROTOR_PURE", "") not in ("1", "true", "yes")


def use_compiled(on: bool) -> None:
    """Turn the compiled loops on (if built) or off."""
    global enabled
    enabled = bool(on) and available


@contextlib.contextmanager
def pure():
    """Run the enclosed code on the pure-Python implementations."""
    global enabled
    saved = enabled
    enabled = False
    try:
        yield
    finally:
        enabled = saved
