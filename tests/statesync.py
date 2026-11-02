"""Moving state between pyxrotor.State and the Fortran reference, and comparing them."""
from __future__ import annotations

import math

import numpy as np

from harness import assert_close
from pyxrotor.state import INDEX_VARS, State

# Blocks that hold scratch space, plotting settings, or I/O unit numbers. Ported routines
# use locals instead of TEMP_X, so these are never compared.
UNCOMPARED_BLOCKS = {"TEMP_X", "IPLO_X", "RPLO_X", "UNIT_X", "CHAR_X"}


def to_python(name: str, fortran_value):
    kind, bounds = State.SPEC[name]
    if kind.startswith("character"):
        return fortran_value
    if bounds:
        arr = np.array(fortran_value)
        if kind == "logical":
            return arr != 0
        if kind == "integer":
            arr = arr.astype(np.int64)
            return arr - 1 if name in INDEX_VARS else arr
        return arr.astype(np.float64)
    if kind == "logical":
        return bool(fortran_value)
    if kind == "integer":
        v = int(fortran_value)
        return v - 1 if name in INDEX_VARS else v
    return float(fortran_value)


def to_fortran(name: str, py_value):
    kind, bounds = State.SPEC[name]
    if kind.startswith("character"):
        return py_value
    if bounds:
        arr = np.asarray(py_value)
        if kind == "logical":
            return arr.astype(np.int32)
        if kind == "integer":
            return (arr + 1 if name in INDEX_VARS else arr).astype(np.int32)
        return arr
    if kind == "logical":
        return int(bool(py_value))
    if kind == "integer":
        return int(py_value) + (1 if name in INDEX_VARS else 0)
    return float(py_value)


def state_from_ref(ref) -> State:
    """A State holding a copy of everything in the reference's COMMON blocks."""
    s = State()
    for name in State.SPEC:
        setattr(s, name.lower(), to_python(name, ref.get(name)))
    return s


def state_to_ref(s: State, ref, names=None) -> None:
    """Write State variables (default: all) into the reference's COMMON blocks."""
    for name in names or State.SPEC:
        ref[name] = to_fortran(name, getattr(s, name.lower()))


def compared_names(names=None):
    return [n for n in (names or State.SPEC) if State.BLOCK_OF[n] not in UNCOMPARED_BLOCKS]


def _same_array(name, got, want, rtol):
    got, want = np.asarray(got), np.asarray(want)
    if got.dtype.kind in "biu" or want.dtype.kind in "biu":
        if not np.array_equal(got, want):
            k = np.argwhere(got != want)[0]
            raise AssertionError(f"{name}: first difference at {tuple(k)}: got {got[tuple(k)]}, want {want[tuple(k)]}")
        return
    nan_g, nan_w = np.isnan(got), np.isnan(want)
    if not np.array_equal(nan_g, nan_w):
        k = np.argwhere(nan_g != nan_w)[0]
        raise AssertionError(f"{name}: NaN mismatch at {tuple(k)}: got {got[tuple(k)]}, want {want[tuple(k)]}")
    assert_close(np.where(nan_g, 0.0, got), np.where(nan_w, 0.0, want), rtol=rtol, name=name)


def assert_states_match(s: State, ref, names=None, rtol: float = 0.0) -> None:
    """Every compared variable of s equals the reference's (NaN == NaN)."""
    for name in compared_names(names):
        got = getattr(s, name.lower())
        want = to_python(name, ref.get(name))
        if isinstance(got, np.ndarray):
            _same_array(name, got, want, rtol)
        elif isinstance(got, float):
            if math.isnan(want):
                assert math.isnan(got), f"{name}: got {got!r}, want nan"
            else:
                assert_close(got, want, rtol=rtol, name=name)
        else:
            assert got == want, f"{name}: got {got!r}, want {want!r}"
