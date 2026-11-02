"""Comparison helpers for differential tests against the Fortran references."""
from __future__ import annotations

import numpy as np

# Default tolerance for a port routine against the 7.55 reference. Routines that are
# straight-line arithmetic should match exactly (rtol=0); iterative solvers may need a
# little room, but anything looser than this needs a comment explaining why.
RTOL = 1e-12


def max_rel_error(got, want, floor: float = 1e-300) -> float:
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    if got.shape != want.shape:
        raise AssertionError(f"shape mismatch: got {got.shape}, want {want.shape}")
    same = (got == want) | (np.isnan(got) & np.isnan(want))
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        rel = np.where(same, 0.0, np.abs(got - want) / np.maximum(np.abs(want), floor))
    return float(np.max(np.where(np.isnan(rel), np.inf, rel), initial=0.0))


def assert_close(got, want, rtol: float = RTOL, atol: float = 0.0, name: str = "value") -> None:
    """Assert |got - want| <= atol + rtol*|want| elementwise, reporting the worst element.
    Identical values (including equal infinities) and NaN vs NaN count as equal."""
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    if got.shape != want.shape:
        raise AssertionError(f"{name}: shape mismatch: got {got.shape}, want {want.shape}")
    same = (got == want) | (np.isnan(got) & np.isnan(want))
    with np.errstate(invalid="ignore", over="ignore"):
        err = np.where(same, 0.0, np.abs(got - want))
        bad = ~same & ~(err <= atol + rtol * np.abs(want))
    if np.any(bad):
        k = np.unravel_index(np.argmax(np.where(bad, err, -1)), err.shape) if err.ndim else ()
        raise AssertionError(
            f"{name}: {int(np.sum(bad))} of {err.size} elements outside rtol={rtol:g} atol={atol:g}; "
            f"worst at {k}: got {got[k]!r}, want {want[k]!r} (rel err {max_rel_error(got, want):.3g})"
        )


def assert_state_close(got: dict, want: dict, names, rtol: float = RTOL, atol: float = 0.0) -> None:
    """Compare selected variables from two state snapshots (e.g. Ref755.snapshot())."""
    for n in names:
        assert_close(got[n], want[n], rtol=rtol, atol=atol, name=n)
