"""Run only by test_stop_guard.py (its name keeps it out of normal collection): a reference
call that hits a Fortran STOP, which exits the process."""
import numpy as np

from reference.ref755 import Ref755


def test_fortran_stop():
    ref = Ref755()
    r = np.array([0.3, 0.2, 0.9])                # not increasing: REFARBI stops
    ref.call("REFARBI", 2, 50.0, 1.0, 0.1, 3, r, np.full(3, 0.1), np.full(3, 30.0))


def test_after_stop():
    pass
