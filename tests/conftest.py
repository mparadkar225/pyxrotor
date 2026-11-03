import os

import pytest
from hypothesis import settings

from reference import ref755, refxp

# HYPOTHESIS_PROFILE=deep for long bug hunts; the default keeps CI under a minute.
settings.register_profile("default", max_examples=400)
settings.register_profile("deep", max_examples=20000)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))

# For tests that run whole solves: 1/8 of the profile's examples (50 by default, 2500 deep).
EXPENSIVE = settings(deadline=None, max_examples=max(25, settings().max_examples // 8))


# A Fortran STOP inside a reference call (e.g. a shim's input check) calls the C library's
# exit(0): pytest dies mid-run, its buffered output is lost, and the exit status is 0, so a
# truncated run looks like a pass. A C-level exit handler turns that into a failure. It
# runs on every exit, so it only fails if pytest did not get to pytest_unconfigure.
_SESSION = {"finished": False}
STOP_EXIT = 13                                 # distinct from pytest's own exit codes


def _install_stop_guard():
    import ctypes
    import sys

    @ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_void_p)
    def guard(status, arg):
        if not _SESSION["finished"]:
            capman = _SESSION.get("capman")
            if capman is not None:                 # pytest is capturing fd 2: let it go
                try:
                    capman.suspend_global_capture(in_=False)
                except Exception:
                    pass
            sys.stderr.write("\nERROR: the process exited inside a test (a Fortran STOP in a "
                             "reference call?). The run is incomplete.\n")
            sys.stderr.flush()
            os._exit(STOP_EXIT)

    libc = ctypes.CDLL(None)
    if hasattr(libc, "on_exit"):               # glibc; elsewhere the guard is not installed
        _SESSION["guard"] = guard              # keep the callback alive
        libc.on_exit(guard, None)


def pytest_unconfigure(config):
    _SESSION["finished"] = True


def pytest_configure(config):
    _SESSION["capman"] = config.pluginmanager.getplugin("capturemanager")
    _install_stop_guard()
    config.addinivalue_line("markers", "reference: checks the compiled reference libraries themselves")
    config.addinivalue_line("markers", "differential: compares a ported routine with the 7.55 reference")
    config.addinivalue_line("markers", "measurement: compares against wind-tunnel or flight data")
    config.addinivalue_line("markers", "model: checks new (non-XROTOR) models against exact limits and independent methods")


@pytest.fixture
def ref():
    """A fresh XROTOR 7.55 instance in its start-up state."""
    if not ref755.available():
        pytest.skip("XROTOR 7.55 reference not built (run reference/build.sh)")
    return ref755.Ref755()


@pytest.fixture(scope="module")
def ref_module():
    """One XROTOR 7.55 instance shared by a module, for routines whose results depend only on
    the inputs a test sets explicitly (much faster under Hypothesis)."""
    if not ref755.available():
        pytest.skip("XROTOR 7.55 reference not built (run reference/build.sh)")
    return ref755.Ref755()


@pytest.fixture
def xrotor_python():
    if not refxp.available():
        pytest.skip("xrotor-python reference not built (run reference/build.sh)")
    return refxp.load()
