import pytest

from reference import ref755, refxp


def pytest_configure(config):
    config.addinivalue_line("markers", "reference: checks the compiled reference libraries themselves")
    config.addinivalue_line("markers", "differential: compares a ported routine with the 7.55 reference")
    config.addinivalue_line("markers", "measurement: compares against wind-tunnel or flight data")


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
