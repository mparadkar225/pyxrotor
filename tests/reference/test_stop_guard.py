"""A Fortran STOP in a reference call exits the process with status 0, which would make a
truncated test run look like a pass. conftest's exit guard must turn it into a failure."""
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import STOP_EXIT

from reference import ref755

HERE = Path(__file__).resolve().parent


@pytest.mark.skipif(not ref755.available(), reason="XROTOR 7.55 reference not built")
def test_fortran_stop_fails_the_run():
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                        str(HERE / "stop_probe.py")], cwd=HERE.parents[1],
                       capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=300)
    assert p.returncode == STOP_EXIT, (p.returncode, p.stdout[-500:], p.stderr[-500:])
    assert "run is incomplete" in p.stderr
