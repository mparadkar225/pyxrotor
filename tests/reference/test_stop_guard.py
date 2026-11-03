"""A Fortran STOP in a reference call exits the whole process with status 0. In a single
process (-n0) that is before pytest writes its summary or JUnit report, so a truncated run
looks like a pass, and CI checks that the report exists. Under xdist only the worker dies,
and the controller reports the test as failed. Both behaviors are pinned here."""
import subprocess
import sys
from pathlib import Path

import pytest

from reference import ref755

HERE = Path(__file__).resolve().parent


pytestmark = pytest.mark.skipif(not ref755.available(), reason="XROTOR 7.55 reference not built")


def run_probe(report, workers):
    return subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                           "-n", workers, f"--junitxml={report}", str(HERE / "stop_probe.py")],
                          cwd=HERE.parents[1], capture_output=True, text=True,
                          stdin=subprocess.DEVNULL, timeout=300)


def test_fortran_stop_exits_zero_without_a_report(tmp_path):
    report = tmp_path / "junit.xml"
    p = run_probe(report, "0")
    assert p.returncode == 0
    assert "passed" not in p.stdout                      # no summary line
    assert not report.exists()


def test_fortran_stop_under_xdist_fails_the_test(tmp_path):
    report = tmp_path / "junit.xml"
    p = run_probe(report, "2")
    assert p.returncode == 1, (p.stdout[-800:], p.stderr[-800:])
    assert "test_fortran_stop" in p.stdout and "crashed" in p.stdout
    assert report.exists()
