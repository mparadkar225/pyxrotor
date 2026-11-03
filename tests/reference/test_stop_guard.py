"""A Fortran STOP in a reference call exits the whole process with status 0, before pytest
writes its summary or JUnit report: a truncated run looks like a pass. CI therefore checks
that the report exists. This pins the behavior that check relies on."""
import subprocess
import sys
from pathlib import Path

import pytest

from reference import ref755

HERE = Path(__file__).resolve().parent


@pytest.mark.skipif(not ref755.available(), reason="XROTOR 7.55 reference not built")
def test_fortran_stop_exits_zero_without_a_report(tmp_path):
    report = tmp_path / "junit.xml"
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                        f"--junitxml={report}", str(HERE / "stop_probe.py")],
                       cwd=HERE.parents[1], capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, timeout=300)
    assert p.returncode == 0
    assert "passed" not in p.stdout                      # no summary line
    assert not report.exists()
