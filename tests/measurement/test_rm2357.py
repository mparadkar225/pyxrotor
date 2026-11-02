"""ARC R&M 2357: the Rotol wooden Spitfire IX propeller in the RAE 24-ft tunnel.

The data checks run now. The comparison with pyxrotor runs once OPER is ported; its
tolerance is an engineering one, because the section polars are estimated.
"""
import json
import math
from pathlib import Path

import pytest

pytestmark = pytest.mark.measurement
DATA = Path(__file__).resolve().parents[2] / "validation" / "rm2357"


def tables():
    return json.loads((DATA / "tests.json").read_text())["tables"]


def test_efficiency_consistent_with_kt_kq():
    """eta = J*KT/(2*pi*KQ). Two rows in the original's Table 1 disagree; everything else agrees."""
    known_bad = {(10.1, 0.302), (10.1, 0.400)}
    for tab in tables():
        for j, mt, kq, kt, eta in tab["rows"]:
            if eta is None:
                continue
            calc = j * kt / (2 * math.pi * kq)
            if (tab["theta07_deg"], j) in known_bad:
                assert abs(calc - eta) > 0.01
            else:
                assert calc == pytest.approx(eta, abs=6e-3), (tab["theta07_deg"], j)


def test_digitized_geometry_matches_report():
    geo = json.loads((DATA / "geometry.json").read_text())["stations"]
    r = [g["r"] for g in geo]
    cr = [g["cR"] for g in geo]
    area = sum(0.5 * (cr[i] + cr[i + 1]) * (r[i + 1] - r[i]) for i in range(len(r) - 1))
    assert area / math.pi == pytest.approx(0.0338, rel=0.02)          # solidity per blade
    tc07 = next(g["tc"] for g in geo if g["r"] == 0.7)
    assert tc07 == pytest.approx(0.093, abs=0.002)                    # 9.3% thick at 0.7R


@pytest.mark.skip(reason="needs pyxrotor OPER")
def test_pyxrotor_against_tunnel_data():
    """Target once OPER exists: KT within ~3% over the normal operating range (J >= 0.4,
    eta >= 0.5), matching what XROTOR 7.55 itself achieves with the same estimated sections."""
