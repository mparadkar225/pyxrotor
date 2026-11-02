"""ARC R&M 2357: the Rotol wooden Spitfire IX propeller in the RAE 24-ft tunnel.

Data-integrity checks, and pyxrotor against the tunnel data at engineering tolerances (the
section polars are estimated, since the report gives only the section shape).
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


def test_pyxrotor_against_tunnel_data():
    """pyxrotor (= XROTOR 7.55, bit for bit) vs the tunnel, at each tabulated J and tip Mach.

    Normal operating range: eta >= 0.5 and KT >= 0.02 (relative error is meaningless near
    zero thrust). With the estimated Clark Y sections: thrust agrees to ~3% on average;
    torque is underpredicted (the section model lacks the blade's real profile and
    high-lift drag), so efficiency is high. These bounds pin that known behavior: a change
    that moves the port away from it should be noticed, whichever direction it goes.
    """
    import numpy as np

    from pyxrotor.aero import AeroSection
    from pyxrotor.constants import PI
    from pyxrotor.rotor import Rotor

    sections = [  # r/R, a0 (deg), dCL/da, CLmax, CLmin, CDmin, CL@CDmin, dCD/dCL^2, Mcrit
        (0.0, -5.0, 5.5, 1.10, -0.5, 0.020, 0.30, 0.010, 0.55),
        (0.4, -4.5, 5.9, 1.30, -0.6, 0.011, 0.35, 0.008, 0.66),
        (0.7, -3.5, 6.1, 1.35, -0.6, 0.0085, 0.40, 0.007, 0.74),
        (1.0, -2.5, 6.2, 1.30, -0.5, 0.0075, 0.35, 0.007, 0.79),
    ]
    geo = json.loads((DATA / "geometry.json").read_text())["stations"]
    radius = 10.75 * 0.3048 / 2
    e_kt, e_kq, e_eta = [], [], []
    for tab in tables():
        for j, mt, kq, kt, eta in tab["rows"]:
            if eta is None or eta < 0.5 or kt < 0.02:
                continue
            r = Rotor()
            r.atmosphere(0.0)
            r.set_sections([AeroSection(xisect=x, a0=a0 * PI / 180, dclda=d, clmax=cmx, clmin=cmn,
                                        dclda_stall=0.1, dcl_stall=0.1, cdmin=cd, cldmin=cl0,
                                        dcdcl2=k, cmcon=-0.1, mcrit=mc, reref=2.0e6, rexp=-0.15)
                            for x, a0, d, cmx, cmn, cd, cl0, k, mc in sections])
            r.arbi(2, 50.0, radius, 0.2 * radius, [g["r"] for g in geo], [g["cR"] for g in geo],
                   [tab["theta07_deg"] + g["dth"] for g in geo])
            n = mt * r.s.vso / (2 * math.pi * radius)
            r.velocity(j * n * 2 * radius)
            assert r.solve_rpm(n * 60.0)
            p = r.performance()
            e_kt.append(p["CT"] / kt - 1)
            e_kq.append(p["CP"] / (2 * math.pi) / kq - 1)
            e_eta.append(p["efficiency"] - eta)
    e_kt, e_kq, e_eta = map(np.array, (e_kt, e_kq, e_eta))
    assert len(e_kt) == 24
    assert np.mean(np.abs(e_kt)) < 0.04 and np.max(np.abs(e_kt)) < 0.15
    assert np.all((e_kq < 0.0) & (e_kq > -0.25))          # torque low, by 3-19%
    assert np.all((e_eta > 0.0) & (e_eta < 0.15))         # efficiency high, by 0.03-0.13
