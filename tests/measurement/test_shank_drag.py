"""Shank drag (pyxrotor.shank, TR-852's cylinder) against three tests:
  * NACA TR-852's own shank-fairing sections (NACA 16-025 and 16-040 in flight, Mach
    0.4-0.82): the compressible drag rise of the 25% section follows the calibrated Mcrit
    within 0.015. For the 40% section it follows 0.85 - 1.2 t + 0.15 without the 0.5 floor
    within 0.03 above Mach 0.63. The floor that pyxrotor.sections keeps starts it too late,
    by 0.1 or more in CD at Mach 0.74-0.79;
  * the YP-39 flight test (Hamilton Standard 3155-6, round shanks to 0.25R, faired by 0.30R):
    the entire propeller's efficiency within 0.035 in every tip-Mach band, against +0.11
    without the shank model. The modeled shank loss grows from 0.02 to 0.20 with speed, as
    the test's does (0.02 to 0.24);
  * RM L6L27a, two full-scale propellers differing only in shank, at constant power: the
    efficiency difference within 0.7 points of the test's up to 375 mph (tip Mach 0.72).
    Without the shank model it is 3.7-6.5 points too small from 325 mph. Above 400 mph the oval
    shank drawn as a cylinder of its thickness gives too much (10.6-19 against 8.4-11.7).
"""
import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from pyxrotor import sections
from pyxrotor.aero import AeroSection, clcdcm

pytestmark = pytest.mark.measurement
VAL = Path(__file__).resolve().parents[2] / "validation"


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, VAL / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _cd(sec, mach):
    return clcdcm(0.0, mach / 0.5, -1.0, 0.5, 1.0, sec).cdrag


def test_tr852_shank_sections():
    rows = list(csv.DictReader(open(VAL / "naca_tr852" / "shank_sections_fig16.csv")))
    by = {}
    for r in rows:
        by.setdefault(r["section"], []).append((float(r["M"]), float(r["cd_on_thickness"]) * float(r["t"])))
    # 16-025: rise above each curve's low-speed level
    cal = sections.naca16(0.25, 0.3, 0.0)
    lo_t, lo_m = by["NACA 16-025"][0][1], _cd(cal, 0.5)
    for m, cd in by["NACA 16-025"]:
        assert abs((_cd(cal, m) - lo_m) - (cd - lo_t)) < 0.015, m
    # 16-040: floorless rule fits, the floored one is late
    floored = sections.naca16(0.40, 0.3, 0.0)
    floorless = AeroSection(**{**floored.__dict__, "mcrit": 0.85 - 1.2 * 0.40 + 0.15})
    lo_t = by["NACA 16-040"][0][1]
    for m, cd in by["NACA 16-040"]:
        rise = cd - lo_t
        if m > 0.63:
            assert abs((_cd(floorless, m) - _cd(floorless, 0.4)) - rise) < 0.03, m
        if m > 0.74:
            assert rise - (_cd(floored, m) - _cd(floored, 0.4)) > 0.1, m


def test_yp39_entire_propeller():
    yp39 = _load("shank_yp39_model", "naca_yp39/model.py")
    test_loss, model_loss, err = [], [], {}
    for row in yp39.runs():
        o = yp39.run(row, mcrit_shift=0.15, shank=(0.25, 0.30))
        assert o["converged"]
        err.setdefault(min((b for b in (0.99, 1.035, 1.07) if o["Mt"] < b)), []).append(
            o["eta"] - float(row["eta_entire"]))
        test_loss.append(float(row["eta_entire"]) - float(row["eta_working"]))
        model_loss.append(o["eta"] - o["eta_blades"])
    assert all(abs(np.mean(v)) < 0.035 for v in err.values()), err
    assert model_loss[0] > -0.03 and model_loss[-1] < -0.18
    assert np.max(np.abs(np.array(model_loss) - np.array(test_loss))) < 0.045


def test_l6l27a_shank_difference():
    cmp = _load("shank_l6l27a_compare", "naca_l6l27a/compare.py")
    cmp.SETTINGS = {"calibrated": {}, "no shank": {"shank": None}}
    for r in cmp.fig28():
        d, d0, t = r["calibrated"]["percent"], r["no shank"]["percent"], r["test_percent"]
        if r["airspeed_mph"] <= 375:
            assert abs(d - t) < 0.7, r["airspeed_mph"]
        else:
            assert d > t + 2.0, r["airspeed_mph"]
        if r["airspeed_mph"] >= 325:
            assert t - d0 > 2.0, r["airspeed_mph"]
