"""XROTOR against NASA's SR-2 (straight) and SR-3 (45 deg swept) eight-blade propfan models in
the NASA Lewis 8x6 ft tunnel, Mach 0.6-0.85 (validation/nasa_propfan; NASA TM-79124,
CR-3505, TM X-73612). What these tests pin:
  * the traced blade geometry reproduces the published activity factors and integrated
    design lift coefficients;
  * SR-2 at the design point (J 3.06, CP 1.7): 7.55 with Mcrit = 0.85 - 1.2 t loses 0.3-0.75
    in efficiency where the test loses under 0.05; Mcrit + 0.20 (as for the thin 16-series
    propellers of RM L51L28) follows the measured loss within 0.025 of the low-speed bias,
    Mcrit + 0.15 (the P-47 correction) loses 0.07-0.09 too much at Mach 0.8-0.85;
  * SR-3: XROTOR, a straight lifting line, needs 5-8.5 deg less blade angle than the test for
    the design power (Hamilton Standard's lifting-line programs were 4-5 deg low);
  * SR-3 at fixed pitch: the test's CP falls 9% and 25% from Mach 0.70 to 0.80 and 0.85
    (CR-3505); XROTOR's rises 15-40% with every compressibility setting (its Prandtl-Glauert
    lift slope, with the local Mach number capped at 0.995).
"""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.measurement
HERE = Path(__file__).resolve().parents[2] / "validation" / "nasa_propfan"
_spec = importlib.util.spec_from_file_location("nasa_propfan_model", HERE / "model.py")
model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(model)


def _integrals(prop):
    g = model.GEO[prop]
    x = np.array(g["r_over_R"])
    xx = np.linspace(x[0], 1.0, 400)
    af = 100000.0 / 16.0 * np.trapezoid(np.interp(xx, x, g["b_over_D"]) * xx ** 3, xx)
    cli = 4.0 * np.trapezoid(np.interp(xx, x, g["design_cl"]) * xx ** 3, xx)
    return af, cli


@pytest.mark.parametrize("prop, af, cli", [("SR-2", 203, 0.08), ("SR-3", 235, 0.214)])
def test_traced_geometry_matches_published_integrals(prop, af, cli):
    """Activity factor (integrated from the spinner, 0.24R) within 1% and integrated design
    CL within 0.006 of the published values; the twist reads 0 at 0.75R."""
    a, c = _integrals(prop)
    assert abs(a / af - 1.0) < 0.01, a
    assert abs(c - cli) < 0.006, c
    g = model.GEO[prop]
    assert abs(np.interp(0.75, g["r_over_R"], g["twist_deg_from_0.75R"])) < 0.15


@pytest.fixture(scope="module")
def sr2_design():
    rows = [r for r in model.table("design_point_eta.csv")
            if r["prop"] == "SR-2" and r["source"] == "TM-79124 fig. 14"]
    out = {}
    for name, kw in {"7.55": {}, "shift15": {"mcrit_shift": 0.15},
                     "shift20": {"mcrit_shift": 0.20}}.items():
        res = {}
        for r in rows:
            o = model.run("SR-2", float(r["M"]), 3.06, cp=1.7, **kw)
            assert o["usable"], (name, r["M"])
            res[float(r["M"])] = (o["eta"], float(r["eta"]) / 100.0, o["beta75"])
        out[name] = res
    return out


def test_sr2_755_overpredicts_losses(sr2_design):
    for M, (eta, test, _) in sr2_design["7.55"].items():
        if M >= 0.7:
            assert eta - test < -0.25, M


def test_sr2_mcrit_plus_020_follows_the_loss(sr2_design):
    bias = sr2_design["shift20"][0.6][0] - sr2_design["shift20"][0.6][1]
    assert 0.0 < bias < 0.02, bias
    for M, (eta, test, _) in sr2_design["shift20"].items():
        assert abs(eta - test - bias) < 0.025, M
    for M in (0.8, 0.85):
        eta, test, _ = sr2_design["shift15"][M]
        assert eta - test - bias < -0.06, M


def test_sr3_blade_angle_too_low():
    """CR-3505: the design CP and J at Mach 0.80 need beta(0.75R) = 61.7 deg (Table III is
    one degree lower). XROTOR needs 5-8.5 deg less."""
    for kw in ({"mach_scale": 0.01}, {"mcrit_shift": 0.15}, {"sweep": "cosine", "mcrit_shift": 0.10}):
        o = model.run("SR-3", 0.80, 3.06, cp=1.695, **kw)
        assert o["usable"], kw
        assert 5.0 < 61.7 - o["beta75"] < 8.5, (kw, o["beta75"])


def test_sr3_fixed_pitch_power_rises_where_the_test_falls():
    """At fixed pitch (set to give the test's CP 1.64 at Mach 0.70, JREF 3.1) the model's CP
    rises with Mach number by 15% or more at 0.80 and 0.85; the test's falls to 0.909 and
    0.750 of its Mach 0.70 value."""
    for kw in ({"mcrit_shift": 0.15}, {"mcrit_shift": 0.20}, {"sweep": "cosine", "mcrit_shift": 0.10}):
        ref = model.run("SR-3", 0.70, 3.1, cp=1.64, **kw)
        assert ref["usable"]
        cp70 = model.run("SR-3", 0.70, 3.1, beta75=ref["beta75"], **kw)["CP"]
        for M in (0.80, 0.85):
            o = model.run("SR-3", M, 3.1, beta75=ref["beta75"], **kw)
            assert o["usable"], (kw, M)
            assert o["CP"] / cp70 > 1.14, (kw, M, o["CP"] / cp70)
    assert 1.49 / 1.64 < 0.92 and 1.23 / 1.64 < 0.76


SR3_TABLE_I = {0.6: 0.802, 0.7: 0.805, 0.8: 0.788, 0.85: 0.770}


def test_sr3_core_sweep_effective_sweep_is_the_published():
    """The calibrated quarter-chord line gives the core the published aerodynamic sweep at
    the design point, within 0.5 deg."""
    g = model.GEO["SR-3"]
    o = model.run("SR-3", 0.8, 3.06, cp=1.695, sweep="core", mach_scale=0.01)
    r = o["rotor"]
    x = np.array(r.s.xi[: r.s.ii])
    want = np.abs(np.interp(x, g["r_over_R"], g["sweep_deg"]))
    assert np.max(np.abs(r.effective_sweep_deg() - want)) < 0.5


def test_sr3_core_sweep_design_point():
    """With the core's sweep model (pyxrotor.sweep) and Mcrit + 0.10 every design-point solve
    converges, and the efficiency follows the test's within 0.02 of the Mach 0.6 bias to
    Mach 0.85; without a shift it loses 0.2 too much at 0.85."""
    e10 = {M: model.run("SR-3", M, 3.06, cp=1.695, sweep="core", mcrit_shift=0.10)
           for M in SR3_TABLE_I}
    assert all(o["usable"] for o in e10.values())
    bias = e10[0.6]["eta"] - SR3_TABLE_I[0.6]
    for M, t in SR3_TABLE_I.items():
        assert abs(e10[M]["eta"] - t - bias) < 0.02, M
    e0 = model.run("SR-3", 0.85, 3.06, cp=1.695, sweep="core")
    assert e0["eta"] - SR3_TABLE_I[0.85] - bias < -0.15


def test_sr3_core_sweep_moves_both_discrepancies_part_way():
    """Stages 1-2 of the sweep model (section-level sweep from the blade's geometry) move the
    blade angle and the fixed-pitch Mach trend the right way, not all the way: the blade angle
    at the test's CP is still 4-5 deg low (6 deg straight), and at fixed pitch CP still rises
    6-14% from Mach 0.70 to 0.85 (straight: 22%; test: falls 25%). What is left needs the
    induction of the swept lifting line and a transonic lift loss."""
    out = {}
    for name, kw in (("straight", {"mcrit_shift": 0.15}),
                     ("core", {"sweep": "core", "mcrit_shift": 0.10})):
        ref = model.run("SR-3", 0.70, 3.1, cp=1.64, **kw)
        b = ref["beta75"]
        c70 = model.run("SR-3", 0.70, 3.1, beta75=b, **kw)["CP"]
        c85 = model.run("SR-3", 0.85, 3.1, beta75=b, **kw)["CP"]
        out[name] = (b - 60.5, c85 / c70)
    assert -5.0 < out["core"][0] < -4.0 and out["straight"][0] < -5.5
    assert 1.05 < out["core"][1] < 1.15 and out["straight"][1] > 1.18


def test_sr3_curved_blade_induction_and_lift_loss():
    """Stage 3 of the sweep model on SR-3 (fixed pitch, J 3.1): the curved blade's own
    induction (sweep="core-ind") lowers the blade angle at the test's CP by 0.75 deg (the
    swept-back tip is loaded by upwash), widening the deficit to 5.3 deg. Hamilton Standard's
    sweep-aware lifting lines were 4-5 deg low too. The fixed-pitch power still rises with
    Mach number (+5-11% at 0.85). The transonic lift loss (lift_loss=0.10) barely acts on the
    swept blade, whose normal Mach numbers stay low. On the straight blade it drives the
    power down far past the test (CP at Mach 0.85 below 0.1 of Mach 0.70's, test 0.75),
    through the thick inboard sections. Neither model explains SR-3's trend; blade
    deflection under load (not modeled) remains."""
    def trend(**kw):
        b = model.run("SR-3", 0.70, 3.1, cp=1.64, **kw)["beta75"]
        cp = {M: model.run("SR-3", M, 3.1, beta75=b, **kw)["CP"] for M in (0.70, 0.85)}
        return b - 60.5, cp[0.85] / cp[0.70]
    db, ratio = trend(sweep="core-ind", mcrit_shift=0.10)
    assert -5.6 < db < -5.0 and 1.05 < ratio < 1.15
    db_ll, ratio_ll = trend(sweep="core-ind", mcrit_shift=0.10, lift_loss=0.10)
    assert abs(db_ll - db) < 0.2 and ratio_ll > 1.05
    _, straight_ll = trend(mcrit_shift=0.15, lift_loss=0.15)
    assert straight_ll < 0.1
