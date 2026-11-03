"""OPER case sweeps (SETCAS + CASPLT's solves), ADDC/CLRC/CASE, CPUT/CGET against the reference."""
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor.cases import (add_case, case_lines, cases_from_lines, clear_cases, rerun_case,
                            run_cases, setcas)
from pyxrotor.constants import ICASX, PI
from rotors import perturbed_states, solved_state
from statesync import assert_states_match, state_to_ref
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)
from conftest import EXPENSIVE  # noqa: E402

FIX = {"A": 0, "T": 1, "Q": 2, "P": 3}


@st.composite
def sweep_args(draw, max_cases=40):
    """(itype, first, last, step, options) roughly around the Rotol operating point."""
    itype = draw(st.sampled_from([1, 2, 3, 4]))
    lo, hi = {1: (0.05, 0.5), 2: (300.0, 2500.0), 3: (5.0, 40.0), 4: (10.0, 120.0)}[itype]
    first = draw(finite(lo, hi))
    last = draw(st.one_of(st.just(first), finite(lo, hi)))
    n = draw(st.integers(1, max_cases))
    step = draw(st.one_of(st.just(0.0), st.just((last - first) / n), finite(-hi, hi).map(lambda v: v / 7)))
    opts = dict(engine_line=draw(st.booleans()), fix=draw(st.sampled_from("ATQP")),
                constant_speed=draw(st.booleans()),
                rpm=draw(st.one_of(st.none(), finite(500.0, 2000.0))),
                power=draw(st.one_of(st.none(), finite(1e4, 1e6))))
    return itype, first, last, step, opts


def ref_setcas(ref, itype, first, last, step, o):
    ref.call("REFSETCAS", itype, first, last, step, o["engine_line"], FIX[o["fix"]],
             o["constant_speed"], o["rpm"] or 0.0, o["power"] is not None, o["power"] or 0.0)


@SETTINGS
@given(s=perturbed_states(), args=sweep_args(), ncase0=st.sampled_from([0, 0, 3, ICASX - 2, ICASX]),
       lpwrvar=st.booleans(), pspec=st.sampled_from([0.0, -1.0, 5e4]))
def test_setcas(ref_module, s, args, ncase0, lpwrvar, pspec):
    itype, first, last, step, o = args
    s.ncase, s.lpwrvar, s.pspec = ncase0, lpwrvar, pspec
    state_to_ref(s, ref_module)
    ref_setcas(ref_module, itype, first, last, step, o)
    setcas(s, itype, first, last, step, **o)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=perturbed_states(), first=finite(0.1, 0.4), step=st.sampled_from([1e-9, 1e-12, 1e-300, 5e-324]))
def test_setcas_huge_count(ref_module, s, first, step):
    """A tiny increment overflows NP = INT(...) + 1 (INTEGER*4): no cases are added."""
    s.ncase = 0
    state_to_ref(s, ref_module)
    ref_setcas(ref_module, 1, first, first + 0.5, step, dict(engine_line=False, fix="A",
               constant_speed=False, rpm=None, power=None))
    setcas(s, 1, first, first + 0.5, step)
    assert_states_match(s, ref_module)


def rotol_with_engine_line(s, rng):
    rpm = s.vel / (s.rad * s.adv * PI / 30.0)
    dim = s.rho * s.vel**3 * s.rad**2
    n = 4
    s.lpwrvar, s.npwrvar = True, n
    s.rpmvar[:n] = np.sort(rpm * rng.uniform(0.5, 1.5, n)) + np.arange(n)
    s.pwrvar[:n] = abs(s.ptot) * dim * rng.uniform(0.7, 1.3, n)
    s.xpwrvar[:n] = rng.normal(0.0, abs(s.ptot) * dim / max(rpm, 1.0), n)


@EXPENSIVE
@given(sweeps=st.lists(sweep_args(max_cases=4), min_size=1, max_size=2), addc=st.booleans(),
       engine=st.booleans(), seed=st.integers(0, 2**32 - 1), theta=st.sampled_from([15.0, 22.0, 27.0]),
       free=st.booleans())
def test_run_cases(ref_module, sweeps, addc, engine, seed, theta, free):
    """Sweeps as OPER runs them: SETCAS then CASPLT's solves (all stored cases, appended),
    with ADDC cases (no type code) mixed in."""
    from rotors import _solved
    s = _solved(theta, 0.55, 0.41).copy()
    s.free = free
    if engine:
        rotol_with_engine_line(s, np.random.default_rng(seed))
    s.ncase = 0
    s.caspar[0, :] = 0.0
    state_to_ref(s, ref_module)
    for itype, first, last, step, o in sweeps:
        if addc:
            ref_module.call("REFADDC")
            add_case(s)
        ref_setcas(ref_module, itype, first, last, step, o)
        setcas(s, itype, first, last, step, **o)
        assert_states_match(s, ref_module)
        if s.ncase > 6:                        # keep the solves affordable
            continue
        ref_module.call("REFCASRUN")
        run_cases(s)
        assert_states_match(s, ref_module)


@EXPENSIVE
@given(args=sweep_args(max_cases=3), icase=st.integers(0, 4), clear=st.booleans())
def test_case_commands(ref_module, args, icase, clear):
    s = solved_state()
    s.ncase = 0
    state_to_ref(s, ref_module)
    itype, first, last, step, o = args
    ref_setcas(ref_module, itype, first, last, step, o)
    setcas(s, itype, first, last, step, **o)
    ref_module.call("REFCASE", icase)
    rerun_case(s, icase)
    assert_states_match(s, ref_module)
    if clear:
        clear_cases(s)
        ref_module["NCASE"], ref_module["KCASE"] = 0, 0      # CLRC is two assignments
        assert_states_match(s, ref_module)


@EXPENSIVE
@given(args=sweep_args(max_cases=5), name=st.sampled_from(["", "Rotol", "x" * 40]))
def test_cput_cget(ref_module, tmp_path_factory, args, name):
    """CPUT writes exactly SHOCAS's text; CGET reads it back (losing results, as 7.55 does)."""
    d = tmp_path_factory.mktemp("cases")
    s = solved_state()
    s.ncase, s.name = 0, name
    state_to_ref(s, ref_module)
    itype, first, last, step, o = args
    ref_setcas(ref_module, itype, first, last, step, o)
    setcas(s, itype, first, last, step, **o)
    if 0 < s.ncase <= 5:
        ref_module.call("REFCASRUN")
        run_cases(s)
    fref = d / "ref.cas"
    ref_module["FNAME"] = str(fref)
    ref_module.call("REFCPUT")
    lines = case_lines(s.caspar, s.ncase, s.rad, s.name)
    assert fref.read_text() == "\n".join(lines) + "\n"

    ref_module["NCASE"] = 0
    ref_module.call("REFCGET")
    cases_from_lines(s, lines)
    assert_states_match(s, ref_module)


def test_cget_full_lines_fix(tmp_path):
    """7.55's CGET drops power/thrust/torque/efficiency (128-character records);
    full_lines=True keeps them, to the 5 digits the file holds."""
    s = solved_state()
    s.ncase = 0
    setcas(s, 1, 0.3, 0.5, 0.1)
    run_cases(s)
    lines = case_lines(s.caspar, s.ncase, s.rad, "t")
    solved = s.caspar[:, : s.ncase].copy()
    assert (solved[8:12] != 999.0).all()
    cases_from_lines(s, lines)
    assert (s.caspar[8:12, : s.ncase] == 999.0).all()
    cases_from_lines(s, lines, full_lines=True)
    np.testing.assert_allclose(s.caspar[8:12, : s.ncase], solved[8:12], rtol=1e-4)
    np.testing.assert_allclose(s.caspar[1:8, : s.ncase], solved[1:8], rtol=1e-4)


def test_rotor_sweeps_and_files(tmp_path):
    """The Rotor wrappers: a session of sweeps against the reference, then CPUT/SAVE/LOAD."""
    from test_rotor import rotol_both
    ref, rotor = rotol_both(22.0)
    ref.call("REFADVA", 0.3)
    rotor.solve_advance_ratio(0.3)
    for itype, a, b, c, opts, method in [
            (1, 0.2, 0.4, 0.1, {}, lambda: rotor.sweep_advance_ratio(0.2, 0.4, 0.1)),
            (2, 900.0, 1100.0, 100.0, {"fix": "T"}, lambda: rotor.sweep_rpm(900, 1100, 100, hold="thrust")),
            (3, 20.0, 24.0, 2.0, {}, lambda: rotor.sweep_blade_angle(20, 24, 2)),
            (4, 40.0, 60.0, 10.0, {"constant_speed": True}, lambda: rotor.sweep_velocity(40, 60, 10, constant_speed=True))]:
        o = dict(engine_line=False, fix="A", constant_speed=False, rpm=None, power=None) | opts
        ref_setcas(ref, itype, a, b, c, o)
        ref.call("REFCASRUN")
        new = method()
        assert len(new) == 3
        assert_states_match(rotor.s, ref)
        if itype != 4:
            rotor.clear_cases()
            ref["NCASE"], ref["KCASE"] = 0, 0
    assert all(c["efficiency"] is not None for c in new)
    rotor.write_cases(tmp_path / "c.cas")
    assert rotor.read_cases(tmp_path / "c.cas") == 3
    rotor.save(tmp_path / "r.xr")
    r2 = __import__("pyxrotor").Rotor()
    assert r2.load(tmp_path / "r.xr")
    assert abs(r2.performance()["efficiency"] - rotor.performance()["efficiency"]) < 1e-3
