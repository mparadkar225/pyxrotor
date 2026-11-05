"""The transonic solves of pyxrotor.robust. Checks against independent results:
  * where XROTOR's POWE converges, solve_power is POWE itself, bit for bit;
  * the fixed-pitch scan's interpolated blade angle and efficiency match POWE's converged
    solution where both exist;
  * the Mach continuation ends at a converged XROTOR solution of the requested power, with
    the speed of sound restored;
  * a stalled Newton iteration (no spread, but no solution: RLX 0) is not accepted as a
    limit cycle, while the small limit cycles of the SR-3 propfan at Mach 0.85 are;
  * a request nothing can solve leaves the rotor as it was.
Cases come from the measurement models in validation/ (transonic propellers)."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from pyxrotor import robust
from pyxrotor.constants import PI

pytestmark = pytest.mark.model
VAL = Path(__file__).resolve().parents[2] / "validation"


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, VAL / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


yp47m = _load("robust_yp47m_model", "naca_yp47m/model.py")
l51 = _load("robust_l51l28_model", "naca_l51l28/model.py")
propfan = _load("robust_propfan_model", "nasa_propfan/model.py")


def yp47m_case(prop, M0, J, cp, shift):
    guess = np.degrees(np.arctan(J / (0.75 * PI))) + 5.0
    r, n, D = yp47m._rotor(prop, M0, J, float(guess), shift, 1.0)
    return r, cp * yp47m.RHO_ALT * n ** 3 * D ** 5, 60.0 * n


def same_state(a, b):
    for k, v in vars(a).items():
        w = getattr(b, k)
        if isinstance(v, np.ndarray):
            if not np.array_equal(v, w, equal_nan=True):
                return False
        elif k != "fixes" and v != w and not (v != v and w != w):
            return False
    return True


CASE_POWE = ("curtiss_836-14C2-18R1", 0.4, 1.7, 0.492, 0.15)


def test_powe_when_it_converges():
    r, power, rpm = yp47m_case(*CASE_POWE)
    ref = r.s.copy()
    sol = robust.solve_power(r, power, rpm)
    r2, _, _ = yp47m_case(*CASE_POWE)
    assert r2.solve_power(power, fix="rpm", rpm=rpm)
    assert sol.method == "powe" and sol.exact and sol.usable and sol.solves == 1
    assert same_state(r.s, r2.s)
    assert not same_state(r.s, ref)
    assert sol.efficiency == r2.performance()["efficiency"]


@pytest.mark.parametrize("case", [CASE_POWE, ("curtiss_732-1C2-0", 0.5, 2.2, 0.616, 0.0),
                                  ("aeroproducts_H20C-162", 0.4, 1.75, 0.5, 0.15)])
def test_scan_interpolation_matches_powe(case):
    """0.5 deg grid: blade angle within 0.03 deg, efficiency within 0.002 of POWE."""
    r, power, rpm = yp47m_case(*case)
    exact = robust.solve_power(r, power, rpm)
    assert exact.method == "powe"
    r, power, rpm = yp47m_case(*case)
    sol = robust.scan_power(r, power, rpm, exact.beta75_deg + 0.2)
    assert sol.method == "scan" and sol.usable and not sol.exact
    assert sol.power == power and sol.bracket_deg == pytest.approx(0.5)
    assert abs(sol.beta75_deg - exact.beta75_deg) < 0.03
    assert abs(sol.efficiency - exact.efficiency) < 0.002
    assert abs(sol.thrust / exact.thrust - 1.0) < 0.003
    # the rotor is left at the lower bracketing solve, a converged fixed-pitch solution
    assert r.s.conv and robust.beta75_deg(r) == pytest.approx(sol.beta75_deg - 0.5, abs=0.5)


def test_continuation_reaches_a_converged_solution():
    """Curtiss 836, Mach 0.5, J 1.75 (tip Mach 1.03), Mcrit + 0.20: POWE from the starting
    blade angle fails, POWE continued in Mach number converges."""
    r, power, rpm = yp47m_case("curtiss_836-14C2-18R1", 0.5, 1.75, 0.4, 0.20)
    vso = r.s.vso
    probe = r.s.copy()
    assert not yp47m_case("curtiss_836-14C2-18R1", 0.5, 1.75, 0.4, 0.20)[0].solve_power(
        power, fix="rpm", rpm=rpm)
    sol = robust.solve_power(r, power, rpm)
    assert sol.method == "continuation" and sol.exact
    assert r.s.conv and r.s.vso == vso == probe.vso
    p = r.performance()
    assert abs(p["power"] / power - 1.0) < 1e-6 and abs(p["rpm"] / rpm - 1.0) < 1e-12
    assert sol.efficiency == p["efficiency"]


def test_stalled_iteration_is_not_a_limit_cycle():
    """RM L51L28 symmetric blade at -2 deg, static, tip Mach 1.0: XROTOR's iteration stops
    moving (RLX 0) far from a solution. Its iterates don't change, but it isn't accepted."""
    o = l51.static(0.0, -2.0, 1.0)
    assert not o["converged"]
    r = o["rotor"]
    sol = robust.solve_fixed_pitch(r, r.performance()["rpm"])
    assert not sol.usable and sol.method is None
    assert sol.limit_cycle["dpower_rel"] < 1e-3 and sol.limit_cycle["min_rlx"] < 1e-6
    assert sol.limit_cycle["max_step_rms"] > 0.1


def test_propfan_limit_cycles_are_accepted():
    """SR-3 at Mach 0.85, J 3.06, Mcrit + 0.15: RPM doesn't converge near its design blade
    angle, but settles into a cycle with efficiency within 0.005 and power within 2%."""
    ref = propfan._solve("SR-3", 0.85, 3.06, 1.695, 58.0, 0.15, None, True, 0.01)
    o = propfan._solve("SR-3", 0.85, 3.06, None, ref["beta75"], 0.15, None, True, 1.0)
    assert not o["converged"]
    r = o["rotor"]
    sol = robust.solve_fixed_pitch(r, r.performance()["rpm"])
    assert sol.usable and sol.method == "cycle" and not sol.exact
    c = sol.limit_cycle
    assert c["deta"] < 0.005 and c["dpower_rel"] < 0.02 and c["min_rlx"] >= 0.05


def test_failure_restores_the_rotor():
    """A power no blade angle in the scan can absorb: nothing usable, rotor unchanged."""
    r, power, rpm = yp47m_case(*CASE_POWE)
    start = r.s.copy()
    sol = robust.solve_power(r, 1e3 * power, rpm, continuation=None, scan=(-1.0, 1.0, 0.5))
    assert not sol.usable and sol.method is None and sol.solves > 1
    assert same_state(r.s, start)
