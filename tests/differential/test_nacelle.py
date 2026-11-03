"""NACE (NACELL, SRCVEL) and OPER ADW (ADWFCTR) against the reference."""
import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from pyxrotor.fortio import fmt_g
from pyxrotor.nacelle import nacelle, srcvel
from rotors import RADIUS, perturbed_states, solved_state
from statesync import assert_states_match, state_to_ref
from test_rotor import rotol_both
from strategies import finite
from test_restart import fresh, short_dir  # noqa: F401  (fixture)

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
from conftest import EXPENSIVE  # noqa: E402


@SETTINGS
@given(zb=st.lists(finite(-3.0, 3.0), min_size=0, max_size=20), z=finite(-2.0, 2.0),
       r=st.one_of(st.just(0.0), finite(0.0, 2.0)), data=st.data())
def test_srcvel(ref_module, zb, z, r, data):
    zb = np.sort(np.array(zb, dtype=float))
    n = len(zb)
    src = np.array(data.draw(st.lists(finite(-0.5, 0.5), min_size=n, max_size=n)), dtype=float)
    out = ref_module.call("SRCVEL", n, zb if n else np.zeros(1), src if n else np.zeros(1),
                          z, r, 0.0, 0.0)
    u, v = srcvel(zb, src, z, r)
    np.testing.assert_array_equal([u, v], [out[5], out[6]])


@st.composite
def body_files(draw):
    """A spinner/cowling: rotor near the nose, area growing then constant (or closing)."""
    n = draw(st.integers(0, 40))
    z = np.sort(np.array(draw(st.lists(finite(-0.5, 3.0), min_size=n, max_size=n, unique=True))))
    rmax = RADIUS * draw(finite(0.1, 0.4))
    shape = np.clip(z / max(draw(finite(0.2, 1.0)), 1e-3), 0.0, 1.0)
    radius = rmax * np.sqrt(shape) * draw(finite(0.8, 1.2))
    as_radius = draw(st.booleans())
    area = -radius if as_radius else np.pi * radius * radius
    zprop = draw(finite(-0.2, 0.5))
    awake = -(rmax * draw(finite(0.0, 1.0))) if draw(st.booleans()) else np.pi * (rmax * 0.5) ** 2
    style = draw(st.sampled_from(["repr", "g12", "split", "extra", "blanks", "commas"]))
    return zprop, awake, z, area, style


def write_body(path, zprop, awake, z, area, style):
    num = (lambda v: fmt_g(v, 12, 5)) if style == "g12" else repr
    lines = []
    head = [num(float(zprop)), num(float(awake))]
    lines += [head[0], head[1]] if style == "split" else [" ".join(head)]
    for zi, ai in zip(z, area):
        a, b = num(float(zi)), num(float(ai))
        if style == "extra":
            lines.append(f"{a} {b}  0.0 comment")
        elif style == "commas":
            lines.append(f"{a},{b}")
        elif style == "split":
            lines += [a, b]
        else:
            lines.append(f"  {a}   {b}")
        if style == "blanks":
            lines.append("")
    path.write_text("\n".join(lines) + "\n")


@SETTINGS
@given(s=perturbed_states(), body=body_files(), fix=st.booleans())
def test_nacelle_file(ref_module, short_dir, s, body, fix):  # noqa: F811
    f = fresh(short_dir, ".nac")
    write_body(f, *body)
    state_to_ref(s, ref_module)
    ref_module["FNAME"] = str(f)
    ref_module.call("REFNACELL", 0.0, fix)
    assert nacelle(s, f, fix_source_positions=fix)
    assert_states_match(s, ref_module)


@SETTINGS
@given(s=perturbed_states(), name=st.sampled_from(["CONST", "const", "Constant ubody"]),
       u=finite(-0.2, 0.2))
def test_nacelle_constant(ref_module, s, name, u):
    state_to_ref(s, ref_module)
    ref_module["FNAME"] = name
    ref_module.call("REFNACELL", u, False)
    assert nacelle(s, name, ubody=u)
    assert_states_match(s, ref_module)


BAD = {
    "missing": None,
    "bad_header": "0.1 abc\n0.0 0.0\n",
    "bad_station": "0.1 0.01\n0.0 0.0\n0.2 0.005\n0.3 x\n",
    "header_only": "0.1 0.01\n",
    "odd_count": "0.1 0.01\n0.0 0.0\n0.2\n",
}


@pytest.mark.parametrize("case", sorted(BAD))
def test_nacelle_bad_files(short_dir, case):  # noqa: F811
    from reference.ref755 import Ref755
    f = fresh(short_dir, ".nac")
    if BAD[case] is not None:
        f.write_text(BAD[case])
    s = solved_state()
    s.ubody[: s.ii] = 0.05
    ref = Ref755()
    state_to_ref(s, ref)
    ref["FNAME"] = str(f)
    ref.call("REFNACELL", 0.0, False)
    ok = nacelle(s, f)
    assert ok == (case in ("header_only", "odd_count"))
    assert_states_match(s, ref)


def test_nacelle_empty_file_stops(short_dir):  # noqa: F811
    """READ(LU,*,ERR=90) ZPROP, AWAKE has no END=: an empty file stops XROTOR."""
    import subprocess
    import sys
    from pathlib import Path
    f = fresh(short_dir, ".nac")
    f.write_text("")
    with pytest.raises(EOFError):
        nacelle(solved_state(), f)
    code = ("import sys; sys.path[:0] = ['tests/differential', 'tests', 'src', '.']\n"
            "from rotors import solved_state\nfrom statesync import state_to_ref\n"
            "from reference.ref755 import Ref755\n"
            f"r = Ref755(); state_to_ref(solved_state(), r); r['FNAME'] = {str(f)!r}\n"
            "r.call('REFNACELL', 0.0, False)\nprint('survived')")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, cwd=Path(__file__).resolve().parents[2])
    assert p.returncode != 0 and "End of file" in p.stderr and "survived" not in p.stdout


def spinner(fix=False):
    """A P-51-like spinner and cowling: 0.3 m spinner radius growing to a 0.45 m cowl."""
    z = np.linspace(0.0, 1.5, 16)
    r = 0.45 * np.sqrt(np.clip(z / 0.9, 0.0, 1.0))
    return dict(rotor_z=0.15, wake_area=-0.30, z=z, area=-r, fix_source_positions=fix)


@pytest.mark.parametrize("fix", [False, True])
def test_nacelle_then_solve(short_dir, fix):  # noqa: F811
    """A spinner body file, then an OPER rpm solve, matches; and the arrays API equals the file."""
    ref, rotor = rotol_both(22.0)
    b = spinner(fix)
    f = fresh(short_dir, ".nac")
    write_body(f, b["rotor_z"], b["wake_area"], b["z"], b["area"], "repr")
    ref["FNAME"] = str(f)
    ref.call("REFNACELL", 0.0, fix)
    rotor.nacelle(**b)
    assert_states_match(rotor.s, ref)
    assert (rotor.s.ubody[: rotor.s.ii] < 0.0).any()       # the body slows the root flow
    ref.call("REFRPM", 1000.0)
    rotor.solve_rpm(1000.0)
    assert_states_match(rotor.s, ref)


@EXPENSIVE
@given(factor=st.one_of(st.just(1.0), finite(0.8, 1.25)), free=st.booleans(),
       adv=finite(0.15, 0.35), theta=st.sampled_from([15.0, 22.0, 27.0]))
def test_adw_factor(factor, free, adv, theta):
    """OPER ADW sets ADWFCTR; an ADVA solve then matches (free wake uses it, rigid ignores it)."""
    ref, rotor = rotol_both(theta)
    ref["FREE"] = rotor.s.free = free
    ref["ADWFCTR"] = factor
    rotor.set_wake_advance_factor(factor)
    ref.call("REFADVA", adv)
    rotor.solve_advance_ratio(adv)
    assert_states_match(rotor.s, ref)
