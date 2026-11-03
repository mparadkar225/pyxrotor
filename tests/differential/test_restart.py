"""SAVE and LOAD restart files against the reference: byte-identical files, identical State."""
import contextlib
import os
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from pyxrotor.restart import LoadError, load, save, save_lines
from pyxrotor.rotor import Rotor
from reference.ref755 import Ref755
from rotors import perturbed_states, solved_state
from statesync import assert_states_match, state_to_ref
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
from conftest import EXPENSIVE  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "data"
_COUNT = iter(range(10**9))


@pytest.fixture(scope="module")
def short_dir():
    """XROTOR's FNAME is CHARACTER*80: keep paths short, and never reuse a name (SAVE
    prompts on stdin before overwriting)."""
    import shutil
    import tempfile
    d = tempfile.mkdtemp(prefix="xr", dir="/tmp")
    yield Path(d)
    shutil.rmtree(d, ignore_errors=True)


def fresh(d: Path, suffix=".xr") -> Path:
    p = d / f"{next(_COUNT)}{suffix}"
    assert len(str(p)) <= 80
    return p


@contextlib.contextmanager
def fortran_stdin(text: str):
    """Answer a Fortran READ(*) prompt: a pipe holding `text` on file descriptor 0."""
    r, w = os.pipe()
    os.write(w, text.encode())
    saved = os.dup(0)
    os.dup2(r, 0)
    try:
        yield
    finally:
        os.dup2(saved, 0)
        for fd in (saved, r, w):
            os.close(fd)


@st.composite
def saveable_states(draw):
    s = draw(perturbed_states())
    s.name = draw(st.sampled_from(["", "Rotol", "Designed blade", "a" * 40]))
    s.wind = draw(st.booleans())
    if draw(st.booleans()):                             # LVDUCT false: ADW = ADV*URDUCT
        s.adw = s.adv * s.urduct
    n = draw(st.sampled_from([0, 1, 2, 7]))
    s.nadd = n
    if n:
        s.radd[:n] = np.linspace(0.0, s.rad, n)
        s.uadd[:n] = draw(finite(-3.0, 3.0))
        s.vadd[:n] = draw(finite(-3.0, 3.0))
    return s


@SETTINGS
@given(s=saveable_states())
def test_save(ref_module, short_dir, s):
    """SAVE: the same bytes, and the same State afterwards (NAME default, GREEK)."""
    fr, fp = fresh(short_dir), fresh(short_dir)
    state_to_ref(s, ref_module)
    ref_module["FNAME"] = str(fr)
    ref_module.call("REFSAVE")
    save(s, fp)
    assert fp.read_bytes() == fr.read_bytes()
    assert_states_match(s, ref_module)
    assert ref_module.get("NAME") == s.name[:32].rstrip()      # NAME is CHARACTER*32


@EXPENSIVE
@given(src=saveable_states(), ii=st.sampled_from([None, None, 20, 45]),
       answer=st.sampled_from([None, 25, 30, 61]), free=st.booleans())
def test_load(short_dir, src, ii, answer, free):
    """LOAD of a file SAVE wrote, into a fresh session (start-up state, optionally another II).
    When the station counts differ, 7.55 prompts; `answer` (None: <Return>) is fed to it."""
    f = fresh(short_dir)
    save(src, f)
    ref = Ref755()
    rotor = Rotor()
    s = rotor.s
    if ii is not None:
        ref["II"], ref["IINF"] = ii, ii + ii // 2
        s.ii, s.iinf = ii, ii + ii // 2
    ref["FREE"] = s.free = free
    ref["FNAME"] = str(f)
    prompt = s.ii != src.ii
    with fortran_stdin(("" if answer is None else str(answer)) + "\n") if prompt else contextlib.nullcontext():
        ref.call("REFLOAD")
    load(s, f, stations=answer)
    assert_states_match(s, ref)


def test_load_xrotor_file(short_dir):
    """A file written by interactive XROTOR 7.55 (the P-51D design), loaded both ways."""
    import shutil
    f = fresh(short_dir)
    shutil.copy(DATA / "p51d_prop.xr", f)
    ref, rotor = Ref755(), Rotor()
    ref["FNAME"] = str(f)
    ref.call("REFLOAD")
    assert load(rotor.s, f)
    assert_states_match(rotor.s, ref)
    p = rotor.performance()
    assert round(p["efficiency"], 4) == 0.7978        # design.out, after LOAD
    assert rotor.s.free                               # from the file's LVDUCT flag


EDGE_FILES = {
    "comments_and_blanks": lambda t: t.replace("\n!", "\n\n# extra comment\n!", 3),
    "no_duct_velocity": lambda t: t[: t.index("!      URDuct")],
    "slipstream_short": lambda t: t + "!Nadd\n 5\n 0.1 1.0 0.0\n 0.5 1.1 0.1\n 0.9 1.2 0.2\n",
    "slipstream_one_bad": lambda t: t + "!Nadd\n 4\n 0.1 1.0 0.0\n x y z\n",
    "slipstream_long": lambda t: t + "!Nadd\n 3\n" + "".join(f" {0.1*k} 1.0 0.0\n" for k in range(1, 4)),
    "end_after_blade": lambda t: t + "END\n 3\n",
    "commas": lambda t: t.replace("   1.7020       112.00", "   1.7020  ,     112.00"),
}


@pytest.mark.parametrize("edit", sorted(EDGE_FILES))
def test_load_edge_files(short_dir, edit):
    text = EDGE_FILES[edit]((DATA / "p51d_prop.xr").read_text())
    f = fresh(short_dir)
    f.write_text(text)
    ref, rotor = Ref755(), Rotor()
    ref["FNAME"] = str(f)
    ref.call("REFLOAD")
    load(rotor.s, f)
    assert_states_match(rotor.s, ref)


def test_load_end_marker_aborts(short_dir):
    """A line reading END where the duct velocity should be: LOAD skips URDUCT, then reads the
    URDUCT value 1.0000 as the integer NADD and stops with a Fortran runtime error (that
    READ has END= but no ERR=). The port raises LoadError."""
    import subprocess
    import sys
    text = (DATA / "p51d_prop.xr").read_text().replace("!      URDuct\n", "END\n")
    f = fresh(short_dir)
    f.write_text(text)
    with pytest.raises(LoadError):
        load(Rotor().s, f)
    code = ("import sys; sys.path[:0] = ['.']\n"
            "from reference.ref755 import Ref755\n"
            f"r = Ref755(); r['FNAME'] = {str(f)!r}; r.call('REFLOAD')\n"
            "print('survived')")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, cwd=Path(__file__).resolve().parents[2])
    assert p.returncode != 0 and "Bad integer" in p.stderr and "survived" not in p.stdout


def test_save_load_round_trip_free_flag(tmp_path):
    """7.55: LOAD reads SAVE's LVDUCT as FREE; keep_free_wake=True keeps the session's flag."""
    s = solved_state()
    s.adw = s.adv * s.urduct                        # LVDUCT = .FALSE.
    save(s, tmp_path / "r.xr")
    r = Rotor()
    r.s.free = True
    load(r.s, tmp_path / "r.xr")
    assert not r.s.free
    r2 = Rotor()
    r2.s.free = True
    load(r2.s, tmp_path / "r.xr", keep_free_wake=True)
    assert r2.s.free


def test_load_errors(tmp_path):
    f = tmp_path / "bad.xr"
    f.write_text("XROTOR VERSION:  7.55\nname\n 1.0 2.0\n")
    with pytest.raises(LoadError):
        load(Rotor().s, f)
    f.write_text("old format first line\n")
    with pytest.raises(NotImplementedError):
        load(Rotor().s, f)
