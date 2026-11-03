"""SETCLD, DESINI, DEITER/DESMIL, DESGEN: the complete compared State must match after each call."""
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor.constants import PI
from pyxrotor.design import deiter, desgen, desini, setcld
from rotors import RADIUS, perturbed_states
from statesync import assert_states_match, state_to_ref
from strategies import finite

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)
from conftest import EXPENSIVE  # noqa: E402


@SETTINGS
@given(s=perturbed_states(), cl1=finite(-0.4, 1.0), cl2=finite(-0.4, 1.0))
def test_setcld(ref_module, s, cl1, cl2):
    state_to_ref(s, ref_module)
    ref_module.call("SETCLD", cl1, cl2)
    setcld(s, cl1, cl2)
    assert_states_match(s, ref_module)


@st.composite
def design_points(draw, allow_vortex=True, clean=False):
    """A state with DESI parameters set around the current operating point (Rotol blade).

    clean: no added inflow or duct velocity, and DESINI always runs (as on a first design), so
    most cases converge; otherwise many iterations diverge, which the port must follow too."""
    s = draw(perturbed_states(duct=not clean, added_inflow=not clean))
    s.nblds = draw(st.sampled_from([2, 3, 4, 5]))
    s.raddes = RADIUS * draw(finite(0.8, 1.2))
    s.veldes = s.vel * draw(finite(0.5, 2.0))
    rpm = s.vel / (s.rad * s.adv * PI / 30.0)
    if draw(st.booleans()):
        s.advdes, s.rpmdes = 0.0, rpm * draw(finite(0.7, 1.3))
    else:
        s.advdes = s.adv * draw(finite(0.6, 1.6))
    s.r0des = s.raddes * draw(finite(0.1, 0.3))
    s.rwdes = s.r0des * draw(finite(0.0, 1.0))
    dim = s.rho * s.vel**2 * s.rad**2
    s.dest = draw(st.booleans())
    s.desp = not s.dest
    s.tddes = abs(s.ttot) * dim * draw(finite(0.3, 2.0))
    s.pddes = abs(s.ptot) * dim * s.vel * draw(finite(0.3, 2.0))
    cl1, cl2 = draw(finite(0.2, 0.9)), draw(finite(0.2, 0.9))
    s.cldes[: s.ii] = (1.0 - s.xi[: s.ii]) * cl1 + s.xi[: s.ii] * cl2
    s.free = draw(st.booleans())
    s.fast = draw(st.booleans())
    s.vrtx = allow_vortex and not s.fast and draw(st.booleans())
    s.duct = False
    s.ldesini = clean or draw(st.booleans())
    s.lrotor = draw(st.booleans())
    s.niterd = draw(st.sampled_from([0, 1, 3, 40, 40, 40]))
    if s.vrtx:
        s.niterd = min(s.niterd, 5)   # vortex-wake designs cost ~0.1 s per iteration
    if draw(st.booleans()):
        s.iaero[: s.ii] = -1          # INIT's state: sections never assigned to stations
    return s


@EXPENSIVE
@given(s=design_points())
def test_desgen(ref_module, s):
    """Mostly harsh inputs: about 2/3 of full-length runs fail to converge or reach NaN."""
    state_to_ref(s, ref_module)
    ref_module.call("REFDESGEN", True, False)
    desgen(s, betz_ok=True)
    assert_states_match(s, ref_module)


@EXPENSIVE
@given(s=design_points(clean=True))
def test_desgen_clean(ref_module, s):
    state_to_ref(s, ref_module)
    ref_module.call("REFDESGEN", True, False)
    desgen(s, betz_ok=True)
    assert_states_match(s, ref_module)


@EXPENSIVE
@given(s=design_points(allow_vortex=False), icase=st.sampled_from([2, 3]))
def test_desini_deiter(ref_module, s, icase):
    """DESINI and DEITER separately, from the setup DESGEN would do."""
    s.niterd = min(s.niterd, 3)
    desgen_setup = s.copy()
    desgen_setup.niterd = 0
    desgen_setup.lrotor = True
    desgen_setup.ldesini = False
    desgen(desgen_setup, betz_ok=True)        # grid, CL, EFFINV and ADW, one DEITER pass
    desgen_setup.niterd = s.niterd
    for name, fn in (("DESINI", lambda st_: desini(st_)), ("DEITER", lambda st_: deiter(st_, icase))):
        state_to_ref(desgen_setup, ref_module)
        if name == "DESINI":
            ref_module.call("DESINI")
        else:
            ref_module.call("DEITER", icase)
        fn(desgen_setup)
        assert_states_match(desgen_setup, ref_module)


@pytest.mark.parametrize("betz_ok", [False, True])
def test_desgen_betz_limit(ref, betz_ok):
    """Windmill thrust past the Betz limit: returns without designing unless confirmed."""
    import rotors
    s = rotors.solved_state()
    s.nblds, s.raddes, s.veldes, s.advdes = 3, RADIUS, 50.0, 0.8
    s.r0des, s.rwdes = 0.2 * RADIUS, 0.0
    s.dest, s.desp = True, False
    s.tddes = -0.6 * PI * s.rho * 50.0**2 * RADIUS**2
    s.cldes[: s.ii] = 0.5
    s.niterd = 40
    state_to_ref(s, ref)
    ref.call("REFDESGEN", betz_ok, False)
    assert desgen(s, betz_ok=betz_ok) is betz_ok
    assert_states_match(s, ref)


def test_desgen_stall_limit(ref):
    import rotors
    s = rotors.solved_state()
    s.cldes[: s.ii] = 0.5
    s.cldes[3] = 5.0
    s.raddes, s.veldes, s.advdes, s.dest, s.tddes = RADIUS, 50.0, 0.8, True, 500.0
    state_to_ref(s, ref)
    ref.call("REFDESGEN", True, False)
    assert desgen(s) is False
    assert_states_match(s, ref)


# ------------------------------------------------------------------ DESI commands (Rotor)
from pyxrotor.aero import AeroSection  # noqa: E402
from pyxrotor.design import BetzLimitError  # noqa: E402
from pyxrotor.rotor import Rotor  # noqa: E402
from reference.ref755 import Ref755  # noqa: E402
from rotors import SECTIONS  # noqa: E402
from strategies import load_sections, section_sets  # noqa: E402

ROTOL_SECTIONS = [
    AeroSection(xisect=x, a0=a0 * PI / 180, dclda=dcla, clmax=clmax, clmin=clmin,
                dclda_stall=0.1, dcl_stall=0.1, cdmin=cdmin, cldmin=cldmin, dcdcl2=dcdcl2,
                cmcon=-0.1, mcrit=mcrit, reref=2.0e6, rexp=-0.15)
    for x, a0, dcla, clmax, clmin, cdmin, cldmin, dcdcl2, mcrit in SECTIONS]


@st.composite
def design_inputs(draw):
    rad = draw(finite(0.8, 2.2))
    d = dict(nblds=draw(st.integers(2, 6)), tip_radius=rad, hub_radius=rad * draw(finite(0.1, 0.25)),
             vel=draw(finite(30.0, 150.0)), cl=draw(st.one_of(st.just(0.0), finite(0.3, 0.8))))
    d["wake_hub_radius"] = d["hub_radius"] * draw(st.sampled_from([0.0, 0.5, 1.0]))
    if draw(st.booleans()):
        d["adv"] = draw(finite(0.2, 0.9))
        omega_r = d["vel"] / d["adv"]
    else:
        d["rpm"] = draw(finite(800.0, 3000.0))
        omega_r = d["rpm"] * PI / 30.0 * rad
    q = 0.5 * 1.226 * (d["vel"] ** 2 + omega_r ** 2) * PI * rad * rad      # rough disk scale
    if draw(st.booleans()):
        d["thrust"] = q * draw(finite(0.005, 0.1))
    else:
        d["power"] = q * d["vel"] * draw(finite(0.005, 0.15))
    return d


@st.composite
def edits(draw, base):
    """A few EDIT keywords, in order."""
    out = {}
    for key in draw(st.lists(st.sampled_from(
            ["nblds", "tip_radius", "vel", "adv", "rpm", "thrust", "power", "cl", "cl_root_tip"]),
            min_size=1, max_size=3, unique=True)):
        if key == "nblds":
            out[key] = draw(st.integers(2, 6))
        elif key == "tip_radius":
            out[key] = base["tip_radius"] * draw(finite(0.9, 1.1))
        elif key == "vel":
            out[key] = base["vel"] * draw(finite(0.8, 1.2))
        elif key == "adv":
            out[key] = draw(finite(0.2, 0.9))
        elif key == "rpm":
            out[key] = draw(finite(800.0, 3000.0))
        elif key in ("thrust", "power"):
            b = base.get(key) or base.get("power", 0.0) / base["vel"] or base.get("thrust", 0.0) * base["vel"]
            out[key] = abs(b) * draw(finite(0.7, 1.3))
        elif key == "cl":
            out[key] = draw(finite(0.3, 0.8))
        else:
            out[key] = (draw(finite(0.3, 0.9)), draw(finite(0.3, 0.9)))
    return out


def ref_redesign(ref, changes, fix):
    """EDIT keywords on the reference's COMMON, then <Return> (LDESINI, DESGEN)."""
    for key, v in changes.items():
        if key == "nblds":
            ref["NBLDS"] = v
        elif key == "tip_radius":
            ref["RADDES"] = v
        elif key == "vel":
            ref["VELDES"] = v
        elif key == "adv":
            ref["ADVDES"], ref["RPMDES"] = v, 0.0
        elif key == "rpm":
            ref["RPMDES"], ref["ADVDES"] = v, 0.0
        elif key == "thrust":
            ref["TDDES"], ref["DEST"], ref["DESP"], ref["PDDES"] = v, True, False, 0.0
        elif key == "power":
            ref["PDDES"], ref["DEST"], ref["DESP"], ref["TDDES"] = v, False, True, 0.0
        elif key == "cl":
            ref["CLDES0"] = v
            ref.call("SETCLD", v, v)
        elif key == "cl_root_tip":
            ref.call("SETCLD", *v)
            ref["CLDES0"] = 0.5 * (v[0] + v[1])
    ref["LDESINI"] = True
    ref.call("REFDESGEN", False, fix)


def both(fn_rotor, fn_ref):
    try:
        fn_rotor()
    except (BetzLimitError, ValueError):
        pass
    fn_ref()


@st.composite
def design_sessions(draw):
    secs = draw(st.one_of(st.just(ROTOL_SECTIONS), st.just(ROTOL_SECTIONS[:1]), section_sets()))
    inputs = draw(design_inputs())
    model = draw(st.sampled_from(["pot", "pot", "grad", "vrtx"]))
    free = draw(st.booleans())
    ii = draw(st.sampled_from([None, None, 15, 41]))
    steps = draw(st.lists(st.one_of(edits(inputs).map(lambda e: ("edit", e)),
                                    st.integers(8, 50).map(lambda n: ("n", n)),
                                    finite(0.2, 0.9).map(lambda a: ("adva", a))), max_size=3))
    return secs, inputs, model, free, ii, steps


def run_design_session(session, fix):
    secs, inputs, model, free, ii, steps = session
    ref, rotor = Ref755(), Rotor(refresh_design_sections=fix)
    s = rotor.s
    if ii is not None:
        ref["II"], ref["IINF"] = ii, ii + ii // 2
        s.ii, s.iinf = ii, ii + ii // 2
    load_sections(ref, secs)
    ref.call("SETIAERO")
    rotor.set_sections(secs)
    rotor.design_formulation(model, toggle_free_wake=not free)
    ref["VRTX"], ref["FAST"], ref["FREE"] = s.vrtx, s.fast, s.free
    if model == "vrtx":
        s.niterd = 15                          # keep vortex-wake designs affordable
        ref["NITERD"] = 15
    d = inputs
    both(lambda: rotor.design(betz_ok=False, **d),
         lambda: ref.call("REFDINPU", d["nblds"], d["tip_radius"], d["hub_radius"],
                          d["wake_hub_radius"], d["vel"], d.get("adv", 0.0), d.get("rpm", 0.0),
                          d.get("thrust", 0.0), d.get("power", 0.0), d["cl"], False, fix))
    assert_states_match(s, ref)
    for kind, arg in steps:
        if kind == "edit":
            both(lambda: rotor.redesign(**arg), lambda: ref_redesign(ref, arg, fix))
        elif kind == "n":
            ref.call("REFDESN", arg)
            rotor.set_stations(arg)
        else:
            ref.call("REFADVA", arg)
            rotor.solve_advance_ratio(arg)
        assert_states_match(s, ref)
    return rotor


@EXPENSIVE
@given(session=design_sessions())
def test_design_sessions(session):
    """INPU, EDIT, N and an OPER solve, as 7.55 does them (stale section indices included)."""
    run_design_session(session, fix=False)


@EXPENSIVE
@given(session=design_sessions())
def test_design_sessions_fixed_sections(session):
    """The same with refresh_design_sections=True, against DESGEN plus SETIAERO after SETX."""
    run_design_session(session, fix=True)


def p51d_design(fix):
    """The P-51D exercise: 4 blades, 1.702 m, 112 m/s at 7.62 km, 1437 rpm, 1118.55 kW, CL 0.5,
    three sections defined before the first design."""
    from dataclasses import replace
    d2r = PI / 180
    s1 = AeroSection(xisect=0.0, a0=-2.0 * d2r, clmax=1.4, clmin=-0.6, cdmin=0.010, cldmin=0.40,
                     dcdcl2=0.006, reref=3.0e6, rexp=-0.15, mcrit=0.65)
    s2 = replace(s1, xisect=0.7, a0=-1.5 * d2r, cdmin=0.0075, cldmin=0.30, mcrit=0.80)
    r = Rotor(refresh_design_sections=fix)
    r.atmosphere(7.62)
    r.set_sections([s1, s2, replace(s2, xisect=1.0)])
    assert r.design(nblds=4, tip_radius=1.702, hub_radius=0.20, wake_hub_radius=0.15, vel=112.0,
                    rpm=1437.0, power=1118550.0, cl=0.5)
    return r


def test_p51d_stale_sections():
    """7.55 (interactive run, design.out) reports 85.27% / induced 85.96%, with negative CD on
    the outer stations. Re-analysing that blade with the right sections (as after SAVE/LOAD)
    gives 79.80% at 1.177 MW; designing with the right sections gives 80.26%."""
    from pyxrotor.aero import setiaero
    from pyxrotor.rotor import sections_of
    bug = p51d_design(fix=False)
    p = bug.performance()
    assert (round(p["efficiency"], 4), round(p["induced_efficiency"], 4)) == (0.8527, 0.8596)
    assert bug.stations()["cd"].min() < 0.0
    s = bug.s
    s.iaero[: s.ii] = setiaero(s.xi[: s.ii], sections_of(s))
    assert bug.solve_rpm(1437.0)
    p = bug.performance()
    assert (round(p["efficiency"], 4), round(p["power"], -3)) == (0.7980, 1177000.0)

    fixed = p51d_design(fix=True)
    p = fixed.performance()
    assert round(p["efficiency"], 4) == 0.8026 and fixed.stations()["cd"].min() > 0.0
    assert abs(p["power"] / 1118550.0 - 1.0) < 1e-5
