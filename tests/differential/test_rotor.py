"""Start-up, ATMO, ARBI and the OPER commands end to end against the reference."""
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyxrotor.aero import AeroSection
from pyxrotor.constants import PI
from pyxrotor.rotor import Rotor, atmo, startup
from reference.ref755 import Ref755
from rotors import RADIUS, SECTIONS, geometry
from statesync import assert_states_match
from strategies import finite, section_sets

pytestmark = pytest.mark.differential
SETTINGS = settings(deadline=None)
from conftest import EXPENSIVE  # noqa: E402


def test_startup_matches_reference(ref):
    assert_states_match(startup(), ref)


@SETTINGS
@given(alt=st.one_of(finite(-2.0, 90.0), st.sampled_from([-1.0, 0.0, 11.0, 40.0, 75.0, 80.0])))
def test_atmo(ref_module, alt):
    _, vso, rho, rmu = ref_module.call("ATMO", alt, 0.0, 0.0, 0.0)
    assert atmo(alt) == (vso.item(), rho.item(), rmu.item())


@st.composite
def blades(draw):
    n = draw(st.integers(3, 15))
    hub = draw(finite(0.05, 0.35))
    tip = draw(finite(0.3, 3.0))
    # ARBI needs XI0 <= r/R, with XI0 = ROOT/RAD = (hub*tip)/tip, which can round to one ulp
    # above hub; REFARBI STOPs (exiting the test process) on a station below it
    xi0 = (hub * tip) / tip
    r = np.sort(draw(st.lists(finite(xi0, 1.0), min_size=n, max_size=n, unique=True)))
    c = np.array(draw(st.lists(finite(0.02, 0.3), min_size=n, max_size=n)))
    b = np.array(draw(st.lists(finite(5.0, 70.0), min_size=n, max_size=n)))
    return dict(nblds=draw(st.integers(1, 6)), vel=draw(finite(5.0, 150.0)),
                tip_radius=tip, hub=hub, r=r, c=c, b=np.sort(b)[::-1].copy())


def setup_both(blade, secs, ii=None, ixspac=2):
    ref = Ref755()
    rotor = Rotor()
    for obj in (ref, rotor.s):
        if ii is not None:
            if obj is ref:
                ref["II"], ref["IINF"] = ii, ii + ii // 2
            else:
                obj.ii, obj.iinf = ii, ii + ii // 2
    ref["IXSPAC"] = rotor.s.ixspac = ixspac
    from strategies import load_sections
    load_sections(ref, secs)
    rotor.set_sections(secs)
    ref.call("SETIAERO")
    rad = blade["tip_radius"]
    ref.call("REFARBI", blade["nblds"], blade["vel"], rad, blade["hub"] * rad, len(blade["r"]),
             blade["r"], blade["c"], blade["b"])
    rotor.arbi(blade["nblds"], blade["vel"], rad, blade["hub"] * rad, blade["r"], blade["c"], blade["b"])
    return ref, rotor


@EXPENSIVE
@given(blade=blades(), secs=section_sets(), ii=st.sampled_from([None, 12, 45]),
       ixspac=st.sampled_from([1, 2]))
def test_arbi(blade, secs, ii, ixspac):
    ref, rotor = setup_both(blade, secs, ii, ixspac)
    assert_states_match(rotor.s, ref)


def rotol_both(theta07):
    geo = geometry()
    secs = [AeroSection(xisect=x, a0=a0 * PI / 180, dclda=dcla, clmax=clmax, clmin=clmin,
                        dclda_stall=0.1, dcl_stall=0.1, cdmin=cdmin, cldmin=cldmin, dcdcl2=dcdcl2,
                        cmcon=-0.1, mcrit=mcrit, reref=2.0e6, rexp=-0.15)
            for x, a0, dcla, clmax, clmin, cdmin, cldmin, dcdcl2, mcrit in SECTIONS]
    blade = dict(nblds=2, vel=50.0, tip_radius=RADIUS, hub=0.2,
                 r=np.array([g["r"] for g in geo]), c=np.array([g["cR"] for g in geo]),
                 b=np.array([theta07 + g["dth"] for g in geo]))
    ref, rotor = setup_both(blade, secs)
    ref.call("REFATMO", 0.0)
    rotor.atmosphere(0.0)
    return ref, rotor


@st.composite
def oper_sessions(draw):
    """A short OPER session: a few commands in sequence, as a user would type them."""
    cmds = []
    for _ in range(draw(st.integers(1, 4))):
        kind = draw(st.sampled_from(["velo", "angl", "adva", "rpm", "thru", "torq", "powe", "flags"]))
        if kind == "velo":
            cmds.append(("velo", draw(finite(10.0, 80.0))))
        elif kind == "angl":
            cmds.append(("angl", draw(finite(-4.0, 4.0))))
        elif kind == "adva":
            cmds.append(("adva", draw(finite(0.05, 0.4))))
        elif kind == "rpm":
            cmds.append(("rpm", draw(finite(500.0, 1500.0))))
        elif kind == "flags":
            cmds.append(("flags", draw(st.booleans()), draw(st.booleans()), draw(st.booleans())))
        else:
            cmds.append((kind, draw(finite(0.6, 1.4)), draw(st.sampled_from([1, 2])),
                         draw(st.one_of(st.just(0.0), finite(600.0, 1400.0)))))
    return draw(st.sampled_from([15.0, 22.0, 27.0])), cmds


def run_session(ref, rotor, cmds):
    s = rotor.s
    for c in cmds:
        if c[0] == "velo":
            ref.call("REFVELO", c[1]); rotor.velocity(c[1])
        elif c[0] == "angl":
            ref.call("REFANGL", c[1]); rotor.change_pitch(c[1])
        elif c[0] == "adva":
            ref.call("REFADVA", c[1]); rotor.solve_advance_ratio(c[1])
        elif c[0] == "rpm":
            ref.call("REFRPM", c[1]); rotor.solve_rpm(c[1])
        elif c[0] == "flags":
            for name, val in zip(("FREE", "FAST", "LOPRINI"), c[1:]):
                ref[name] = val
                setattr(s, name.lower(), val)
        else:
            ispec = {"thru": 1, "torq": 2, "powe": 3}[c[0]]
            scale, ifix, rpm = c[1], c[2], c[3]
            dim = s.rho * s.vel * s.vel * s.rad * s.rad
            # targets: a multiple of the current value, or a nominal one before any solve
            base = {1: s.ttot * dim or 300.0, 2: s.qtot * dim * s.rad or 150.0,
                    3: s.ptot * dim * s.vel or 3.0e4}[ispec]
            spec = base * scale
            ref.call("REFTPQ", ispec, spec, ifix, rpm)
            fix = "pitch" if ifix == 1 else "rpm"
            getattr(rotor, {1: "solve_thrust", 2: "solve_torque", 3: "solve_power"}[ispec])(
                spec, fix=fix, rpm=rpm if rpm > 0 else None)
        assert_states_match(s, ref)


@EXPENSIVE
@given(session=oper_sessions())
def test_oper_sessions(session):
    theta07, cmds = session
    ref, rotor = rotol_both(theta07)
    assert_states_match(rotor.s, ref)
    run_session(ref, rotor, cmds)
