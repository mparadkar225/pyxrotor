"""The transonic lift loss in the compiled core's CLCDCM (pyxrotor.transonic), by limits and
independent checks:
  * off (State.mlift None, or mlift 0 in aero.clcdcm) is CLCDCM; an MLIFT far above the
    Mach number is CLCDCM bit for bit (the loss is cut off below M_L - 30 W);
  * the C implementation equals the formula written out independently here in NumPy;
  * its CL derivatives match finite differences (they include M_L's dependence on CL);
  * the formula reproduces the lift of CR-176883's Clark Y bank (the R&M 2036 data it was
    fitted to) within 0.035 rms, 0.06 above Mach 0.75, where XROTOR's Prandtl-Glauert lift
    is off by more than 0.5 (on this sample of the fit's range);
  * with a swept blade the loss starts at the normal Mach number;
  * a rotor with the loss on converges where its tips pass M_L.
"""
import math

import numpy as np
import pytest

from pyxrotor import aero, databank, sections, transonic
from pyxrotor.sweep import Sweep
from pyxrotor import _speedups

SEC = sections.clark_y(0.09, 0.7, mcrit_shift=0.15)


def cl(alf, mach, mlift, sec=SEC):
    return aero.clcdcm(alf, mach / 0.5, 0.0, 0.5, 1.0, sec, mlift)


def formula(alf, mach, mlift, sec=SEC):
    """The attached lift of pyxrotor.transonic, before XROTOR's stall model."""
    pg = lambda m: 1.0 / math.sqrt(1.0 - min(m, 0.995) ** 2)
    da = alf - sec.a0
    clpg = sec.dclda * pg(mach) * da
    ml = mlift - transonic.PEN * abs(clpg - sec.cldmin)
    x = transonic.W * np.logaddexp(0.0, (mach - ml) / transonic.W)
    cll = sec.dclda * pg(mach - x) * da
    return math.copysign(max(abs(cll) - (transonic.A * x + transonic.B * x * x), 0.0), cll)


def test_off_and_far_above_are_clcdcm():
    for alf in (-0.05, 0.02, 0.08):
        for mach in (0.3, 0.7, 0.9, 0.98):
            base = cl(alf, mach, 0.0)
            assert cl(alf, mach, 5.0) == base          # bit for bit


@pytest.mark.parametrize("mach", [0.6, 0.75, 0.82, 0.88, 0.95])
@pytest.mark.parametrize("alf", [0.0, 0.03, 0.07])
def test_core_is_the_formula(mach, alf):
    """Below XROTOR's stall limits CLCDCM's CL is the attached lift: check it equals the
    formula (the stall model is XROTOR's and acts on top)."""
    sec = sections.clark_y(0.09, 0.7, mcrit_shift=0.5)       # stall limits out of the way
    sec = aero.AeroSection(**{**sec.__dict__, "clmax": 5.0, "clmin": -5.0})
    got = cl(alf, mach, 0.80, sec).clift
    want = formula(alf, mach, 0.80, sec)
    assert got == pytest.approx(want, rel=1e-6, abs=1e-9)   # (XROTOR's stall smoothing leaks ~1e-9)


@pytest.mark.parametrize("mach", [0.7, 0.8, 0.86, 0.92])
@pytest.mark.parametrize("alf", [-0.03, 0.04, 0.09])
def test_derivatives(mach, alf):
    """(XROTOR's own CL_W leaves out the Mach dependence of its compressible stall limit
    CLMAX = (MCRIT + 0.215 - M)/0.25 + CLDMIN, so that limit is kept out of the way.)"""
    sec = aero.AeroSection(**{**sections.clark_y(0.09, 0.7, mcrit_shift=0.5).__dict__,
                              "clmax": 5.0, "clmin": -5.0})
    h = 1e-6
    r = cl(alf, mach, 0.78, sec)
    d_alf = (cl(alf + h, mach, 0.78, sec).clift - cl(alf - h, mach, 0.78, sec).clift) / (2 * h)
    w = mach / 0.5
    hw = 1e-6
    d_w = (aero.clcdcm(alf, w + hw, 0.0, 0.5, 1.0, sec, 0.78).clift
           - aero.clcdcm(alf, w - hw, 0.0, 0.5, 1.0, sec, 0.78).clift) / (2 * hw)
    assert r.cl_alf == pytest.approx(d_alf, rel=1e-5, abs=1e-7)
    assert r.cl_w == pytest.approx(d_w, rel=1e-5, abs=1e-7)


def test_reproduces_the_clark_y_bank():
    pg = lambda m: 1.0 / math.sqrt(1.0 - min(m, 0.995) ** 2)
    err, err_hi, err_pg = [], [], []
    for t in (0.05, 0.07, 0.09, 0.11, 0.12):
        for a in (2.0, 4.0, 6.0, 8.0):
            c0 = databank.clcdcy(a, t, 0.05).cl / pg(0.05)
            sec = aero.AeroSection(xisect=0.0, a0=0.0, clmax=9.0, clmin=-9.0, dclda=c0,
                                   dclda_stall=0.1, dcl_stall=0.1, cdmin=0.01, cldmin=0.22,
                                   dcdcl2=0.0, cmcon=0.0, reref=1e6, rexp=0.0, mcrit=0.8)
            for m in np.arange(0.3, 0.951, 0.025):
                bank = databank.clcdcy(a, t, m).cl
                e = formula(1.0, m, transonic.mlift(t), sec) - bank
                err.append(e)
                if m > 0.75:
                    err_hi.append(e)
                    err_pg.append(c0 * pg(m) - bank)
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
    assert rms(err) < 0.035 and rms(err_hi) < 0.06 and rms(err_pg) > 0.5


def _station(cos_l, mach, mlift):
    import sys
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from test_sweep import rotol_like
    r = rotol_like()
    s = r.s
    s.vel, s.vso = 100.0, 100.0 / mach
    r.set_transonic_lift([mlift] * s.naero)
    s.sweep = Sweep([0.0, 1.0], [0.0, 0.0])
    s.sweep_cos[:] = cos_l
    return s


def test_swept_loss_at_normal_mach():
    i, alf = 5, 0.05
    s = _station(0.7, 0.95, 0.70)
    swept = _speedups.getclcdcm(s, i, alf, 1.0, 0.0)[0]
    s.mlift = None
    s.conv = False
    swept_off = _speedups.getclcdcm(s, i, alf, 1.0, 0.0)[0]
    assert swept < swept_off                      # Mn = 0.665 is near M_L = 0.70 - pen
    s2 = _station(0.7, 0.70, 0.70)                # Mn = 0.49: well below M_L
    a = _speedups.getclcdcm(s2, i, alf, 1.0, 0.0)[0]
    s2.mlift = None
    b = _speedups.getclcdcm(s2, i, alf, 1.0, 0.0)[0]
    assert a == pytest.approx(b, rel=1e-4)        # (the smooth onset's tail: 1e-5)


def test_rotor_with_loss_converges():
    import sys
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from test_sweep import rotol_like
    ml = [transonic.mlift(t) for t in (0.25, 0.10, 0.07, 0.05)]
    r = rotol_like()
    r.atmosphere(9.0)
    r.set_transonic_lift(ml)
    assert r.solve_rpm(1700.0)                    # helical tip Mach 1.03
    p = r.performance()
    r.set_transonic_lift(None)
    assert r.solve_rpm(1700.0)
    assert p["power"] < r.performance()["power"]  # less lift at the tip, less power
    # far past the local-Mach cap (tip Mach 1.35) 7.55's model has no solution; with the
    # lift loss the lift falls instead of rising, and the solve converges
    r = rotol_like()
    r.atmosphere(9.0)
    assert not r.solve_rpm(2200.0)
    r = rotol_like()
    r.atmosphere(9.0)
    r.set_transonic_lift(ml)
    assert r.solve_rpm(2200.0)
