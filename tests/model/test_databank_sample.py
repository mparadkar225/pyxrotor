"""CR-176883's data banks against the report's own sample output.

The transcription (reference/cr176883), compiled as the listings were written (default
REAL*4), must print each sample case's CL and CD to within one unit of the last printed
digit (1e-6): the 1986 machine rounded single precision differently. The double-precision
port (pyxrotor.databank, bit-identical to the REAL*8 build) must agree as closely.
Every printed digit was read from the page images, so a transcription error in a table
or constant shows up here.
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from pyxrotor import databank

CR = Path(__file__).resolve().parents[2] / "reference" / "cr176883"
TOL = 1.5e-6


def sample(name):
    return np.loadtxt(CR / f"sample_{name}.txt", comments="#", ndmin=2)


def run_fortran(tmp_path, name, sources):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not installed")
    exe = tmp_path / name
    subprocess.run(["gfortran", "-std=legacy", "-w", *[str(CR / s) for s in sources],
                    str(CR / f"sample_{name}.f"), "-o", str(exe)], check=True)
    out = subprocess.run([str(exe)], check=True, capture_output=True, text=True).stdout
    return np.array([[float(v) for v in line.split()] for line in out.splitlines() if line.strip()])


def test_naca16_sample(tmp_path):
    want = sample("naca16")                  # mach, alpha, CL, CD
    got = run_fortran(tmp_path, "naca16", ["naca16.f", "atmcon.f"])
    assert np.abs(got[:, 2:4] - want[:, 2:4]).max() < TOL
    for m, a, cl, cd in want:
        atm = databank.atmosphere(0.0)
        c = databank.clcd16(m, a + 2.0, 1, 1.0, 0.04, 0.228150, 0.5, atm.wa, atm.vis, atm.den)
        assert abs(c.cl - cl) < TOL and abs(c.cd - cd) < TOL, (m, a)
    # the printed Reynolds numbers (millions)
    for m, re in ((0.45, 3.27645), (0.70, 5.09669)):
        atm = databank.atmosphere(0.0)
        assert abs(databank.clcd16(m, 2, 1, 1.0, .04, .22815, .5, atm.wa, atm.vis, atm.den).re
                   - re) < 2e-5


def test_naca4_sample(tmp_path):
    want = sample("naca4")                   # alpha, CL, CD
    got = run_fortran(tmp_path, "naca4", ["naca4.f"])
    assert np.abs(got[:, 1:3] - want[:, 1:3]).max() < TOL
    for a, cl, cd in want:
        c = databank.naca4(a, 0.5, 0.12, cld=0.2, altitude_ft=10000.0)
        assert abs(c.cl - cl) < TOL and abs(c.cd - cd) < TOL, a
        assert abs(c.re - 0.274529e7) < 20


def test_naca16_limits():
    """Independent checks of the 16-series bank's structure: lift is linear in the angle
    from zero lift with the Prandtl-Glauert-corrected slope 0.1096 (1 - t)/sqrt(1 - (M+t)^2)
    per degree; drag grows with Mach past the drag-divergence table; the incompressible
    switch (icl 0) returns friction drag only."""
    atm = databank.atmosphere(0.0)
    c1 = databank.clcd16(0.5, 1.0, 1, 1.0, 0.06, 0.3, 0.7, atm.wa, atm.vis, atm.den)
    assert c1.cl == pytest.approx(0.1096 * 0.94 / np.sqrt(1 - 0.56 ** 2))
    c0 = databank.clcd16(0.5, 1.0, 0, 1.0, 0.06, 0.3, 0.7, atm.wa, atm.vis, atm.den)
    assert c0.cd == c0.cd_friction
    cds = [databank.naca16(2.0, m, 0.06, 0.3).cd for m in (0.6, 0.7, 0.8, 0.85)]
    assert np.all(np.diff(cds) > 0)


def test_naca4_limits():
    """4-digit: zero angle gives CL = CLD exactly; minimum drag inside the angle range.
    CLDVSZ's two tables (camber from CLD, CLD from camber) are only roughly inverse:
    CLD 0.4 at 12% thickness comes back as 0.445."""
    c = databank.naca4(0.0, 0.3, 0.12, cld=0.4)
    assert c.cl == 0.4
    z = databank.naca4(2.0, 0.3, 0.12, z=c.z)
    assert z.cld == pytest.approx(0.4, abs=0.05)
    cds = [databank.naca4(a, 0.3, 0.12, cld=0.4).cd for a in np.arange(-8.0, 9.0, 1.0)]
    assert 0 < int(np.argmin(cds)) < len(cds) - 1


def test_fit_section_reproduces_bank():
    """fit_section's XROTOR section follows the bank it was fitted to at the fit point: the
    same lift line (through XROTOR's own CLCDCM, Prandtl-Glauert at the fit Mach), the same
    fitted minimum drag within 0.0004 of the bank's, and drag within 0.0012 over CL 0-0.6 (the bank's drag due to lift is not
    quadratic: its coefficient varies with CL - CLD). Its MCRIT is where the bank's drag at
    CLDMIN, at constant Re, has risen 0.002."""
    from pyxrotor.aero import clcdcm
    from pyxrotor.constants import PI
    bank = databank.bank16(0.06, 0.3)
    sec = databank.fit_section(bank, clmax=1.15, clmin=-0.55, mach=0.3, re=2e6)
    for a in (-1.0, 0.0, 2.0, 4.0):
        cl, cd = bank(a, 0.3, 2e6)
        # CLCDCM(alpha, W, Re, VEL, VSO): local Mach W*VEL/VSO = 0.3
        out = clcdcm(a * PI / 180, 0.3 / 0.5, 2e6, 0.5, 1.0, sec)
        assert out.clift == pytest.approx(cl, abs=2e-3), a
        if 0.0 <= cl <= 0.6:
            assert abs(out.cdrag - cd) < 1.2e-3, a
    # the parabola's minimum sits below the bank's (its drag due to lift isn't quadratic)
    raw = min(bank(a / 10, 0.3, 2e6)[1] for a in range(-40, 80))
    assert 0 < raw - sec.cdmin < 4e-4
    def cd_at_cldmin(m):                      # the bank's lift is linear in the angle
        l1, l2 = bank(-2.0, m, 2e6)[0], bank(4.0, m, 2e6)[0]
        return bank(-2.0 + (sec.cldmin - l1) * 6.0 / (l2 - l1), m, 2e6)[1]
    assert cd_at_cldmin(sec.mcrit) - cd_at_cldmin(0.3) == pytest.approx(0.002, abs=1e-6)


def test_bank_vs_thickness_estimate():
    """What the 16-series bank says about pyxrotor.sections.naca16 (a record, not a
    requirement on either): its drag-rise Mach number at CLDMIN is within -0.08/+0.05 of the
    plain thickness estimate 0.85 - 1.2 t at 3-10% thickness and design CL 0.2-0.7 (highest
    at CL 0.3: the camber-drag fit is a 6-point polynomial in CL), far below the calibrated
    +0.15-0.20 (a lumped propeller correction, not a section property); its drag due to lift
    is about 8 times the estimate's DCDCL2 = 0.007; its minimum drag up to 45% lower."""
    from pyxrotor import sections
    for t in (0.03, 0.06, 0.10):
        for cld in (0.2, 0.3, 0.5, 0.7):
            f = databank.naca16_section(t, 0.9, cld)
            e = sections.naca16(t, 0.9, cld, mcrit_shift=0.0)
            assert -0.08 < f.mcrit - e.mcrit < 0.05, (t, cld)
            assert f.mcrit < sections.mcrit(t) - 0.1
            assert 6 < f.dcdcl2 / e.dcdcl2 < 10
            assert 0.55 < f.cdmin / e.cdmin < 1.0, (t, cld, f.cdmin / e.cdmin)


def test_bank_sections_drop_in():
    """clark_y_section / naca16_section take sections' arguments and keep its stall and
    moment estimates; mcrit_shift adds to the bank's MCRIT (None counts as 0); thick
    sections fall back to the estimate; below CLD 0.2 the 16-series drag terms come from
    CLD 0.2 with the zero-lift angle and CLDMIN following the design CL."""
    from pyxrotor import sections
    for t in (0.06, 0.12):
        b, e = databank.clark_y_section(t, 0.5), sections.clark_y(t, 0.5)
        assert (b.clmax, b.clmin, b.cmcon, b.reref, b.rexp) == (e.clmax, e.clmin, e.cmcon,
                                                                 e.reref, e.rexp)
        assert databank.clark_y_section(t, 0.5, 0.1).mcrit == pytest.approx(b.mcrit + 0.1)
        assert databank.clark_y_section(t, 0.5, None).mcrit == b.mcrit
    assert databank.clark_y_section(0.4, 0.2, 0.05) == sections.clark_y(0.4, 0.2, 0.05)
    assert databank.naca16_section(0.45, 0.2, 0.3) == sections.naca16(0.45, 0.2, 0.3, 0.0)
    lo, c2 = databank.naca16_section(0.06, 0.8, 0.0), databank.naca16_section(0.06, 0.8, 0.2)
    assert (lo.cdmin, lo.dcdcl2, lo.mcrit) == (c2.cdmin, c2.dcdcl2, c2.mcrit)
    assert lo.a0 == 0.0 and lo.cldmin == pytest.approx(c2.cldmin - 0.2)


def test_clark_y_bank_vs_estimate():
    """The Clark Y bank (R&M 2036 data) against pyxrotor.sections.clark_y, a record: at 6-12%
    thickness its minimum drag is within 0.0011 of the estimate (or below it), DCDCL2 within
    0.002 of 0.007, and its drag-rise Mach number at CLDMIN between the plain thickness
    estimate - 0.02 and + 0.08: like the 16-series bank, nowhere near the calibrated +0.15.
    Its minimum drag sits at CL 0.21-0.23, not the estimate's 0.35."""
    from pyxrotor import sections
    for t in (0.06, 0.09, 0.12):
        f = databank.fit_section(databank.bank_clark_y(t), clmax=1.3, clmin=-0.5)
        e = sections.clark_y(t, 0.7, mcrit_shift=0.0)
        assert f.cdmin - e.cdmin < 0.0011
        assert abs(f.dcdcl2 - 0.007) < 0.002
        assert -0.02 < f.mcrit - e.mcrit < 0.08
        assert f.mcrit < sections.mcrit(t) - 0.07
        assert 0.2 < f.cldmin < 0.24
