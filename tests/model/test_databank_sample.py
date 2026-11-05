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
    minimum drag, and drag within 0.0012 over CL 0-0.6 (the bank's drag due to lift is not
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
    assert sec.cdmin == pytest.approx(min(bank(a / 10, 0.3, 2e6)[1] for a in range(-40, 80)),
                                      abs=1e-6)
    a = sec.a0 * 180 / PI + sec.cldmin / (sec.dclda / np.sqrt(1 - sec.mcrit ** 2) * PI / 180)
    rise = bank(a, sec.mcrit, 2e6)[1] - bank(sec.a0 * 180 / PI + sec.cldmin / (
        sec.dclda / np.sqrt(1 - 0.09) * PI / 180), 0.3, 2e6)[1]
    assert rise == pytest.approx(0.002, abs=2e-4)


def test_bank_vs_thickness_estimate():
    """What the 16-series bank says about pyxrotor.sections.naca16 (a record, not a
    requirement on either): its drag-rise Mach number at CLDMIN is 0.02-0.08 above the plain
    thickness estimate 0.85 - 1.2 t at 4-9% thickness, far below the calibrated +0.15-0.20
    (which is a lumped propeller correction, not a section property); its drag due to lift
    is about 8 times the estimate's DCDCL2 = 0.007; its minimum drag 15-50% lower."""
    from pyxrotor import sections
    for t, cld in ((0.04, 0.3), (0.06, 0.3), (0.09, 0.5)):
        f = databank.fit_section(databank.bank16(t, cld), clmax=cld + .85, clmin=cld - .85)
        e = sections.naca16(t, 0.7, cld, mcrit_shift=0.0)
        assert 0.0 < f.mcrit - e.mcrit < 0.1 or t >= 0.09
        assert f.mcrit < sections.mcrit(t) - 0.08
        assert 6 < f.dcdcl2 / e.dcdcl2 < 10
        assert 0.5 < f.cdmin / e.cdmin < 0.9


def test_clark_y_sample(tmp_path):
    """Clark Y, both sample cases (t/c 0.117 at Mach 0-0.4, alpha -1 to 5 deg; t/c 0.30 at
    Mach 0.6 and 0.85, alpha -2 to 10 deg; 10000 ft, chord 0.75 ft): CL and CD within 2e-6
    (the 1986 values sit at or just below the double-precision ones, CL 1.110322 printed
    1.110320). Covers CLTMAX, the thin-section lift and drag, CD11F/CD12F (ALPHAO > 5 at
    11.7%), THKSLP and CD30F/CD35F. The printed Reynolds numbers follow from ATMCON."""
    tol = 2.5e-6
    want = sample("clarky")                  # case, mach, re, alpha, CL, CD
    got = run_fortran(tmp_path, "clarky", ["clarky.f", "atmcon.f"])
    assert np.abs(got[:, 4:6] - want[:, 4:6]).max() < tol
    # 6 printed digits; the 1986 Reynolds numbers are 1-3 units of the 6th higher
    assert np.all(np.abs(got[:, 2] - want[:, 2]) <= 5e-6 * want[:, 2])
    for case, m, re, a, cl, cd in want:
        c = databank.clark_y(a, m, 0.117 if case == 1 else 0.30)
        assert abs(c.cl - cl) < tol and abs(c.cd - cd) < tol, (case, m, a)
    d = databank.cltmax(0.117)
    assert abs(d.cld - 0.505381) < 1e-6 and abs(d.alphlo + 5.186451) < 2e-6
    d = databank.cltmax(0.30)
    assert abs(d.cld - 1.016961) < 1e-6 and abs(d.alphlo + 12.712020) < 2e-6


def test_clark_y_bank_vs_estimate():
    """The Clark Y bank (R&M 2036 data) against pyxrotor.sections.clark_y, a record: at 6-12%
    thickness its minimum drag is within 0.0011 of the estimate (or below it), DCDCL2 within
    0.002 of 0.007, and its drag-rise Mach number at CLDMIN between the plain thickness
    estimate - 0.02 and + 0.08: like the 16-series bank, nowhere near the calibrated +0.15.
    Its minimum drag sits at CL 0.19, not the estimate's 0.35."""
    from pyxrotor import sections
    for t in (0.06, 0.09, 0.12):
        f = databank.fit_section(databank.bank_clark_y(t), clmax=1.3, clmin=-0.5)
        e = sections.clark_y(t, 0.7, mcrit_shift=0.0)
        assert f.cdmin - e.cdmin < 0.0011
        assert abs(f.dcdcl2 - 0.007) < 0.002
        assert -0.02 < f.mcrit - e.mcrit < 0.08
        assert f.mcrit < sections.mcrit(t) - 0.07
        assert abs(f.cldmin - 0.19) < 0.01
