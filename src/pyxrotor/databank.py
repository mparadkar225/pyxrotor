"""NASA CR-176883's airfoil data banks: section CL and CD of the NACA 16-series, the NACA
4-digit family (maximum camber at mid-chord) and the Clark Y, as functions of thickness,
camber, angle of attack, Mach and Reynolds number.

Korkan, Camba & Morris, "Aerodynamic data banks for Clark-Y, NACA 4-digit and NACA
16-series airfoil families", Texas A&M, NASA CR-176883, 1986 (NTRS 19860021221).

Not XROTOR: a port of the report's Fortran listings, which are transcribed in
reference/cr176883/ (the oracle; tests/differential/test_databank.py compares them bit for
bit, built with REAL*8 like the XROTOR reference, and tests/model/test_databank_sample.py
reproduces the report's own sample output with the REAL*4 build). Operation order,
integer powers (gfortran's multiplication chains, `_powi`) and libm calls follow the
compiled Fortran. The limitation warnings the programs print are not reproduced; the
ranges they check are listed with each function.

    atmosphere(altitude_ft)             ATMCON: the banks' standard atmosphere (English units)
    naca16(alpha_deg, mach, tc, cld)    16-series, at the main program's angle convention
    clcd16(...)                         its CLCD, angle from zero lift
    naca4(alpha_deg, mach, tc, ...)     4-digit, angle from the design attitude
    clark_y(alpha_deg, mach, tc)        Clark Y, at the main program's angle convention
    clcdcy(alphao, t, hm)               its CLCD, angle from zero lift
    cltmax(t, ncld)                     its design CL from t/c (or t/c from CL), ALPHLO
    bank16, bank4, bank_clark_y         the banks as cl_cd(alpha_deg, mach, re)
    fit_section(cl_cd, ...)             XROTOR section parameters fitted to a bank
    clark_y_section, naca16_section     drop-in for sections.clark_y / naca16, from the banks
"""
from __future__ import annotations

import dataclasses
import functools
import math
from typing import NamedTuple

from . import fmath


def _powi(x: float, n: int) -> float:
    """x**n for a constant integer n as gfortran -O1 compiles it (GCC's powi table):
    x2 = x*x, x3 = x2*x, x4 = x2*x2, x5 = x3*x2, x6 = x3*x3."""
    if n == 1:
        return x
    x2 = x * x
    if n == 2:
        return x2
    x3 = x2 * x
    if n == 3:
        return x3
    if n == 4:
        return x2 * x2
    if n == 5:
        return x3 * x2
    if n == 6:
        return x3 * x3
    raise ValueError(n)


def _lagrange(xs, ys, z):
    """The 16-series CLCD's inline Lagrange interpolation over all points."""
    total = 0.0
    for i, (xi, s) in enumerate(zip(xs, ys)):
        for j, xj in enumerate(xs):
            if j != i:
                s = fmath.div(s * (z - xj), xi - xj)
        total = total + s
    return total


# ----------------------------------------------------------------------------- ATMCON
class Atmosphere(NamedTuple):
    temp: float      # deg R
    pres: float      # lb/ft^2
    den: float       # slug/ft^3
    vis: float       # lb s/ft^2
    wa: float        # speed of sound, ft/s


def atmosphere(altitude_ft: float) -> Atmosphere:
    """ATMCON (all three banks; the 4-digit one's arguments differ, its arithmetic is the
    same). Troposphere lapse rate, isothermal above 11 km, Sutherland viscosity."""
    alt = float(altitude_ft) * .3048
    slope = -.0065
    temp = 288.16 + slope * alt
    e = 2.718281828
    tsea = 288.16
    psea = 101325.
    pres1 = 22700.
    grav = 9.82
    r = 288.
    power = grav / (slope * r)
    pres = psea * fmath.pow(1. / (temp / tsea), power)
    if not alt <= 11000.:
        pres = pres1 * fmath.pow(e, -grav / (r * 216.66) * (alt - 11000.))
        temp = 216.66
    wa = fmath.sqrt(1.4 * r * temp) / .3048
    temp = (temp - 273.15) * 9. / 5. + 491.67
    vis = .35 * fmath.pow(temp / 492., 1.5) * (690. / (temp + 198.)) * 1.E-6
    pres = pres * 2116.2 / psea
    den = pres / (1718. * temp)
    return Atmosphere(temp, pres, den, vis, wa)


# ----------------------------------------------------------------------------- 16-series
_A = (1.325276, .1024011, .2890272, .5288957, 1.177606, 2.125449)
_B = (-5.133061, -.0976595, -.2461832, -.0205697, -.2953677, .1974199)
_C = (8.412967, -.2752136, -.8075892, -2.787452, -5.423342, -11.9608)
_D = (-7.325109, .5574539, 1.520081, 4.335364, 8.694619, 18.21422)
_E = (3.522405, -.3924382, -1.018046, -2.761126, -5.583289, -11.50479)
_F = (-.8777567, .1226128, .3042805, .8093843, 1.64263, 3.357388)
_G = (.087806, -.0140459, -.0334875, -.0889465, -.1814412, -.3697504)
_TB = (.031, .0625, .1, .1375, .175, .2125, .25, .2875, .325, .3625, .4)
_CDO = (0.0, .00025, .00065, .0013, .0022, .0033, .00465, .0061, .0078, .01, .013)
_DMACH = (.9, .84, .79, .75, .72, .7, .68, .665, .65, .6375, .62)
_AA = (-5.834527753, -1.45916514, 8.3327673, 2.6038473, -.10428842, .208297374)
_BB = (37.09086074, 9.58857686, -51.4547602, -16.2479854, .78201755, -1.2497732)
_CC = (-94.01887327, -25.1902339, 126.241026, 40.198068, -2.2988001, 2.94734747)
_DD = (118.8401952, 33.12056567, -153.736709, -49.206188, 3.3320883, -3.39928973)
_EE = (-74.98024889, -21.858250, 92.8384573, 29.7250413, -2.3960093, 1.903809351)
_FF = (18.94459392, 5.83350647, -22.1932814, -7.0552833, .689991999, -.408891283)
_CLIT = tuple(float(i) / 10. + .1 for i in range(1, 7))
_TB1 = tuple(float(i) / 100. for i in range(1, 7))


class Coefficients16(NamedTuple):
    cl: float
    cd: float
    re: float            # Reynolds number, millions (XRN on return)
    cd_friction: float   # CDFRIC
    cd_lift: float       # CDLIFT
    cd_camber: float     # CDCAMB
    cd_zero: float       # CDZERO (thickness and compressibility)


def clcd16(xm, alpha, icl, chord, xtc, xcli, xr, wa, vis, den) -> Coefficients16:
    """The 16-series bank's SUBROUTINE CLCD (report pp. 143-146; CLCD16 in the transcription).

    xm Mach number; alpha (deg) from zero lift; icl 0 gives the incompressible lift
    slope and friction drag only; chord (ft); xtc thickness ratio; xcli design CL; xr
    the station's r/R (sets the laminar fraction PLAM of the friction law); wa, vis, den
    from `atmosphere`. ICD is unused in the listing and is left out.

    Ranges the main program warns outside of: t/c 0.02-0.30 and Mach 0.3-0.75 for the
    moment, CLD 0-0.5; the drag-rise data cover t/c to 0.4."""
    xm, alpha, chord, xtc, xcli, xr = map(float, (xm, alpha, chord, xtc, xcli, xr))
    wa, vis, den = float(wa), float(vis), float(den)
    plam = (33. + 1. / 3.) * (xr * xr) + (8. + 1. / 3.) * xr + 17.0
    if not xm < 1.0:
        dclda = fmath.div(.0698, fmath.sqrt(xm * xm - 1))
    elif icl == 0:
        dclda = fmath.div(.1096 * (1 - xtc), fmath.sqrt(1 - xtc * xtc))
    else:
        dclda = fmath.div(.1096, fmath.sqrt(1 - xm * xm))
        if not (xm + xtc) >= 1.0:
            s = xm + xtc
            dclda = fmath.div(.1096 * (1. - xtc), fmath.sqrt(1 - s * s))
    cl = alpha * dclda
    xrn = fmath.div(den * xm * wa * chord, vis)
    delta = cl - xcli
    xmxt = fmath.div(xm, .9 - xtc)
    d2 = delta * delta
    delcd = (xmxt * xmxt - xmxt - .1765) * (-.4397 * d2 + .1173 * delta - .0938)
    testm = (delcd - .03) * (.9 - .875) / (.03 - .08) + .9
    if xmxt > testm:
        slope = (-117.218415 * _powi(delta, 6) + 321.890357 * _powi(delta, 5)
                 - 332.087766 * _powi(delta, 4) + 151.8614153 * _powi(delta, 3)
                 - 24.597998 * d2 - 2.6451615 * delta + 1.48)
        delcd = slope * (xmxt - 1.011111) + .0911111
    if xm >= 1.0:
        delcd = fmath.sqrt(xm * xm - 1) / 4.
        delta = cl
    cdlift = delcd * (delta * delta)
    if cdlift < 0.0:
        cdlift = 0.0
    sq = fmath.sqrt(xrn)
    lg = _log10(xrn)
    cdfric = (fmath.div(2.65, sq) + (100 - plam) / 100.
              * (fmath.div(.944, fmath.pow(lg, 2.6)) - fmath.div(2.65, sq)))
    cdcamb = 0.0
    if not xmxt > 1.0:
        cdd = [_A[i] * _powi(xmxt, 6) + _B[i] * _powi(xmxt, 5) + _C[i] * _powi(xmxt, 4)
               + _D[i] * _powi(xmxt, 3) + _E[i] * (xmxt * xmxt) + _F[i] * xmxt + _G[i]
               for i in range(6)]
        cdcamb = _lagrange(_CLIT, cdd, xcli)
        if not cdcamb > 0.0:
            cdcamb = 0.0
    if xm >= .95:
        cdo1 = [_AA[i] * _powi(xm, 5) + _BB[i] * _powi(xm, 4) + _CC[i] * _powi(xm, 3)
                + _DD[i] * (xm * xm) + _EE[i] * xm + _FF[i] for i in range(6)]
        cdzero = _lagrange(_TB1, cdo1, xtc)
        if cdzero < 0.0:
            cdzero = 0.0
    elif xtc <= .03:
        cdzero = 0.0
    else:
        cdzero = _lagrange(_TB, _CDO, xtc)
        divmac = _lagrange(_TB, _DMACH, xtc)
        if not xm <= divmac:
            xc = .7225 + (xm - divmac)
            ychang = cdzero - .0023845
            cdzero = (-37732.79544 * _powi(xc, 5) + 143059.1078 * _powi(xc, 4)
                      - 216835.3262 * _powi(xc, 3) + 164241.3211 * (xc * xc)
                      - 62170.09813 * xc + 9408.610638 + ychang)
        if cdzero < 0.0:
            cdzero = 0.0
    cd = cdlift + cdfric + cdcamb + cdzero
    if icl == 0:
        cd = cdfric
    return Coefficients16(cl, cd, xrn / 1000000., cdfric, cdlift, cdcamb, cdzero)


def naca16_alpha0(tc: float, cld: float) -> float:
    """Zero-lift angle (deg) of the 16-series main program: -CLD / (0.1096 sqrt((1+t)/(1-t)))."""
    tc, cld = float(tc), float(cld)
    return -cld / (0.1096 * math.sqrt((1 + tc) / (1 - tc)))


def naca16(alpha_deg, mach, tc, cld, *, chord_ft=1.0, altitude_ft=0.0, x=0.5,
           icl=1) -> Coefficients16:
    """16-series section at geometric angle alpha_deg (the main program's ANGLE): CLCD at
    ANGLE - ALPHAO. The Reynolds number follows from Mach, chord and altitude."""
    atm = atmosphere(altitude_ft)
    ang = float(alpha_deg) - naca16_alpha0(tc, cld)
    return clcd16(mach, ang, icl, chord_ft, tc, cld, x, atm.wa, atm.vis, atm.den)


# ----------------------------------------------------------------------------- 4-digit
def linear(xin, yin, xout):
    """SUBROUTINE LINEAR (4-digit bank): straight line through the bracketing pair, the
    end pair outside the table."""
    nu = len(xin)
    nm1 = nu - 1
    f = 1
    if not xout < xin[0]:
        f = nm1
        if not xout > xin[nu - 1]:
            for ii in range(2, nm1 + 1):
                f = ii - 1
                if xout <= xin[ii - 1]:
                    break
            else:
                f = nm1
    y1, x1, y2, x2 = yin[f - 1], xin[f - 1], yin[f], xin[f]
    slp = fmath.div(y2 - y1, x2 - x1)
    b2 = y1 - x1 * slp
    return slp * xout + b2


_ZI = (0.02, 0.04, 0.06, 0.08)
_CLDI = (0.2, 0.4, 0.6, 0.8)
_TCI = (0.04, 0.06, 0.09, 0.12, 0.18)
_CLDATA = ((.225537, .235576, .246819, .256079, .265394),
           (.455984, .48389, .547405, .52539, .489006),
           (.672912, .707067, .786416, .691206, .6380),
           (.892327, .899246, .906358, .930514, .843705))
_ZDATA = ((.0195, .0186, .0175, .0168, .0160),
          (.0370, .0360, .0350, .0340, .0330),
          (.0550, .0540, .0530, .0520, .0530),
          (.0730, .0720, .0710, .0700, .0760))


def cldvsz(camp, tc, cld=0.0, z=0.0):
    """CLDVSZ: design CL from camber z (camp 1) or camber from CLD (camp 0). Returns (cld, z)."""
    tc = float(tc)
    if camp == 1:
        cldta = [linear(_TCI, c, tc) for c in _CLDATA]
        return linear(_ZI, cldta, float(z)), float(z)
    zadata = [linear(_TCI, c, tc) for c in _ZDATA]
    return float(cld), linear(_CLDI, zadata, float(cld))


_Z1 = (0.0, 0.02, 0.04, 0.06, 0.08)
_TC1 = (0.03, .06, .09, .12, .15, .18, .21)
_LDATA = ((1.4, 1.3, 1.4, 1.4, 1.4, 1.4, 1.4),
          (1.5, 1.4, 1.5, 1.5, 1.4, 1.4, 1.3),
          (1.65, 1.6, 1.7, 1.7, 1.65, 1.6, 1.55),
          (1.8, 1.9, 1.9, 1.8, 1.7, 1.8, 1.65),
          (2.0, 2.1, 2.1, 2.0, 1.9, 1.8, 1.7))


def maxcl(z, tc):
    """MAXCL: CLmax (McCormick) for camber z, thickness tc."""
    return linear(_Z1, [linear(_TC1, c, float(tc)) for c in _LDATA], float(z))


_ZP = (0.02, 0.04, 0.06, 0.08)
_TCP = (0.04, .06, .09, .12, .18)
_X2 = (-0.3, -.2, -.1, 0., .1, .2, .3, .4, .5)
_X4 = (-0.5, -.4, -.3, -.2, -.1, 0., .1, .2, .3, .4, .5)
_X6 = (-0.7, -.6, -.5, -.4, -.3, -.2, -.1, 0., .1, .2, .3, .4)
# CR-176883: the listing's X8 continuation line is "O1,.2,.3,.4,.5" (column 6 holds the
# continuation mark), so X8(13) = 1. as compiled; .1 was meant (fix_x8).
_X8 = (-1.1, -1.0, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0., 1., .2, .3, .4, .5)
_X8_FIXED = _X8[:12] + (.1,) + _X8[13:]
_DCD = (
    ((.00065, .0003, 0., 0., .0003, .00065, .0011, .0017, .0025),
     (.0009, .00035, .0001, 0., .00025, .0008, .0018, .0032, .0046),
     (.0014, .0006, .0002, .0002, .0006, .0012, .0026, .0042, .0058),
     (.0014, .0008, .0003, .0002, .0001, 0., .0015, .0036, .0058),
     (.0021, .0018, .0013, 0., .001, .0015, .0019, .0023, .0041)),
    ((.0008, .00045, .0002, .0001, 0., .00015, .0004, .0008, .0015, .0022, .0030),
     (.0016, .0013, .001, .00065, .0002, 0., .0006, .0017, .003, .0044, .0059),
     (.001, .0007, .0005, .0002, 0., 0., 0., .0008, .0021, .0036, .0052),
     (.0012, .0007, .0004, .0002, .0005, 0., .0001, .0003, .0022, .0040, .0066),
     (.0042, .0037, .0032, .0026, .002, .0008, .0001, .0001, .0018, .0032, .0051)),
    ((.0013, .0011, .001, .0011, .0013, .0016, .0021, .0028, .0037, .0047, .0057, .0067),
     (.0015, .0014, .00135, .0013, .0013, .0013, .00135, .0018, .0027, .0042, .0062, .008),
     (.0005, .0004, .00045, .0006, .0007, .0003, .0003, .0009, .001, .0018, .0032, .005),
     (.0020, .0018, .0016, .0013, .001, .0008, .0005, .00045, .001, .0028, .0039, .0053),
     (.0021, .0019, .0015, .0011, .0009, .001, .0015, .0012, .0015, .0020, .0040, .0054)),
    ((.0030, .0029, .0028, .0027, .00275, .003, .0034, .004, .0049, .0057, .0066, .0074,
      .0082, .009, .0098, .0106, .0113),
     (.0025, .0023, .0024, .0028, .0032, .0034, .003, .0025, .0024, .003, .0045, .006,
      .0075, .009, .0105, .012, .0135),
     (.0005, .0006, .0007, .0008, .001, .0012, .0016, .0019, .0023, .0027, .0031, .0035,
      .0038, .0041, .0044, .0047, .0049),
     (.0013, .0012, .001, .0008, .0007, .0006, .0005, .00045, .0004, .00035, .0004, .0007,
      .0023, .0038, .0041, .0044, .0047),
     (.0005, .0004, .0003, .0002, .0001, 0., 0., 0., .0001, .0002, .00025, .0006, .0012,
      .003, .0035, .0038, .0047)),
)


def cdl(x, z, tc, fix_x8=False):
    """CDL: drag due to lift (McCormick fig. 2-52) at X = (CL - CLopt)/(CLmax - CLopt)."""
    xs = (_X2, _X4, _X6, _X8_FIXED if fix_x8 else _X8)
    x, z, tc = float(x), float(z), float(tc)
    dlcd = [linear(_TCP, [linear(xs[i], d, x) for d in _DCD[i]], tc) for i in range(4)]
    return linear(_ZP, dlcd, z)


_Z3 = (0.0, .02, .04, .06, .08, .10)
_TC2 = (.02, .04, .06, .08, .10, .12, .14, .16, .18, .20, .22, .24)
_DMDATA = (
    (.00465, .00478, .0049, .00515, .0055, .00565, .00578, .0064, .0069, .0074, .00765, .00865),
    (.005, .0051, .0054, .0058, .0062, .0067, .0069, .0076, .0083, .0086, .0089, .0092),
    (.0056, .00569, .00545, .00661, .0071, .00785, .0081, .0084, .0088, .009, .0092, .00953),
    (.0048, .0052, .006, .0065, .0073, .0081, .0083, .0084, .0091, .0094, .0097, .0103),
    (.0046, .0049, .0057, .0070, .00793, .00851, .0089, .0092, .00973, .0103, .0107, .0111),
    (.0085, .009, .00903, .00948, .00952, .010, .01003, .01053, .01103, .01107, .01203, .01253),
)


def cdmn(z, tc):
    """CDMN: minimum drag for 6e5 < Re < 5e6 (McCormick fig. 2-53)."""
    return linear(_Z3, [linear(_TC2, c, float(tc)) for c in _DMDATA], float(z))


_R = (1.E5, 1.5E5, 2.E5, 3.E5, 4.E5, 5.E5, 6.E5, 6.E6, 8.E6)
_ER = (3.8, 2.8, 2.2, 1.6, 1.3, 1.1, 1.0, 1.0, .95)


def recd(re):
    """RECD: Reynolds number factor on minimum drag (McCormick fig. 2-54)."""
    return linear(_R, _ER, float(re))


class Coefficients4(NamedTuple):
    cl: float
    cd: float
    re: float
    cld: float
    z: float          # maximum camber y/c (at x/c 0.5)
    clmax: float      # CLMX (MAXCL)
    clopt: float      # optimum CL (minimum drag)
    x: float          # (CL - CLopt)/(CLmax - CLopt)
    cdmin: float
    delcd: float


def naca4(alpha_deg, mach, tc, *, cld=None, z=None, chord_ft=1.0, altitude_ft=0.0,
          re=None, fix_x8=False) -> Coefficients4:
    """The 4-digit main program's CL/CD for one point (report pp. 103-104; NACA4 in the
    transcription). Give the design CL `cld` (camber from CLDVSZ) or the camber `z`
    (y/c). alpha_deg is measured from the design attitude: CL = CLD + dCL/dalpha alpha,
    Prandtl-Glauert. re None (or <= 0): from Mach, chord and altitude.

    Ranges the program warns outside of: Re 4e5-1e7, camber to 0.09. Its tables span
    t/c 0.04-0.18 (CL/CD), camber 0.02-0.08 and X -1.1 to 0.5; LINEAR extrapolates."""
    if (cld is None) == (z is None):
        raise ValueError("give exactly one of cld and z")
    alpha, mach, tc = float(alpha_deg), float(mach), float(tc)
    if re is None or not re > 0.0:
        atm = atmosphere(altitude_ft)
        m = 0.05 if mach == 0.0 else mach
        vel = m * atm.wa
        re = fmath.div(atm.den * vel * float(chord_ft), atm.vis)
    re = float(re)
    camp = 0 if z is None else 1
    cld, z = cldvsz(camp, tc, 0.0 if cld is None else cld, 0.0 if z is None else z)
    cre = -14.0 * (_log10(re) - 6.0) + 24.250
    dtc = -3.84375 * tc + 1
    clopt = z * cre * dtc
    clmx = maxcl(z, tc)
    if alpha == 0.0:
        delcl = 0.0
    else:
        k = (.1321 if tc >= .20 else .1216 if tc >= .17 else .1175 if tc >= .14 else
             .1134 if tc >= .11 else .1043 if tc >= .08 else .1039 if tc >= .05 else .1021)
        dcldal = fmath.div(k * (1 - tc), fmath.sqrt(1 - tc * tc))
        dcldal = fmath.div(dcldal, fmath.sqrt(1 - mach * mach))
        delcl = dcldal * alpha
    cl = cld + delcl
    x = fmath.div(cl - clopt, clmx - clopt)
    delcd = cdl(x, z, tc, fix_x8)
    cdmin = cdmn(z, tc) * recd(re)
    return Coefficients4(cl, cdmin + delcd, re, cld, z, clmx, clopt, x, cdmin, delcd)


def _log10(x):
    """ALOG10 with IEEE results for x <= 0 and nan."""
    if x > 0.0:
        return math.log10(x)
    return -math.inf if x == 0.0 else math.nan


# ----------------------------------------------------------------------------- XROTOR sections
def _re_chord16(re, mach, atm):
    """Chord (ft) giving Reynolds number re at Mach `mach` in the 16-series bank."""
    return re * atm.vis / (atm.den * mach * atm.wa)


def bank16(tc, cld, *, x=0.75, altitude_ft=0.0):
    """The 16-series bank as cl_cd(alpha_deg, mach, re) for `fit_section`."""
    atm = atmosphere(altitude_ft)
    a0 = naca16_alpha0(tc, cld)

    def cl_cd(alpha_deg, mach, re):
        c = clcd16(mach, alpha_deg - a0, 1, _re_chord16(re, mach, atm), tc, cld, x,
                   atm.wa, atm.vis, atm.den)
        return c.cl, c.cd
    return cl_cd


def bank4(tc, *, cld=None, z=None, fix_x8=False):
    """The 4-digit bank as cl_cd(alpha_deg, mach, re) for `fit_section` (angle from the
    design attitude, as in `naca4`)."""
    def cl_cd(alpha_deg, mach, re):
        c = naca4(alpha_deg, mach, tc, cld=cld, z=z, re=re, fix_x8=fix_x8)
        return c.cl, c.cd
    return cl_cd


def fit_section(cl_cd, xisect=0.0, *, mach=0.3, re=2.0e6, clmax, clmin, cmcon=-0.1,
                mcrit_default=None, alpha_window=(-2.0, 8.0), cl_window=0.35):
    """XROTOR AeroSection parameters fitted to a data bank `cl_cd(alpha_deg, mach, re)`:

      a0, dclda   the lift line at `mach`, with XROTOR's Prandtl-Glauert factor taken out
                  (dclda per rad, incompressible);
      cdmin, cldmin, dcdcl2   a least-squares parabola CD = CDMIN + DCDCL2 (CL - CLDMIN)^2
                  through the bank's drag within cl_window of its lowest point over the
                  angles alpha_window (deg from zero lift, CL about -0.2 to 0.85). The banks'
                  drag is
                  not a parabola (the Clark Y thick tables are coarse steps, and its drag
                  dips where the tables extrapolate past ALPHAO 9-11 deg), so the fit, not
                  the raw minimum, is what XROTOR's law can follow;
      reref, rexp the fit's Re, and the exponent from CD(CLDMIN) at 2 re;
      mcrit       the Mach number where CD at CLDMIN (at constant Re) has risen 0.002,
                  which is where XROTOR's drag rise CDC = 10 (M - Mcrit_eff)^3 equals 0.002
                  at CL = CLDMIN. mcrit_default if the bank has no rise below Mach 0.98 (the
                  4-digit bank has no compressibility drag).

    The banks give no stall: clmax and clmin are the caller's (the 4-digit bank's MAXCL is
    one estimate)."""
    from .aero import AeroSection
    deg = math.pi / 180.0
    c1, _ = cl_cd(-2.0, mach, re)
    c2, _ = cl_cd(4.0, mach, re)
    slope = (c2 - c1) / 6.0                              # per degree
    alpha_zero = -2.0 - c1 / slope
    beta = math.sqrt(1.0 - mach * mach)

    def alpha_at(cl, m=mach):
        a1, a2 = -2.0, 4.0
        l1, _ = cl_cd(a1, m, re)
        l2, _ = cl_cd(a2, m, re)
        return a1 + (cl - l1) * (a2 - a1) / (l2 - l1)

    grid = _frange(alpha_zero + alpha_window[0], alpha_zero + alpha_window[1], 0.05)
    cls = [cl_cd(a, mach, re)[0] for a in grid]
    cds = [cl_cd(a, mach, re)[1] for a in grid]
    # least-squares CD = p0 + p1 CL + p2 CL^2 near the lowest drag, sums in order
    low = min(cds)                          # centre: the middle of the lowest (flat) part
    flat = [c for c, d in zip(cls, cds) if d <= low + 1e-6]
    centre = 0.5 * (flat[0] + flat[-1])
    pts = [(c, d) for c, d in zip(cls, cds) if abs(c - centre) <= cl_window]
    cls, cds = [c for c, _ in pts], [d for _, d in pts]
    n = len(cls)
    sx = [sum(c ** k for c in cls) for k in range(5)]
    sy = [sum(d * c ** k for c, d in zip(cls, cds)) for k in range(3)]
    m3 = [[sx[i + j] for j in range(3)] for i in range(3)]
    p0, p1, p2 = _solve3(m3, sy)
    if p2 > 0.0:
        cldmin = -p1 / (2.0 * p2)
        cdmin = p0 + p1 * cldmin + p2 * cldmin * cldmin
        dcdcl2 = p2
    else:                                   # no curvature: flat at the window's minimum
        i = min(range(n), key=cds.__getitem__)
        cldmin, cdmin, dcdcl2 = cls[i], cds[i], 0.0

    def cd_at(m, r=re):
        return cl_cd(alpha_at(cldmin, m), m, r)[1]
    rexp = math.log(cd_at(mach, 2 * re) / cd_at(mach)) / math.log(2.0)
    base = cd_at(mach)
    mcrit = mcrit_default
    lo, hi = mach, None
    for m in _frange(mach, 0.98, 0.01):
        if cd_at(m) - base > 0.002:
            hi = m
            break
        lo = m
    if hi is not None:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if cd_at(mid) - base <= 0.002 else (lo, mid)
        mcrit = 0.5 * (lo + hi)
    return AeroSection(xisect=float(xisect), a0=alpha_zero * deg, clmax=clmax, clmin=clmin,
                       dclda=slope / deg * beta, dclda_stall=0.1, dcl_stall=0.1, cdmin=cdmin,
                       cldmin=cldmin, dcdcl2=dcdcl2, cmcon=cmcon, reref=re, rexp=rexp,
                       mcrit=mcrit if mcrit is not None else 0.8)


def _solve3(m, y):
    """3x3 linear solve by Cramer's rule."""
    def det(a):
        return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    d = det(m)
    out = []
    for k in range(3):
        mk = [[y[i] if j == k else m[i][j] for j in range(3)] for i in range(3)]
        out.append(det(mk) / d)
    return out


def _frange(a, b, step):
    n = int(round((b - a) / step))
    return [a + k * step for k in range(n + 1)]


# ----------------------------------------------------------------------------- Clark Y
# Report chapter I (Haines & Monaghan, R&M 2036; MacDougall, R&M 2224); the tables are the
# listing's DATA statements (reference/cr176883/clarky.f), generated from it.
_CY_AOF_A = (0.1095, 0.1094, 0.1091, 0.1085, 0.1069, 0.1043, 0.1017, 0.0991, 0.0999, 0.0999,
             0.0999, 0.0998, 0.0996, 0.0992, 0.0987, 0.098, 0.0972, 0.0962, 0.095, 0.0936, 0.092,
             0.0902, 0.0881, 0.0858, 0.0832, 0.0804, 0.0773, 0.0738)
_CY_AOF_H = (0.001, 0.002, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15,
             0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
             0.3)
_CY_CLSF_C = (0.01, 0.017, 0.022, 0.025, 0.026, 0.025, 0.021, 0.015, 0.007, -0.002, -0.014,
              -0.027, -0.04, -0.054, -0.07, -0.086, -0.102, -0.119, -0.138, -0.157, -0.175,
              -0.193, -0.212, -0.231, -0.25, -0.27, -0.289, -0.307, -0.325, -0.344, -0.363,
              -0.381, -0.4, -0.417, -0.432, -0.446, -0.459, -0.47, -0.478, -0.485)
_CY_CLSF_S = (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,
              0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28,
              0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4)
_CY_THK_A15 = (9.2, 8.95, 8.7, 8.4, 8.0, 7.9)
_CY_THK_A17 = (9.6, 9.4, 9.2, 8.95, 8.8)
_CY_THK_A19 = (10.1, 10.0, 9.9, 9.8)
_CY_THK_H = (0.4, 0.45, 0.5, 0.55, 0.6, 0.65)
_CY_THK_S = (0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.3, 0.35, 0.4)
_CY_HMLF_A = (3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5,
              11.0)
_CY_HMLF_B = (0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1,
              0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16)
_CY_HMLF_C = (  # C(I, J): [J-1][I-1], I thickness (B), J ALPHAO (A)
    (0.879, 0.858, 0.83125, 0.8125, 0.79125, 0.7775, 0.765, 0.751, 0.742, 0.732, 0.723, 0.713,
     0.704, 0.695, 0.686, 0.677, 0.669, 0.66, 0.651, 0.643, 0.635, 0.627, 0.618, 0.61, 0.603),
    (0.874, 0.852, 0.828, 0.807, 0.79, 0.775, 0.764, 0.751, 0.742, 0.732, 0.723, 0.713, 0.704,
     0.695, 0.686, 0.677, 0.669, 0.66, 0.651, 0.643, 0.635, 0.627, 0.618, 0.61, 0.603),
    (0.867, 0.846, 0.825, 0.804, 0.7875, 0.773, 0.765, 0.749, 0.74, 0.73, 0.721, 0.711, 0.702,
     0.694, 0.685, 0.676, 0.668, 0.658, 0.649, 0.641, 0.633, 0.625, 0.616, 0.608, 0.601),
    (0.859, 0.838, 0.817, 0.798, 0.78125, 0.768, 0.756, 0.745, 0.735, 0.725, 0.716, 0.707, 0.698,
     0.689, 0.68, 0.671, 0.663, 0.654, 0.645, 0.637, 0.628, 0.62, 0.611, 0.603, 0.596),
    (0.848, 0.828, 0.807, 0.789, 0.774, 0.761, 0.749, 0.738, 0.728, 0.718, 0.709, 0.699, 0.689,
     0.68, 0.671, 0.663, 0.654, 0.646, 0.637, 0.627, 0.618, 0.612, 0.602, 0.593, 0.585),
    (0.836, 0.817, 0.797, 0.78, 0.765, 0.752, 0.74, 0.729, 0.719, 0.709, 0.699, 0.689, 0.679,
     0.669, 0.66, 0.65, 0.641, 0.631, 0.621, 0.611, 0.6, 0.589, 0.578, 0.568, 0.559),
    (0.823, 0.805, 0.787, 0.771, 0.755, 0.742, 0.73, 0.719, 0.707, 0.697, 0.687, 0.675, 0.666,
     0.656, 0.646, 0.635, 0.625, 0.614, 0.603, 0.591, 0.579, 0.567, 0.555, 0.543, 0.532),
    (0.808, 0.792, 0.775, 0.76, 0.744, 0.731, 0.718, 0.706, 0.695, 0.684, 0.674, 0.663, 0.652,
     0.641, 0.629, 0.618, 0.607, 0.595, 0.584, 0.572, 0.56, 0.548, 0.536, 0.524, 0.513),
    (0.794, 0.778, 0.762, 0.745, 0.73, 0.717, 0.704, 0.691, 0.679, 0.667, 0.657, 0.646, 0.635,
     0.623, 0.611, 0.6, 0.588, 0.577, 0.565, 0.553, 0.541, 0.53, 0.518, 0.507, 0.496),
    (0.781, 0.765, 0.748, 0.731, 0.716, 0.702, 0.689, 0.676, 0.663, 0.651, 0.64, 0.628, 0.617,
     0.605, 0.593, 0.582, 0.57, 0.558, 0.546, 0.533, 0.523, 0.513, 0.502, 0.492, 0.483),
    (0.768, 0.752, 0.734, 0.717, 0.701, 0.687, 0.673, 0.661, 0.648, 0.635, 0.622, 0.611, 0.599,
     0.587, 0.575, 0.563, 0.551, 0.539, 0.528, 0.517, 0.507, 0.498, 0.489, 0.479, 0.47),
    (0.755, 0.739, 0.72, 0.703, 0.685, 0.67, 0.656, 0.644, 0.631, 0.618, 0.606, 0.593, 0.58,
     0.568, 0.556, 0.545, 0.534, 0.523, 0.513, 0.503, 0.493, 0.484, 0.474, 0.466, 0.457),
    (0.742, 0.726, 0.706, 0.689, 0.667, 0.652, 0.638, 0.626, 0.613, 0.601, 0.588, 0.575, 0.563,
     0.551, 0.539, 0.528, 0.518, 0.508, 0.498, 0.489, 0.479, 0.471, 0.462, 0.454, 0.446),
    (0.729, 0.713, 0.692, 0.675, 0.647, 0.633, 0.62, 0.608, 0.595, 0.583, 0.57, 0.558, 0.546,
     0.535, 0.524, 0.513, 0.503, 0.494, 0.484, 0.475, 0.467, 0.46, 0.453, 0.446, 0.439),
    (0.716, 0.7, 0.678, 0.661, 0.625, 0.613, 0.601, 0.589, 0.576, 0.564, 0.552, 0.541, 0.53,
     0.519, 0.509, 0.499, 0.49, 0.481, 0.473, 0.465, 0.458, 0.451, 0.445, 0.438, 0.432),
    (0.703, 0.687, 0.664, 0.647, 0.601, 0.592, 0.582, 0.57, 0.558, 0.547, 0.536, 0.526, 0.516,
     0.505, 0.496, 0.487, 0.479, 0.471, 0.463, 0.457, 0.45, 0.444, 0.439, 0.432, 0.426),
    (0.69, 0.674, 0.65, 0.633, 0.575, 0.57, 0.564, 0.552, 0.541, 0.531, 0.521, 0.511, 0.502,
     0.493, 0.485, 0.477, 0.47, 0.463, 0.456, 0.45, 0.445, 0.439, 0.433, 0.427, 0.422),
)
_CY_CDFF_A = (0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.19, 0.21, 0.23, 0.25, 0.3,
              0.35)
_CY_BOF_A = (-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8,
             1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2,
             4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6,
             7.8, 8.0)
_CY_BOF_B = (1.131, 1.114, 1.098, 1.082, 1.067, 1.054, 1.041, 1.03, 1.019, 1.009, 1.0, 0.994,
             0.987, 0.983, 0.979, 0.976, 0.974, 0.973, 0.972, 0.971, 0.972, 0.974, 0.976, 0.978,
             0.98, 0.982, 0.985, 0.989, 0.994, 0.999, 1.005, 1.012, 1.02, 1.028, 1.037, 1.047,
             1.057, 1.067, 1.078, 1.09, 1.103, 1.117, 1.132, 1.147, 1.162, 1.179, 1.197, 1.216,
             1.235, 1.255, 1.276)
_CY_COF_A = (0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09,
             0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155,
             0.16)
_CY_COF_B = (0.00745, 0.00753, 0.00761, 0.0077, 0.0078, 0.0079, 0.008, 0.0081, 0.0082, 0.0083,
             0.0084, 0.00851, 0.00862, 0.00874, 0.00886, 0.00898, 0.0091, 0.00921, 0.00932,
             0.00944, 0.00957, 0.0097, 0.00983, 0.00997, 0.0101, 0.01025, 0.0104)
_CY_CDSF_A = (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,
              0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28,
              0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4)
_CY_CDSF_B = (0.0001, 0.0009, 0.0022, 0.004, 0.0063, 0.0091, 0.0122, 0.0154, 0.0187, 0.0221,
              0.0255, 0.0288, 0.0322, 0.0356, 0.0389, 0.0423, 0.0456, 0.049, 0.0524, 0.0557,
              0.0591, 0.0625, 0.0658, 0.0692, 0.0726, 0.0759, 0.0793, 0.0826, 0.086, 0.0894,
              0.0927, 0.0961, 0.0995, 0.1028, 0.1062, 0.1096, 0.1129, 0.1163, 0.1196, 0.123)
_CY_HMDF_A = (-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
              6.5, 7.0)
_CY_HMDF_B = (0.045, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16)
_CY_HMDF_C = (  # C(I, J): [J-1][I-1], I thickness (B), J ALPHAO (A)
    (0.718, 0.705, 0.674, 0.627, 0.62, 0.5915, 0.56, 0.5255, 0.5015, 0.4715, 0.4435, 0.412, 0.378),
    (0.763, 0.749, 0.718, 0.675, 0.6575, 0.629, 0.599, 0.5705, 0.5425, 0.514, 0.487, 0.4575,
     0.4265),
    (0.8, 0.785, 0.755, 0.715, 0.69, 0.6625, 0.634, 0.6065, 0.5795, 0.5525, 0.5265, 0.499, 0.471),
    (0.82, 0.811, 0.782, 0.75, 0.728, 0.709, 0.691, 0.674, 0.658, 0.643, 0.629, 0.614, 0.6),
    (0.81, 0.817, 0.8, 0.774, 0.75, 0.728, 0.709, 0.69, 0.672, 0.656, 0.64, 0.625, 0.613),
    (0.787, 0.805, 0.808, 0.79, 0.768, 0.746, 0.725, 0.705, 0.685, 0.667, 0.651, 0.637, 0.625),
    (0.762, 0.785, 0.806, 0.799, 0.776, 0.755, 0.735, 0.715, 0.696, 0.678, 0.661, 0.647, 0.637),
    (0.735, 0.761, 0.796, 0.8, 0.774, 0.752, 0.732, 0.715, 0.698, 0.682, 0.667, 0.655, 0.646),
    (0.707, 0.735, 0.779, 0.795, 0.766, 0.741, 0.716, 0.699, 0.684, 0.671, 0.662, 0.654, 0.648),
    (0.68, 0.707, 0.756, 0.781, 0.752, 0.724, 0.697, 0.681, 0.667, 0.657, 0.65, 0.643, 0.637),
    (0.651, 0.68, 0.731, 0.761, 0.733, 0.705, 0.678, 0.662, 0.648, 0.639, 0.632, 0.626, 0.621),
    (0.621, 0.65, 0.702, 0.737, 0.713, 0.684, 0.66, 0.644, 0.63, 0.62, 0.611, 0.605, 0.6),
    (0.593, 0.62, 0.671, 0.708, 0.69, 0.672, 0.642, 0.626, 0.612, 0.602, 0.593, 0.586, 0.58),
    (0.565, 0.591, 0.641, 0.679, 0.663, 0.641, 0.625, 0.61, 0.597, 0.587, 0.58, 0.574, 0.57),
    (0.537, 0.563, 0.611, 0.649, 0.6605, 0.6715, 0.6825, 0.6935, 0.7055, 0.7185, 0.734, 0.7495,
     0.7665),
    (0.509, 0.535, 0.583, 0.62, 0.6325, 0.646, 0.6595, 0.671, 0.6835, 0.6965, 0.713, 0.729, 0.747),
    (0.481, 0.507, 0.556, 0.59, 0.601, 0.614, 0.627, 0.636, 0.646, 0.656, 0.6705, 0.684, 0.7),
    (0.453, 0.479, 0.529, 0.564, 0.5735, 0.588, 0.6025, 0.611, 0.6205, 0.6295, 0.644, 0.657,
     0.673),
)
_CY_AOIOF_A = (0.1032, 0.1012, 0.0991, 0.0952, 0.0914, 0.0876)
_CY_AOIOF_H = (0.06, 0.08, 0.1, 0.14, 0.18, 0.22)
# CDnnF: (AO, ((HM, CD) for each AO)), in the order the listing passes them to CDTHKF
_CY_CDTHK = {
    'CD9F': ((5.0, 6.0, 7.0, 8.0), (
        ((0.6545, 0.66364, 0.67273, 0.68182, 0.7, 0.8, 0.9),
         (0.00909, 0.00982, 0.01091, 0.01327, 0.01818, 0.05018, 0.08545)),
        ((0.62727, 0.63636, 0.64545, 0.65455, 0.7, 0.8, 0.9),
         (0.00945, 0.01, 0.01091, 0.01273, 0.02727, 0.06182, 0.1)),
        ((0.60909, 0.61818, 0.62727, 0.63636, 0.7, 0.8, 0.9),
         (0.01073, 0.01109, 0.01273, 0.01455, 0.03655, 0.07636, 0.12273)),
        ((0.6, 0.60909, 0.61818, 0.62727, 0.7, 0.8, 0.9),
         (0.01273, 0.01455, 0.01636, 0.01818, 0.04682, 0.09455, 0.14764)))),
    'CD10F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.7, 0.8, 0.9, 1.0),
         (0.0086, 0.009, 0.01, 0.0109, 0.013, 0.0158, 0.024, 0.058, 0.0906, 0.124)),
        ((0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.7, 0.8, 0.9, 1.0),
         (0.0094, 0.01, 0.0108, 0.012, 0.0138, 0.016, 0.03, 0.065, 0.105, 0.149)),
        ((0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.7, 0.8, 0.9, 1.0),
         (0.01, 0.0104, 0.011, 0.0121, 0.0141, 0.018, 0.04, 0.078, 0.128, 0.184)),
        ((0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.7, 0.8, 0.9, 0.97),
         (0.0112, 0.012, 0.0136, 0.016, 0.019, 0.0222, 0.05, 0.098, 0.154, 0.2)),
        ((0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.7, 0.8, 0.9, 0.91),
         (0.012, 0.0128, 0.0138, 0.016, 0.018, 0.0216, 0.067, 0.1224, 0.192, 0.2)),
        ((0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.7, 0.8, 0.86),
         (0.014, 0.0142, 0.016, 0.017, 0.0198, 0.022, 0.026, 0.03, 0.034, 0.086, 0.151, 0.2)),
        ((0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.6, 0.7, 0.8, 0.82),
         (0.016, 0.01618, 0.0164, 0.0176, 0.02, 0.023, 0.0262, 0.03, 0.045, 0.103, 0.178, 0.2)))),
    'CD11F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0), (
        ((0.62727, 0.63636, 0.64545, 0.65455, 0.66364, 0.67273, 0.7, 0.8, 0.9),
         (0.00903, 0.01, 0.01109, 0.01291, 0.01545, 0.01855, 0.02727, 0.06182, 0.09818)),
        ((0.6, 0.60909, 0.61818, 0.62727, 0.63636, 0.64545, 0.7, 0.8, 0.9),
         (0.00964, 0.01, 0.01091, 0.01209, 0.01455, 0.01636, 0.03327, 0.06909, 0.11018)),
        ((0.57273, 0.58182, 0.59091, 0.6, 0.60909, 0.61818, 0.7, 0.8, 0.9),
         (0.01036, 0.01064, 0.011, 0.01182, 0.01345, 0.01545, 0.04364, 0.08545, 0.134)),
        ((0.56364, 0.57273, 0.58182, 0.59091, 0.6, 0.60909, 0.7, 0.8, 0.9),
         (0.01091, 0.01136, 0.01273, 0.01455, 0.01636, 0.01909, 0.05455, 0.10455, 0.16273)),
        ((0.53636, 0.54545, 0.55454, 0.56364, 0.57273, 0.58182, 0.6, 0.64545, 0.7, 0.8, 0.88182),
         (0.01182, 0.01218, 0.01273, 0.01427, 0.01636, 0.01864, 0.02491, 0.03291, 0.07091,
          0.12909, 0.18)),
        ((0.52727, 0.53636, 0.54545, 0.55454, 0.56364, 0.6, 0.64545, 0.7, 0.8),
         (0.01273, 0.01309, 0.01455, 0.01636, 0.02, 0.036, 0.06, 0.08909, 0.15636)))),
    'CD12F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0), (
        ((0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.7, 0.8, 0.9),
         (0.0099, 0.0102, 0.0114, 0.013, 0.015, 0.017, 0.0198, 0.034, 0.066, 0.1)),
        ((0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.7, 0.8, 0.9),
         (0.0102, 0.011, 0.012, 0.013, 0.0147, 0.0164, 0.018, 0.02, 0.038, 0.074, 0.115)),
        ((0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.65, 0.7, 0.8, 0.9),
         (0.0112, 0.0114, 0.012, 0.0136, 0.0156, 0.0169, 0.019, 0.03, 0.048, 0.0886, 0.14)),
        ((0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.65, 0.7, 0.8, 0.9),
         (0.012, 0.0122, 0.014, 0.015, 0.0179, 0.02, 0.0221, 0.04, 0.06, 0.11, 0.168)),
        ((0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.6, 0.65, 0.7, 0.8, 0.9),
         (0.012, 0.013, 0.014, 0.016, 0.018, 0.02, 0.03, 0.05, 0.074, 0.1302, 0.2)),
        ((0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.6, 0.65, 0.7, 0.8, 0.85),
         (0.014, 0.0144, 0.016, 0.018, 0.02, 0.0237, 0.026, 0.04, 0.0638, 0.092, 0.158, 0.2)))),
    'CD13F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0), (
        ((0.58, 0.59, 0.6, 0.7, 0.9),
         (0.01, 0.0102, 0.011, 0.038, 0.1027)),
        ((0.57, 0.58, 0.59, 0.6, 0.63, 0.7, 0.8, 0.9),
         (0.0103, 0.011, 0.012, 0.014, 0.02, 0.042, 0.078, 0.12)),
        ((0.56, 0.57, 0.58, 0.59, 0.6, 0.7, 0.8, 0.9),
         (0.0112, 0.012, 0.0134, 0.0154, 0.018, 0.05, 0.096, 0.146)),
        ((0.54, 0.55, 0.56, 0.57, 0.58, 0.6, 0.7, 0.8, 0.9),
         (0.012, 0.013, 0.014, 0.016, 0.018, 0.0254, 0.064, 0.114, 0.174)),
        ((0.52, 0.53, 0.54, 0.55, 0.56, 0.6, 0.7, 0.8, 0.9),
         (0.0132, 0.0139, 0.0142, 0.016, 0.018, 0.034, 0.078, 0.1362, 0.2)),
        ((0.51, 0.52, 0.53, 0.54, 0.6, 0.7, 0.8, 0.845),
         (0.014, 0.01416, 0.016, 0.02, 0.042, 0.092, 0.1622, 0.2)))),
    'CD14F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0), (
        ((0.58, 0.59, 0.6, 0.61, 0.7, 0.8, 0.9),
         (0.01, 0.01018, 0.0112, 0.014, 0.04, 0.074, 0.1079)),
        ((0.57, 0.58, 0.59, 0.6, 0.61, 0.7, 0.8, 0.9),
         (0.0108, 0.012, 0.014, 0.0156, 0.018, 0.047, 0.084, 0.128)),
        ((0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.7, 0.8, 0.9),
         (0.0114, 0.012, 0.0132, 0.0156, 0.018, 0.02, 0.058, 0.102, 0.1526)),
        ((0.53, 0.54, 0.55, 0.56, 0.57, 0.6, 0.7, 0.8, 0.9),
         (0.012, 0.0128, 0.014, 0.016, 0.02, 0.0286, 0.07, 0.12, 0.1822)),
        ((0.51, 0.52, 0.53, 0.54, 0.55, 0.6, 0.7, 0.8, 0.88),
         (0.0136, 0.0138, 0.0144, 0.016, 0.02, 0.038, 0.082, 0.1424, 0.2)),
        ((0.5, 0.51, 0.52, 0.53, 0.6, 0.7, 0.8, 0.84),
         (0.014, 0.0144, 0.016, 0.018, 0.048, 0.098, 0.166, 0.2)))),
    'CD15F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0), (
        ((0.57, 0.58, 0.59, 0.6, 0.7, 0.8, 0.9),
         (0.01, 0.0104, 0.0114, 0.013, 0.042, 0.076, 0.11)),
        ((0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.7, 0.8, 0.9),
         (0.0104, 0.0108, 0.012, 0.014, 0.016, 0.018, 0.05, 0.09, 0.1334)),
        ((0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.6, 0.7, 0.8, 0.9),
         (0.0113, 0.012, 0.0124, 0.014, 0.016, 0.02, 0.025, 0.062, 0.106, 0.16)),
        ((0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.6, 0.7, 0.8, 0.9),
         (0.012, 0.0122, 0.014, 0.015, 0.0172, 0.02, 0.022, 0.032, 0.07378, 0.126, 0.19)),
        ((0.5, 0.51, 0.52, 0.53, 0.54, 0.6, 0.7, 0.8, 0.87),
         (0.0134, 0.0136, 0.014, 0.016, 0.02, 0.04, 0.086, 0.148, 0.2)),
        ((0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.6, 0.7, 0.8, 0.835),
         (0.0147, 0.0155, 0.017, 0.02, 0.024, 0.028, 0.0496, 0.1032, 0.172, 0.2)))),
    'CD16F': ((5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.55, 0.56667, 0.5833, 0.6, 0.65, 0.9),
         (0.01083, 0.01233, 0.015, 0.01833, 0.03333, 0.118)),
        ((0.51667, 0.53333, 0.55, 0.56667, 0.6, 0.65, 0.7, 0.8, 0.9),
         (0.01167, 0.01267, 0.015, 0.01833, 0.02667, 0.0425, 0.06, 0.10167, 0.14783)),
        ((0.4833, 0.5, 0.51667, 0.53333, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9),
         (0.01217, 0.01267, 0.015, 0.01742, 0.02167, 0.03333, 0.05167, 0.07033, 0.11833, 0.178)),
        ((0.46667, 0.48333, 0.5, 0.51667, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9),
         (0.01283, 0.013333, 0.01575, 0.01942, 0.02708, 0.04333, 0.06, 0.08267, 0.105, 0.13667,
          0.20192)),
        ((0.46667, 0.48333, 0.5, 0.5333, 0.55, 0.6, 0.7, 0.75, 0.8, 0.89417),
         (0.01333, 0.01583, 0.01933, 0.02667, 0.0315, 0.04833, 0.095, 0.12167, 0.15733, 0.23)),
        ((0.4625, 0.48333, 0.5, 0.525, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85833),
         (0.015, 0.01867, 0.02217, 0.03, 0.03708, 0.05667, 0.0815, 0.11, 0.14217, 0.18133, 0.23)),
        ((0.45, 0.45833, 0.46667, 0.5, 0.525, 0.6, 0.65, 0.7, 0.75, 0.8, 0.81667),
         (0.01667, 0.01717, 0.01917, 0.02667, 0.03667, 0.06617, 0.09367, 0.12633, 0.165, 0.212,
          0.23)))),
    'CD17F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.64, 0.65, 0.7, 0.8),
         (0.01, 0.011, 0.027, 0.068)),
        ((0.59, 0.6, 0.62, 0.63, 0.65, 0.7, 0.8),
         (0.01, 0.011, 0.013, 0.015, 0.021, 0.0275, 0.08)),
        ((0.57, 0.58, 0.6, 0.65, 0.7, 0.8),
         (0.01125, 0.012, 0.01545, 0.0325, 0.0525, 0.095)),
        ((0.54, 0.55, 0.56, 0.57, 0.58, 0.6, 0.61, 0.8),
         (0.01225, 0.0125, 0.015, 0.016, 0.0195, 0.025, 0.02975, 0.11)),
        ((0.52, 0.53, 0.54, 0.55, 0.57, 0.6, 0.75, 0.8),
         (0.014, 0.015, 0.0154, 0.018, 0.025, 0.0355, 0.1035, 0.13)),
        ((0.5, 0.51, 0.52, 0.53, 0.55, 0.6, 0.7, 0.8),
         (0.015, 0.016, 0.017, 0.02, 0.025, 0.044, 0.095, 0.1475)),
        ((0.49, 0.5, 0.51, 0.53, 0.55, 0.57, 0.6, 0.65, 0.7, 0.75, 0.8),
         (0.01625, 0.0175, 0.02, 0.025, 0.0315, 0.04, 0.052, 0.075, 0.1025, 0.13375, 0.1675)),
        ((0.47, 0.48, 0.5, 0.52, 0.55, 0.6, 0.65, 0.7, 0.8),
         (0.0175, 0.01875, 0.022, 0.0275, 0.03875, 0.06, 0.086, 0.1185, 0.19)),
        ((0.46, 0.47, 0.5, 0.52, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8),
         (0.01925, 0.02, 0.026, 0.034, 0.045, 0.07, 0.1, 0.135, 0.1965, 0.2175)))),
    'CD19F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.59, 0.61, 0.63, 0.65, 0.67, 0.8),
         (0.01, 0.0105, 0.01425, 0.02, 0.0275, 0.0925)),
        ((0.56, 0.59, 0.61, 0.63, 0.65, 0.8),
         (0.0105, 0.012, 0.01575, 0.02175, 0.03, 0.1025)),
        ((0.54, 0.56, 0.59, 0.61, 0.8),
         (0.0125, 0.014, 0.0175, 0.0245, 0.11625)),
        ((0.51, 0.54, 0.56, 0.59, 0.61, 0.75, 0.8),
         (0.014, 0.015, 0.0175, 0.025, 0.034, 0.1085, 0.1355)),
        ((0.49, 0.51, 0.54, 0.56, 0.59, 0.75, 0.8),
         (0.015, 0.0154, 0.02, 0.025, 0.035, 0.124, 0.1515)),
        ((0.47, 0.49, 0.51, 0.54, 0.56, 0.59, 0.75, 0.8),
         (0.016, 0.0175, 0.02, 0.025, 0.0315, 0.045, 0.1385, 0.171)),
        ((0.45, 0.47, 0.49, 0.51, 0.54, 0.56, 0.59, 0.75, 0.8),
         (0.0175, 0.01925, 0.02025, 0.025, 0.032, 0.04, 0.055, 0.1505, 0.1875)),
        ((0.45, 0.47, 0.49, 0.51, 0.54, 0.56, 0.59, 0.75, 0.8),
         (0.01925, 0.02125, 0.025, 0.03, 0.04, 0.049, 0.065, 0.17, 0.2085)),
        ((0.45, 0.47, 0.49, 0.51, 0.54, 0.56, 0.59, 0.75, 0.8),
         (0.02, 0.02475, 0.0295, 0.035, 0.0465, 0.059, 0.0754, 0.19125, 0.236)))),
    'CD21F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.57, 0.6, 0.62, 0.65, 0.67, 0.8),
         (0.01425, 0.015, 0.018, 0.026, 0.035, 0.11)),
        ((0.56, 0.58, 0.6, 0.64, 0.72, 0.8),
         (0.01425, 0.0155, 0.01925, 0.0305, 0.075, 0.125)),
        ((0.52, 0.55, 0.58, 0.6, 0.65, 0.8),
         (0.01425, 0.015, 0.02, 0.025, 0.045, 0.1405)),
        ((0.5, 0.55, 0.575, 0.6, 0.65, 0.7, 0.8),
         (0.015, 0.02, 0.025, 0.0335, 0.057, 0.0895, 0.1625)),
        ((0.47, 0.5, 0.55, 0.6, 0.65, 0.8),
         (0.016, 0.018, 0.0255, 0.0425, 0.0705, 0.18)),
        ((0.45, 0.48, 0.5, 0.55, 0.6, 0.65, 0.8),
         (0.0175, 0.01975, 0.021, 0.0345, 0.0515, 0.084, 0.195)),
        ((0.45, 0.48, 0.52, 0.57, 0.62, 0.7, 0.8),
         (0.0195, 0.022, 0.03045, 0.05, 0.075, 0.1305, 0.214)),
        ((0.44, 0.48, 0.5, 0.55, 0.65, 0.8),
         (0.02, 0.025, 0.03, 0.0505, 0.11, 0.231)),
        ((0.43, 0.45, 0.48, 0.55, 0.61, 0.71, 0.8),
         (0.022, 0.0235, 0.03, 0.06, 0.1, 0.175, 0.26)))),
    'CD23F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.55, 0.58, 0.6, 0.63, 0.65, 0.68, 0.7, 0.75, 0.8),
         (0.013, 0.015, 0.017, 0.025, 0.03075, 0.045, 0.056, 0.09, 0.129)),
        ((0.53, 0.56, 0.58, 0.6, 0.63, 0.67, 0.71, 0.8),
         (0.013, 0.015, 0.0175, 0.02075, 0.03, 0.05, 0.075, 0.145)),
        ((0.48, 0.52, 0.55, 0.58, 0.6, 0.64, 0.69, 0.75, 0.8),
         (0.013, 0.015, 0.0165, 0.2, 0.02545, 0.039, 0.075, 0.12, 0.161)),
        ((0.49, 0.52, 0.55, 0.59, 0.6, 0.63, 0.67, 0.735, 0.8),
         (0.015, 0.017, 0.0204, 0.03, 0.034, 0.0495, 0.075, 0.125, 0.181)),
        ((0.4, 0.45, 0.5, 0.54, 0.58, 0.6, 0.66, 0.71, 0.8),
         (0.0155, 0.016, 0.02, 0.025, 0.035, 0.043, 0.08, 0.12, 0.205)),
        ((0.4, 0.45, 0.5, 0.53, 0.55, 0.58, 0.6, 0.71, 0.8),
         (0.016, 0.0175, 0.02425, 0.03, 0.035, 0.045, 0.0535, 0.135, 0.222)),
        ((0.4, 0.45, 0.485, 0.55, 0.57, 0.6, 0.65, 0.71, 0.75),
         (0.0175, 0.02, 0.025, 0.043, 0.0505, 0.06625, 0.10075, 0.152, 0.192)),
        ((0.4, 0.42, 0.45, 0.5, 0.54, 0.585, 0.63, 0.69, 0.73),
         (0.0195, 0.02, 0.02325, 0.035, 0.05, 0.075, 0.105, 0.155, 0.19)),
        ((0.4, 0.43, 0.47, 0.5, 0.55, 0.58, 0.64, 0.72, 0.8),
         (0.0204, 0.025, 0.03, 0.0445, 0.07, 0.09, 0.135, 0.205, 0.281)))),
    'CD25F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.55, 0.58, 0.6, 0.624, 0.65, 0.68, 0.75, 0.8),
         (0.015, 0.0165, 0.02, 0.025, 0.0325, 0.05, 0.1, 0.1454)),
        ((0.52, 0.55, 0.58, 0.6, 0.64, 0.65, 0.68, 0.76, 0.8),
         (0.015, 0.0155, 0.02, 0.02425, 0.0355, 0.04, 0.06, 0.125, 0.166)),
        ((0.5, 0.56, 0.59, 0.62, 0.65, 0.71, 0.75, 0.8),
         (0.015, 0.02, 0.025, 0.035, 0.05, 0.095, 0.135, 0.19)),
        ((0.46, 0.51, 0.55, 0.6, 0.63, 0.67, 0.725, 0.8),
         (0.015, 0.0175, 0.0235, 0.035, 0.05, 0.0775, 0.125, 0.2125)),
        ((0.4, 0.45, 0.5, 0.55, 0.59, 0.63, 0.68, 0.76, 0.8),
         (0.0155, 0.017, 0.0205, 0.03, 0.04, 0.06, 0.1, 0.18, 0.233)),
        ((0.4, 0.45, 0.5, 0.55, 0.59, 0.61, 0.65, 0.74, 0.8),
         (0.0175, 0.02, 0.025, 0.035, 0.05, 0.06, 0.0825, 0.18, 0.25)),
        ((0.4, 0.43, 0.47, 0.5, 0.55, 0.58, 0.61, 0.67, 0.8),
         (0.0195, 0.02, 0.025, 0.03, 0.045, 0.0555, 0.075, 0.126, 0.2675)),
        ((0.4, 0.44, 0.49, 0.53, 0.58, 0.65, 0.73, 0.8),
         (0.02, 0.025, 0.035, 0.04975, 0.076, 0.132, 0.207, 0.285)),
        ((0.4, 0.46, 0.5, 0.55, 0.62, 0.71, 0.8),
         (0.021, 0.035, 0.05, 0.078, 0.13, 0.215, 0.305)))),
    'CD30F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.4, 0.55, 0.585, 0.63, 0.65, 0.69, 0.735, 0.76, 0.78),
         (0.015, 0.018, 0.02, 0.03, 0.037, 0.06, 0.1, 0.13, 0.16)),
        ((0.4, 0.53, 0.59, 0.64, 0.68, 0.72, 0.757, 0.8),
         (0.015, 0.018, 0.024, 0.04, 0.065, 0.102, 0.15, 0.22)),
        ((0.4, 0.51, 0.57, 0.62, 0.65, 0.7, 0.74, 0.8),
         (0.015, 0.018, 0.025, 0.04, 0.055, 0.1, 0.15, 0.25)),
        ((0.4, 0.5, 0.55, 0.6, 0.64, 0.66, 0.7, 0.78, 0.8),
         (0.016, 0.02, 0.025, 0.0395, 0.06, 0.077, 0.12, 0.24, 0.275)),
        ((0.46, 0.53, 0.58, 0.62, 0.65, 0.7, 0.74, 0.8),
         (0.02, 0.03, 0.04, 0.059, 0.081, 0.14, 0.199, 0.3)),
        ((0.4, 0.45, 0.5, 0.55, 0.6, 0.68, 0.665, 0.7, 0.75, 0.8),
         (0.02, 0.021, 0.03, 0.041, 0.061, 0.082, 0.12, 0.165, 0.24, 0.33)),
        ((0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.72, 0.8),
         (0.02, 0.025, 0.033, 0.05, 0.08, 0.133, 0.22, 0.35)),
        ((0.4, 0.46, 0.5, 0.54, 0.59, 0.65, 0.73, 0.8),
         (0.02, 0.03, 0.04, 0.06, 0.1, 0.16, 0.2608, 0.365)),
        ((0.4, 0.45, 0.5, 0.56, 0.61, 0.72, 0.8),
         (0.021, 0.035, 0.06, 0.105, 0.15, 0.277, 0.381)))),
    'CD35F': ((3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0), (
        ((0.57, 0.6, 0.62, 0.65, 0.68, 0.71, 0.75, 0.78, 0.8),
         (0.02, 0.025, 0.0308, 0.042, 0.06, 0.09, 0.15, 0.205, 0.25)),
        ((0.55, 0.6, 0.65, 0.68, 0.705, 0.76, 0.8),
         (0.02, 0.029, 0.05, 0.075, 0.1, 0.19, 0.28)),
        ((0.53, 0.58, 0.63, 0.67, 0.7, 0.73, 0.78, 0.8),
         (0.02, 0.028, 0.05, 0.08, 0.111, 0.16, 0.26, 0.312)),
        ((0.52, 0.565, 0.6, 0.63, 0.67, 0.71, 0.75, 0.8),
         (0.02, 0.03, 0.043, 0.06, 0.1, 0.15, 0.22, 0.342)),
        ((0.48, 0.54, 0.59, 0.62, 0.65, 0.7, 0.74, 0.8),
         (0.02, 0.03, 0.05, 0.07, 0.095, 0.156, 0.23, 0.375)),
        ((0.46, 0.54, 0.59, 0.62, 0.65, 0.7, 0.75, 0.8),
         (0.02, 0.036, 0.06, 0.085, 0.112, 0.188, 0.29, 0.405)),
        ((0.43, 0.46, 0.52, 0.564, 0.6, 0.65, 0.75, 0.8),
         (0.02, 0.025, 0.04, 0.06, 0.088, 0.144, 0.325, 0.43)),
        ((0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8),
         (0.02, 0.03, 0.05, 0.079, 0.12, 0.27, 0.45)),
        ((0.4, 0.44, 0.48, 0.538, 0.587, 0.68, 0.77, 0.8),
         (0.021, 0.04, 0.059, 0.1, 0.15, 0.28, 0.42, 0.47)))),
}


def _ylagn(x, y, n, z):
    """FUNCTION YLAGN (Clark Y bank): Lagrange interpolation through the first n points."""
    total = 0.0                                          # SUM
    for i in range(n):
        s = 1.0
        for j in range(n):
            if j != i:
                s = fmath.div(s * (z - x[j]), x[i] - x[j])
        total = total + s * y[i]
    return total


def _out(y, n, x):
    """SUBROUTINE OUT: (IS, IL), 0-based, of the n points bracketing x, stepping by
    IS = IL; the end group outside the table."""
    nmax = len(y)
    il = n - 1
    if x <= y[il]:
        return 0, il
    il = nmax - 1
    is_ = il - n + 1
    if x >= y[is_]:
        return is_, il
    is_ = 0
    while True:
        il = is_ + n - 1
        if (x - y[is_]) * (x - y[il]) <= 0.:
            return is_, il
        is_ = il


def _bracket(a, n, x):
    """The inline search most Clark Y functions use (first group, last group, else
    stepping IS = IL), as OUT but written out; same result."""
    return _out(a, n, x)


def _table(xs, ys, n, z, key):
    is_, il = _bracket(xs, n, key)
    return _ylagn(xs[is_:il + 1], ys[is_:il + 1], n, z)


def _aof(t):
    """FUNCTION AOF: low-speed lift slope (per deg) against t/c (H & M table 1). Its own
    search: the last pair once IL reaches 28. 0.1096 at t = 0."""
    if t == 0.:
        return 0.1096
    h = _CY_AOF_H
    is_, il = 0, 1
    while True:
        if il >= 27:
            is_, il = 26, 27
            break
        if (t - h[is_]) * (t - h[il]) <= 0.:
            break
        is_ = il
        il = is_ + 1
    return _ylagn(h[is_:il + 1], _CY_AOF_A[is_:il + 1], 2, t)


def _aoiof(t):
    """FUNCTION AOIOF: low-speed lift slope from NACA TR 628 for CLTMAX; 0.1096 below t 0.06."""
    if t < 0.06:
        return 0.1096
    h = _CY_AOIOF_H
    is_, il = 0, 1
    while True:
        if il >= 5:
            is_, il = 4, 5
            break
        if (t - h[is_]) * (t - h[il]) <= 0.:
            break
        is_ = il
        il = is_ + 1
    return _ylagn(h[is_:il + 1], _CY_AOIOF_A[is_:il + 1], 2, t)


def _clsf(h):
    """FUNCTION CLSF: lift change above the lift-critical Mach number, H = M - ML (H & M
    table 4), quadratic through three points."""
    return _table(_CY_CLSF_S, _CY_CLSF_C, 3, h, h)


def _thkslp(t, hm):
    """FUNCTION THKSLP: lift slope (per deg) of sections over 16% thick (H & M table 5,
    inverse slopes against t/c and Mach); AOF below Mach 0.3."""
    if hm < 0.3:
        return _aof(t)
    h = _CY_THK_H
    is_, il = _out(_CY_THK_S, 2, t)
    js, jl = _out(h, 2, hm)
    u = (h[js], h[jl])
    x = (_CY_THK_S[is_], _CY_THK_S[il])
    y = []
    for i in range(is_, il + 1):                     # GO TO (15,17,...,40), I
        tab = (_CY_THK_A15, _CY_THK_A17, _CY_THK_A19)[i] if i < 3 else None
        if tab is not None:
            ans = tab[0]
            if not hm <= h[0]:
                ans = tab[-1]
                if not hm >= h[len(tab) - 1]:
                    ans = _ylagn(u, (tab[js], tab[jl]), 2, hm)
        else:
            ans = (10.5, 10.9, 11.4, 13.6, 19.2, 30.0)[i - 3]
        y.append(ans)
    return fmath.div(1., _ylagn(x, y, 2, t))


def _hmlf(t, alphao):
    """FUNCTION HMLF: lift-critical Mach number ML against t/c and ALPHAO (H & M table 3).
    CR-176883: the switch to MacDougall's drag-critical table below 3 deg is commented out
    (IMPROV = 1), so ALPHAO < 3 uses the 3 deg column and labels 12-20 never run; they
    are not ported. The listing sets ALPHAO = 3 and restores it on return."""
    if alphao <= 3.:
        alphao = 3.
    b = _CY_HMLF_B
    is_, il = _bracket(b, 2, t)
    if not alphao > 3.:
        x = b[is_:il + 1]
        y = tuple(_CY_HMLF_C[0][is_:il + 1])
        return _ylagn(x, y, 2, t)
    a = _CY_HMLF_A
    js, jl = _bracket(a, 2, alphao)
    u = a[js:jl + 1]
    x, y = [], []
    for i in range(is_, il + 1):
        x.append(b[i])
        v = tuple(_CY_HMLF_C[j][i] for j in range(js, jl + 1))
        y.append(_ylagn(u, v, 2, alphao))
    return _ylagn(x, y, 2, t)


def _clff(alphao, t, hm):
    """FUNCTION CLFF: CL. Over 16% thick, ALPHAO times THKSLP; otherwise Prandtl-Glauert
    on AOF's slope up to ML (HMLF), then ML's value plus CLSF(M - ML)."""
    if not t < 0.16:
        return alphao * _thkslp(t, hm)
    ao = _aof(t)
    hml = _hmlf(t, alphao)
    h = hm - hml
    clo = ao * alphao
    clff = math.nan          # unset in the listing if neither branch applies (NaN ML)
    if hm <= hml:
        clff = fmath.div(clo, fmath.sqrt(1. - hm * hm))
    if hm > hml:
        clff = fmath.div(clo, fmath.sqrt(1. - hml * hml)) + _clsf(h)
    return clff


def _bof(alphao):
    """FUNCTION BOF: low-speed drag factor against ALPHAO (H & M table 6), quadratic."""
    return _table(_CY_BOF_A, _CY_BOF_B, 3, alphao, alphao)


def _cof(t):
    """FUNCTION COF: low-speed drag against t/c (H & M table 7)."""
    return _table(_CY_COF_A, _CY_COF_B, 2, t, t)


def _cdsf(hm, h):
    """FUNCTION CDSF: drag rise above the drag-critical Mach number, H = M - MD (MacDougall
    table 2). CR-176883: the interval is searched with HM, the Mach number, and the
    interpolation done at H: kept as listed."""
    return _table(_CY_CDSF_A, _CY_CDSF_B, 2, h, hm)


def _hmdf(t, alphao):
    """FUNCTION HMDF: drag-critical Mach number MD against t/c and ALPHAO (MacDougall
    tables 1A and 1B)."""
    b, a = _CY_HMDF_B, _CY_HMDF_A
    is_, il = _bracket(b, 2, t)
    js, jl = _bracket(a, 2, alphao)
    u = a[js:jl + 1]
    x, y = [], []
    for i in range(is_, il + 1):
        x.append(b[i])
        v = tuple(_CY_HMDF_C[j][i] for j in range(js, jl + 1))
        y.append(_ylagn(u, v, 2, alphao))
    return _ylagn(x, y, 2, t)


def _getf(xhm, ycd, hm):
    """FUNCTION GETF: CD against Mach in one CDnnF curve; its first value below the first Mach."""
    if not hm > xhm[0]:
        return ycd[0]
    js, jl = _out(xhm, 2, hm)
    return _ylagn(xhm[js:jl + 1], ycd[js:jl + 1], 2, hm)


def _cdthkf(name, alphao, hm):
    """FUNCTION CDTHKF with CDnnF's tables: linear in ALPHAO (not below the first AO)
    between the bracketing Mach curves."""
    ao, curves = _CY_CDTHK[name]
    alp = alphao
    if alp <= ao[0]:
        alp = ao[0]
    is_, il = _out(ao, 2, alp)
    y = [_getf(*curves[i], hm) for i in range(is_, il + 1)]
    return _ylagn(ao[is_:il + 1], y, 2, alp)


_CY_CDFF_FUN = ('CD9F', 'CD10F', 'CD11F', 'CD12F', 'CD13F', 'CD14F', 'CD15F', 'CD16F',
                'CD17F', 'CD19F', 'CD21F', 'CD23F', 'CD25F', 'CD30F', 'CD35F')


def _cdff(t, alphao, hm):
    """FUNCTION CDFF: CD. Thin sections at low incidence (t/c <= 0.16, and ALPHAO <= 5 or
    t/c <= 0.09): BOF*COF plus CDSF above MD. Otherwise linear in t/c between the two
    bracketing thickness families' CDnnF (9-35%)."""
    if not t > 0.16 and not (t > 0.09 and alphao > 5.):
        cdo = _bof(alphao) * _cof(t)
        hmd = _hmdf(t, alphao)
        cds = 0.
        h = hm - hmd
        if h > 0.:
            cds = _cdsf(hm, h)
        return cdo + cds
    a = _CY_CDFF_A
    is_, il = _bracket(a, 2, t)
    y = [_cdthkf(_CY_CDFF_FUN[i], alphao, hm) for i in range(is_, il + 1)]
    return _ylagn(a[is_:il + 1], y, 2, t)


class CoefficientsCY(NamedTuple):
    cl: float
    cd: float


def clcdcy(alphao, t, hm, ncl=1, ncd=1) -> CoefficientsCY:
    """The Clark Y bank's SUBROUTINE CLCD (report pp. 40-41; CLCDCY in the transcription):
    CL and CD at ALPHAO (deg, from zero lift), thickness ratio t, Mach hm. ncl/ncd 0 skip
    CL/CD (returned as NaN: the listing leaves them unset).

    Ranges the functions warn outside of: t/c 0.04-0.16 for lift and thin-section drag (to
    0.40 for THKSLP, 0.35 for the thick drag families), ALPHAO -2 to 11 deg (3-11 for the
    thick drag), Mach to 0.8 for THKSLP, 1.0 for CDFF; M - ML and M - MD 0.01-0.40."""
    alphao, t, hm = float(alphao), float(t), float(hm)
    cl = _clff(alphao, t, hm) if ncl != 0 else math.nan
    cd = _cdff(t, alphao, hm) if ncd != 0 else math.nan
    return CoefficientsCY(cl, cd)


class ClarkYDesign(NamedTuple):
    t: float          # thickness ratio (overwritten from CLD when CLD is given)
    cld: float        # design CL (NCLD)
    alphlo: float     # zero-lift angle, deg
    dncld: float      # 0: CLD computed from t; 1: t from CLD


def cltmax(t, ncld=0.0) -> ClarkYDesign:
    """SUBROUTINE CLTMAX (NACA TR 628): design CL from t/c when ncld is 0, else t/c from
    the design CL; and the zero-lift angle ALPHLO = -CLD / AOIOF(t). Warns outside t/c
    0.06-0.22 and CLD 0.30-0.90."""
    t, ncld = float(t), float(ncld)
    if not ncld != 0.:
        ncld = fmath.div(-0.16879 + fmath.sqrt(0.02849 + 0.496405 * t), 0.2482025)
        aoi = _aoiof(t)
        return ClarkYDesign(t, ncld, fmath.div(ncld, -aoi), 0.)
    t = 0.1241 * (ncld * ncld) + 0.16879 * ncld
    aoi = _aoiof(t)
    return ClarkYDesign(t, ncld, fmath.div(ncld, -aoi), 1.)


def clark_y(alpha_deg, mach, tc, *, cld=None) -> CoefficientsCY:
    """Clark Y section at geometric angle alpha_deg (the main program's ALPHA): CLTMAX's
    ALPHLO from tc (or tc from cld, which then replaces tc), then CLCD at ALPHA - ALPHLO."""
    d = cltmax(tc, 0.0 if cld is None else cld)
    return clcdcy(float(alpha_deg) - d.alphlo, d.t, mach)


def bank_clark_y(tc, *, cld=None):
    """The Clark Y bank as cl_cd(alpha_deg, mach, re) for `fit_section`. re is ignored: the
    bank's tables are R&M 2036's (and 2224's) test data at those tests' Reynolds numbers,
    and the program only prints a Reynolds number."""
    d = cltmax(tc, 0.0 if cld is None else cld)

    def cl_cd(alpha_deg, mach, re):
        c = clcdcy(float(alpha_deg) - d.alphlo, d.t, mach)
        return c.cl, c.cd
    return cl_cd


# ----------------------------------------------------------------------------- bank sections
# Drop-in replacements for pyxrotor.sections.clark_y / naca16 (same arguments): XROTOR
# sections fitted to the banks (fit_section at Mach 0.3, Re 2e6). The banks have no stall
# or moment, so CLMAX, CLMIN and CMCON are pyxrotor.sections' estimates. mcrit_shift is
# added to the bank's MCRIT (None counts as 0: the calibrated rule belongs to the estimates).
CLARK_Y_T_MAX = 0.35     # the Clark Y tables' thickest section
NACA16_T_MAX = 0.40      # the 16-series drag-divergence table's
NACA16_CLD_MIN = 0.2     # the 16-series camber-drag fit's lowest design CL (CLIT(1))


@functools.lru_cache(maxsize=None)
def _fit_clark_y(t):
    return fit_section(bank_clark_y(t), clmax=1.3, clmin=-0.5)


@functools.lru_cache(maxsize=None)
def _fit_naca16(t, cld, x):
    return fit_section(bank16(t, cld, x=x), clmax=1.0, clmin=-1.0)


def clark_y_section(t, xisect, mcrit_shift=0.0):
    """Clark Y section from the CR-176883 bank. Above CLARK_Y_T_MAX (shanks) the
    sections.clark_y estimate with its plain Mcrit + mcrit_shift. Re is not in the bank:
    REREF and REXP are the estimate's."""
    from . import sections
    t, shift = float(t), float(mcrit_shift or 0.0)
    est = sections.clark_y(t, xisect, mcrit_shift=shift)
    if t > CLARK_Y_T_MAX:
        return est
    f = _fit_clark_y(t)
    return dataclasses.replace(est, a0=f.a0, dclda=f.dclda, cdmin=f.cdmin, cldmin=f.cldmin,
                               dcdcl2=f.dcdcl2, mcrit=f.mcrit + shift)


def naca16_section(t, xisect, cld, mcrit_shift=0.0):
    """NACA 16-series section from the CR-176883 bank, its laminar fraction at r/R = xisect.
    Below NACA16_CLD_MIN the bank's camber drag extrapolates (it rises with Mach at CLD 0 as
    if supercritical by Mach 0.45), so the drag terms come from CLD 0.2 and only the zero-lift
    angle and CLDMIN follow cld. Above NACA16_T_MAX the sections.naca16 estimate."""
    from . import sections
    t, cld, shift = float(t), float(cld), float(mcrit_shift or 0.0)
    est = sections.naca16(t, xisect, cld, mcrit_shift=shift)
    if t > NACA16_T_MAX:
        return est
    c = max(cld, NACA16_CLD_MIN)
    f = _fit_naca16(t, c, min(max(float(xisect), 0.0), 1.0))
    a0 = f.a0 if c == cld else naca16_alpha0(t, cld) * math.pi / 180.0
    return dataclasses.replace(est, a0=a0, dclda=f.dclda, cdmin=f.cdmin,
                               cldmin=f.cldmin - c + cld, dcdcl2=f.dcdcl2, reref=f.reref,
                               rexp=f.rexp, mcrit=f.mcrit + shift)
