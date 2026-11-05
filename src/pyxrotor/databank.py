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
"""
from __future__ import annotations

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
                mcrit_default=None, cl_window=0.4):
    """XROTOR AeroSection parameters fitted to a data bank `cl_cd(alpha_deg, mach, re)`:

      a0, dclda   the lift line at `mach`, with XROTOR's Prandtl-Glauert factor taken out
                  (dclda per rad, incompressible);
      cdmin, cldmin, dcdcl2   CD's minimum over CL and a least-squares fit of
                  CD - CDMIN = DCDCL2 (CL - CLDMIN)^2 within cl_window of it;
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

    cls = [cl_cd(a, mach, re)[0] for a in _frange(alpha_zero - 4, alpha_zero + 14, 0.05)]
    cds = [cl_cd(a, mach, re)[1] for a in _frange(alpha_zero - 4, alpha_zero + 14, 0.05)]
    i = min(range(len(cds)), key=cds.__getitem__)
    cldmin, cdmin = cls[i], cds[i]
    pts = [(cl - cldmin, cd - cdmin) for cl, cd in zip(cls, cds) if abs(cl - cldmin) <= cl_window]
    dcdcl2 = sum(d * d * e for d, e in pts) / sum(d ** 4 for d, _ in pts)

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


def _frange(a, b, step):
    n = int(round((b - a) / step))
    return [a + k * step for k in range(n + 1)]
