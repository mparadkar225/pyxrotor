# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# cython: initializedcheck=False, nonecheck=False
"""Compiled versions of pyxrotor's hot loops (optional; see pyxrotor._accel).

Each function performs exactly the operations of its pure-Python counterpart, in the same
order, so results are bit-identical (and so equal to XROTOR 7.55's). This needs the C
compiler to keep IEEE semantics: setup.py builds this file with -ffp-contract=off (no fused
multiply-add) and without -ffast-math. cdivision=True gives IEEE division (x/0 -> inf, 0/0
-> nan), which is what pyxrotor.fmath.div returns. NaNs produced by sqrt/log of a negative
number are returned as the positive quiet NaN, like pyxrotor.fmath.

The pure-Python versions stay the reference implementation; tests check the two against
each other (tests/differential/test_speedups.py) as well as against the Fortran.
"""
from libc.math cimport sqrt as c_sqrt, exp as c_exp, log as c_log, pow as c_pow, fabs, NAN, INFINITY
from cython.view cimport array as cvarray


# ------------------------------------------------------------------ fmath equivalents
cdef inline double f_sqrt(double x) noexcept nogil:
    # fmath.sqrt: math.sqrt(x) if x >= 0 else nan
    if x >= 0.0:
        return c_sqrt(x)
    return NAN


cdef inline double f_log(double x) noexcept nogil:
    if x > 0.0:
        return c_log(x)
    if x == 0.0:
        return -INFINITY
    return NAN


cdef inline double f_fmax(double a, double b) noexcept nogil:
    # fmath.fmax: b if b > a else a (NaN -> a)
    return b if b > a else a


cdef inline double f_fmin(double a, double b) noexcept nogil:
    return b if b < a else a


# ------------------------------------------------------------------ CLCDCM
# aero.py's compressibility drag model constants
cdef double CDMFACTOR = 10.0
cdef double CLMFACTOR = 0.25
cdef double MEXP = 3.0
cdef double CDMDD = 0.0020
cdef double CDMSTALL = 0.1000


def clcdcm(double alf, double w, double rey, double vel, double vso, sec):
    """aero.clcdcm (SUBROUTINE CLCDCM); returns the SectionAero fields as a tuple."""
    cdef double a0 = sec.a0, clmax = sec.clmax, clmin = sec.clmin
    cdef double dclda = sec.dclda, dclda_stall = sec.dclda_stall, dcl_stall = sec.dcl_stall
    cdef double cdmin = sec.cdmin, cldmin = sec.cldmin, dcdcl2 = sec.dcdcl2, cmcon = sec.cmcon
    cdef double mcrit = sec.mcrit, reref = sec.reref, rexp = sec.rexp
    cdef double msq, msq_w, pg, pg_w, mach, mach_w, cla, cla_alf, cla_w, dmstall
    cdef double clmaxm, clminm, ecmax, ecmin, cllim, cllim_cla, fstall, clift, cl_alf, cl_w
    cdef double cmom, cm_al, cm_w, rcorr, rcorr_rey, dcl, cdrag, cd_alf, cd_w, cd_rey
    cdef double dcdx, dcd, dcd_alf, dcd_w, dmdd, critmach, critmach_alf, critmach_w
    cdef double cdc, cdc_alf, cdc_w, fac, fac_w
    cdef bint mach_limited = False, stallf

    msq = w * w * (vel * vel) / (vso * vso)
    msq_w = 2.0 * w * (vel * vel) / (vso * vso)
    if msq >= 1.0:
        mach_limited = True
        msq = 0.99
        msq_w = 0.0
    pg = 1.0 / f_sqrt(1.0 - msq)
    pg_w = 0.5 * msq_w * (pg * pg * pg)

    mach = f_sqrt(msq)
    mach_w = 0.0
    if mach != 0.0:
        mach_w = 0.5 * msq_w / mach

    cla = dclda * pg * (alf - a0)
    cla_alf = dclda * pg
    cla_w = dclda * pg_w * (alf - a0)

    dmstall = c_pow(CDMSTALL / CDMFACTOR, 1.0 / MEXP)
    clmaxm = f_fmax(0.0, (mcrit + dmstall - mach) / CLMFACTOR) + cldmin
    clmax = f_fmin(clmax, clmaxm)
    clminm = f_fmin(0.0, -(mcrit + dmstall - mach) / CLMFACTOR) + cldmin
    clmin = f_fmax(clmin, clminm)

    ecmax = c_exp(f_fmin(200.0, (cla - clmax) / dcl_stall))
    ecmin = c_exp(f_fmin(200.0, (clmin - cla) / dcl_stall))
    cllim = dcl_stall * f_log((1.0 + ecmax) / (1.0 + ecmin))
    cllim_cla = ecmax / (1.0 + ecmax) + ecmin / (1.0 + ecmin)

    fstall = dclda_stall / dclda
    clift = cla - (1.0 - fstall) * cllim
    cl_alf = cla_alf - (1.0 - fstall) * cllim_cla * cla_alf
    cl_w = cla_w - (1.0 - fstall) * cllim_cla * cla_w

    stallf = clift > clmax or clift < clmin

    cmom = pg * cmcon
    cm_al = 0.0
    cm_w = pg_w * cmcon

    if rey <= 0:
        rcorr = 1.0
        rcorr_rey = 0.0
    else:
        rcorr = c_pow(rey / reref, rexp)
        rcorr_rey = rexp / rey

    dcl = clift - cldmin
    cdrag = (cdmin + dcdcl2 * (dcl * dcl)) * rcorr
    cd_alf = (2.0 * dcdcl2 * dcl * cl_alf) * rcorr
    cd_w = (2.0 * dcdcl2 * dcl * cl_w) * rcorr
    cd_rey = cdrag * rcorr_rey

    fstall = dclda_stall / dclda
    dcdx = (1.0 - fstall) * cllim / (pg * dclda)
    dcd = 2.0 * (dcdx * dcdx)
    dcd_alf = 4.0 * dcdx * (1.0 - fstall) * cllim_cla * cla_alf / (pg * dclda)
    dcd_w = 4.0 * dcdx * ((1.0 - fstall) * cllim_cla * cla_w / (pg * dclda) - dcd / pg * pg_w)

    dmdd = c_pow(CDMDD / CDMFACTOR, 1.0 / MEXP)
    critmach = mcrit - CLMFACTOR * fabs(clift - cldmin) - dmdd
    critmach_alf = -CLMFACTOR * fabs(cl_alf)
    critmach_w = -CLMFACTOR * fabs(cl_w)
    if mach < critmach:
        cdc = 0.0
        cdc_alf = 0.0
        cdc_w = 0.0
    else:
        cdc = CDMFACTOR * c_pow(mach - critmach, MEXP)
        cdc_w = MEXP * mach_w * cdc / mach - MEXP * critmach_w * cdc / critmach
        cdc_alf = -(MEXP * critmach_alf * cdc / critmach)

    fac = 1.0
    fac_w = 0.0
    cdrag = fac * cdrag + dcd + cdc
    cd_alf = fac * cd_alf + dcd_alf + cdc_alf
    cd_w = fac * cd_w + fac_w * cdrag + dcd_w + cdc_w
    cd_rey = fac * cd_rey

    return (clift, cl_alf, cl_w, stallf, cdrag, cd_alf, cd_w, cd_rey,
            cmom, cm_al, cm_w, clmax, clmin, mach_limited)


# ------------------------------------------------------------------ GAUSS
def gauss(double[:, ::1] z, double[:, ::1] r, int nn):
    """numerics.gauss (SUBROUTINE GAUSS) on C-contiguous float64 arrays, in place.
    r has shape (NSIZ, NRHS)."""
    cdef int nrhs = r.shape[1]
    cdef int np_, np1, n, nx, k, l
    cdef double d, pivot, ztmp, temp
    with nogil:
        for np_ in range(nn - 1):
            np1 = np_ + 1
            # ---- find max pivot index NX (only '>' moves it; NaN moves it)
            nx = np_
            for n in range(np1, nn):
                d = fabs(z[n, np_]) - fabs(z[nx, np_])
                if not (d < 0.0 or d == 0.0):
                    nx = n
            pivot = 1.0 / z[nx, np_]
            # ---- switch pivots
            z[nx, np_] = z[np_, np_]
            # ---- switch rows & normalize pivot row
            for l in range(np1, nn):
                temp = z[nx, l] * pivot
                z[nx, l] = z[np_, l]
                z[np_, l] = temp
            for l in range(nrhs):
                temp = r[nx, l] * pivot
                r[nx, l] = r[np_, l]
                r[np_, l] = temp
            # ---- forward eliminate everything
            for k in range(np1, nn):
                ztmp = z[k, np_]
                for l in range(np1, nn):
                    z[k, l] = z[k, l] - ztmp * z[np_, l]
                for l in range(nrhs):
                    r[k, l] = r[k, l] - ztmp * r[np_, l]
        # ---- solve for last row
        for l in range(nrhs):
            r[nn - 1, l] = r[nn - 1, l] / z[nn - 1, nn - 1]
        # ---- back substitute everything
        for np_ in range(nn - 2, -1, -1):
            for k in range(np_ + 1, nn):
                for l in range(nrhs):
                    r[np_, l] = r[np_, l] - z[np_, k] * r[k, l]


# ------------------------------------------------------------------ HELICO
def helico_harmonics(double[:, ::1] sys_, double[::1] gam, int ii, int iimax, int nn,
                     bint duct, double[:, ::1] vg3, double[::1] va3):
    """The harmonic loop of induced.helico (SUBROUTINE HELICO): solves each harmonic's
    tridiagonal systems in turn and sums them into vg3 (VIND_GAM(3,I,J), Fortran indices,
    shape (ii+1, ii+1)) and va3 (VIND_ADW(3,I)), including the 1/n^2 series extrapolation.
    vg3 and va3 must be zero on entry. sys_ is SYS(1..4, 0..IIMAX) as in the Python."""
    cdef int i, j, jl, n, m
    cdef double rn, rn2, a, b, c, fsum, t
    cdef double[::1] ainv = memview_zeros(iimax + 1)
    cdef double[::1] csav = memview_zeros(iimax + 1)
    cdef double[::1] an_adw = memview_zeros(iimax + 1)
    cdef double[:, ::1] an_gam = memview_zeros2(iimax + 1, ii + 1)
    with nogil:
        for n in range(2, nn + 1, 2):
            rn = <double>n
            rn2 = rn * rn
            # ---- inner BC dAn/dx = 0, interior equations, outer BC
            ainv[0] = 1.0
            csav[0] = -1.0
            for i in range(1, iimax):
                b = sys_[0, i]
                a = -(sys_[0, i] + sys_[1, i] * rn2 + sys_[2, i])
                c = sys_[2, i]
                ainv[i] = 1.0 / (a - b * csav[i - 1])
                csav[i] = c * ainv[i]
            b = sys_[0, iimax]
            ainv[iimax] = 1.0 / (1.0 - b * csav[iimax - 1])

            # ==== An and dAn(i)/dGam(j): right-hand sides
            for i in range(iimax + 1):
                for j in range(ii + 1):
                    an_gam[i, j] = 0.0
            for i in range(1, ii + 1):
                an_gam[i, 0] = sys_[1, i] * rn2 * gam[i - 1]
                an_gam[i, i] = sys_[1, i] * rn2
            # ---- forward elimination, only up to the nonzero elements
            for i in range(1, iimax + 1):
                jl = (i if i < ii else ii) + 1
                for j in range(jl):
                    an_gam[i, j] = (an_gam[i, j] - sys_[0, i] * an_gam[i - 1, j]) * ainv[i]
            for i in range(iimax - 1, -1, -1):
                for j in range(ii + 1):
                    an_gam[i, j] = an_gam[i, j] - csav[i] * an_gam[i + 1, j]

            # ==== dAn(i)/dAdw
            an_adw[0] = 0.0
            for i in range(1, iimax):
                an_adw[i] = sys_[3, i] * rn2 * (gam[i - 1] + an_gam[i, 0])
            an_adw[iimax] = 0.0
            for i in range(1, iimax + 1):
                an_adw[i] = (an_adw[i] - sys_[0, i] * an_adw[i - 1]) * ainv[i]
            for i in range(iimax - 1, -1, -1):
                an_adw[i] = an_adw[i] - csav[i] * an_adw[i + 1]

            # ---- sum potential harmonics into the swirl-circulation coefficients
            for i in range(1, ii + 1):
                for j in range(1, ii + 1):
                    vg3[i, j] = vg3[i, j] + an_gam[i, j]
            for i in range(1, ii + 1):
                vg3[i, i] = vg3[i, i] + 1.0
                va3[i] = va3[i] + an_adw[i]

        # ---- extrapolate the series assuming (An + Gam) ~ 1/n^2 (from the last harmonic)
        if not duct:
            fsum = 0.0
            for m in range(nn + 2, 4 * nn + 1, 2):
                t = <double>nn / <double>m
                fsum = fsum + t * t
            for i in range(1, ii + 1):
                for j in range(1, ii + 1):
                    vg3[i, j] = vg3[i, j] + an_gam[i, j] * fsum
            for i in range(1, ii + 1):
                vg3[i, i] = vg3[i, i] + 1.0 * fsum
                va3[i] = va3[i] + an_adw[i] * fsum


cdef double[::1] memview_zeros(Py_ssize_t n):
    cdef double[::1] a = cvarray(shape=(n,), itemsize=sizeof(double), format="d")
    a[:] = 0.0
    return a


cdef double[:, ::1] memview_zeros2(Py_ssize_t n, Py_ssize_t m):
    cdef double[:, ::1] a = cvarray(shape=(n, m), itemsize=sizeof(double), format="d")
    a[:, :] = 0.0
    return a
