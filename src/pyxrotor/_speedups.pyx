# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# cython: initializedcheck=False, nonecheck=False
"""The compiled core of pyxrotor: the OPER solver loop (APITER, APINIT, TPQ, SETXW, VCALC,
HELICO, with CSCALC, UVADD, SEVAL, GETCLCDCM, GETALF and CLCDCM as C helpers) and GAUSS.

Ported from XROTOR 7.55 with the Fortran's operations in the Fortran's order, so results
are bit-identical to it; the differential tests check this routine by routine against the
compiled Fortran. That needs IEEE arithmetic from the C compiler: setup.py builds this file
with -ffp-contract=off (no fused multiply-add) and without -ffast-math. cdivision=True gives
IEEE division (x/0 -> inf, 0/0 -> nan), as pyxrotor.fmath.div does in the Python modules;
MIN/MAX follow fmath's per-call-site NaN conventions. Before arithmetic, State scalars are
converted to C doubles: Cython does Python float arithmetic on object attributes, which
raises where the Fortran gives inf.

(Until commit 092d52e the Python modules also had pure-Python versions of these routines;
the two were compared bit for bit, including on non-finite inputs.)
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
# CLCDCM's compressibility drag model constants ("HHY 10/23/00")
cdef double CDMFACTOR = 10.0
cdef double CLMFACTOR = 0.25
cdef double MEXP = 3.0
cdef double CDMDD = 0.0020
cdef double CDMSTALL = 0.1000


ctypedef struct SecAero:
    double clift, cl_alf, cl_w
    bint stallf
    double cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w, clmax, clmin
    bint mach_limited


# AERODATA rows (0-based) of the section parameters, as aero.AERODATA_ROW - 1
DEF R_A0 = 0
DEF R_CLMAX = 1
DEF R_CLMIN = 2
DEF R_DCLDA = 3
DEF R_DCLDA_STALL = 4
DEF R_DCL_STALL = 5
DEF R_CDMIN = 6
DEF R_CLDMIN = 7
DEF R_DCDCL2 = 8
DEF R_CMCON = 9
DEF R_REREF = 10
DEF R_REXP = 11
DEF R_MCRIT = 12


def clcdcm(double alf, double w, double rey, double vel, double vso, sec):
    """aero.clcdcm (SUBROUTINE CLCDCM); returns the SectionAero fields as a tuple."""
    cdef double p[13]
    p[R_A0], p[R_CLMAX], p[R_CLMIN] = sec.a0, sec.clmax, sec.clmin
    p[R_DCLDA], p[R_DCLDA_STALL], p[R_DCL_STALL] = sec.dclda, sec.dclda_stall, sec.dcl_stall
    p[R_CDMIN], p[R_CLDMIN], p[R_DCDCL2], p[R_CMCON] = sec.cdmin, sec.cldmin, sec.dcdcl2, sec.cmcon
    p[R_REREF], p[R_REXP], p[R_MCRIT] = sec.reref, sec.rexp, sec.mcrit
    cdef SecAero r = clcdcm_c(alf, w, rey, vel, vso, p)
    return (r.clift, r.cl_alf, r.cl_w, r.stallf, r.cdrag, r.cd_alf, r.cd_w, r.cd_rey,
            r.cmom, r.cm_al, r.cm_w, r.clmax, r.clmin, r.mach_limited)


cdef SecAero clcdcm_c(double alf, double w, double rey, double vel, double vso,
                      const double* p) noexcept nogil:
    """CLCDCM for the parameters p[0..12] (AERODATA(1..13,N) order)."""
    cdef double a0 = p[R_A0], clmax = p[R_CLMAX], clmin = p[R_CLMIN]
    cdef double dclda = p[R_DCLDA], dclda_stall = p[R_DCLDA_STALL], dcl_stall = p[R_DCL_STALL]
    cdef double cdmin = p[R_CDMIN], cldmin = p[R_CLDMIN], dcdcl2 = p[R_DCDCL2], cmcon = p[R_CMCON]
    cdef double mcrit = p[R_MCRIT], reref = p[R_REREF], rexp = p[R_REXP]
    cdef SecAero out
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

    out.clift, out.cl_alf, out.cl_w, out.stallf = clift, cl_alf, cl_w, stallf
    out.cdrag, out.cd_alf, out.cd_w, out.cd_rey = cdrag, cd_alf, cd_w, cd_rey
    out.cmom, out.cm_al, out.cm_w = cmom, cm_al, cm_w
    out.clmax, out.clmin, out.mach_limited = clmax, clmin, mach_limited
    return out


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


# =================================================================== the solver loop
# State access: s is a pyxrotor.state.State. Its arrays are C-contiguous NumPy arrays,
# taken as memoryviews (stall through a uint8 view, iaero as int64); its scalars are read
# once and written back where the Python writes them, so the final State is the same.
import numpy as np
from libc.math cimport atan2 as c_atan2, cos as c_cos, sin as c_sin
from libc.stdint cimport int64_t

cdef double PI_X = 3.141592654             # pyxrotor.constants.PI (XROTOR.INC's PI)


cdef inline double f_div(double a, double b) noexcept nogil:
    return a / b                          # IEEE, = fmath.div


cdef inline double f_fmin2(double a, double b) noexcept nogil:
    # fmath.fmin2: a if a < b else b (NaN -> b)
    return a if a < b else b


cdef inline double f_fmax2(double a, double b) noexcept nogil:
    return a if a > b else b


# ------------------------------------------------------------------ SEVAL, UVADD
cdef int interval(double ss, const double* sv, int n) noexcept nogil:
    """numerics._interval (SEVAL/SEVLIN bisection)."""
    cdef int ilow = 0, i = n - 1, imid
    while i - ilow > 1:
        imid = (i + ilow + 2) // 2 - 1
        if ss < sv[imid]:
            i = imid
        else:
            ilow = imid
    return i


cdef double seval_c(double ss, const double* x, const double* xs, const double* sv,
                    int n) noexcept nogil:
    """numerics.seval (FUNCTION SEVAL)."""
    cdef int i
    cdef double ds, t, cx1, cx2
    if n == 1:
        return x[0]
    i = interval(ss, sv, n)
    ds = sv[i] - sv[i - 1]
    t = (ss - sv[i - 1]) / ds
    cx1 = ds * xs[i - 1] - x[i] + x[i - 1]
    cx2 = ds * xs[i] - x[i] + x[i - 1]
    return t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)


cdef struct Ctx:
    # scalars and array pointers a routine needs from the State
    int ii, nadd, nblds, naero
    bint duct, section_fallback
    double adv, adw, vel, vso, rad, rho, rmu, urduct, rake
    double* xi
    double* dxi
    double* ubody
    double* beta
    double* beta0
    double* ch
    double* gam
    double* cl
    double* cd
    double* cm
    double* re
    unsigned char* stall
    int64_t* iaero
    double* vind          # VIND(3, 100): vind[k*ix + i]
    double* vind_adw
    double* vind_gam      # VIND_GAM(3, 100, 100): [(k*ix + i)*ix + j]
    double* radd
    double* uadd
    double* uaddr
    double* vadd
    double* vaddr
    double* aerodata      # AERODATA(14, NAX): [row*nax + n]
    double* xiaero
    int ix, nax
    # swept-blade model (pyxrotor.sweep; not XROTOR): off = XROTOR's straight blade
    bint sweep
    double* sw            # (5, IX): s, s', n, n', BETA' at the stations: [k*ix + i]
    double* coslam        # State.sweep_cos: cos of the effective sweep, set by CSCALC


cdef double* ptr(arr) except NULL:
    """Data pointer of a C-contiguous float64 array (checked)."""
    cdef double[::1] v = arr.reshape(-1)
    return &v[0]


_FLOAT_ARRAYS = ("xi", "dxi", "ubody", "beta", "beta0", "ch", "gam", "cl", "cd", "cm", "re",
                 "vind", "vind_adw", "vind_gam", "radd", "uadd", "uaddr", "vadd", "vaddr",
                 "aerodata", "xiaero", "xv", "dq", "dgam", "dgamold", "ti_gam", "pi_gam",
                 "tw_gam", "pw_gam", "tv_gam", "pv_gam", "effp", "vwak", "xw", "dxw", "vw_adv",
                 "vw_adw", "dxw_adv", "dxw_adw", "vw_gam", "dxw_gam", "xw_gam", "xw_adv",
                 "xw_adw", "dtii", "dpii", "dtwi", "dpwi", "dtvi", "dpvi")


def state_ok(s):
    """True if the State's arrays have the layouts the compiled routines need (as State()
    creates them)."""
    for name in _FLOAT_ARRAYS:
        a = getattr(s, name)
        if a.dtype != np.float64 or not a.flags.c_contiguous:
            return False
    return (s.stall.dtype == np.bool_ and s.stall.flags.c_contiguous
            and s.iaero.dtype == np.int64 and s.iaero.flags.c_contiguous)


cdef void load_ctx(Ctx* c, s) except *:
    if not state_ok(s):
        raise TypeError("State arrays must be C-contiguous float64 (stall bool, iaero int64)")
    c.ii, c.nadd, c.nblds, c.naero = s.ii, s.nadd, s.nblds, s.naero
    c.duct, c.section_fallback = s.duct, s.fixes.section_fallback
    c.adv, c.adw, c.vel, c.vso = float(s.adv), float(s.adw), float(s.vel), float(s.vso)
    c.rad, c.rho, c.rmu = float(s.rad), float(s.rho), float(s.rmu)
    c.urduct, c.rake = float(s.urduct), float(s.rake)
    c.xi, c.dxi, c.ubody = ptr(s.xi), ptr(s.dxi), ptr(s.ubody)
    c.beta, c.beta0, c.ch, c.gam = ptr(s.beta), ptr(s.beta0), ptr(s.ch), ptr(s.gam)
    c.cl, c.cd, c.cm, c.re = ptr(s.cl), ptr(s.cd), ptr(s.cm), ptr(s.re)
    cdef unsigned char[::1] st = s.stall.view(np.uint8)
    c.stall = &st[0]
    cdef int64_t[::1] ia = s.iaero
    c.iaero = &ia[0]
    c.vind, c.vind_adw, c.vind_gam = ptr(s.vind), ptr(s.vind_adw), ptr(s.vind_gam)
    c.radd, c.uadd, c.uaddr = ptr(s.radd), ptr(s.uadd), ptr(s.uaddr)
    c.vadd, c.vaddr = ptr(s.vadd), ptr(s.vaddr)
    c.aerodata, c.xiaero = ptr(s.aerodata), ptr(s.xiaero)
    c.ix = s.vind.shape[1]
    c.nax = s.aerodata.shape[1]
    c.sweep = s.sweep is not None
    c.coslam = ptr(s.sweep_cos)
    c.sw = NULL
    if c.sweep:
        from .sweep import station_arrays
        if s.sweep_cos.size != c.ix or s.sweep_cos.dtype != np.float64:
            raise TypeError("State.sweep_cos must be float64 of length IX")
        s._sweep_sta = station_arrays(s)            # kept alive on the State
        c.sw = ptr(s._sweep_sta)


cdef inline double VG(Ctx* c, int k, int i, int j) noexcept nogil:
    return c.vind_gam[(k * c.ix + i) * c.ix + j]


cdef void uvadd_c(Ctx* c, double xiw, double* wa, double* wt) noexcept nogil:
    """wake.uvadd (SUBROUTINE UVADD)."""
    cdef double rdim
    cdef int n
    wa[0] = 0.0
    wt[0] = 0.0
    if c.nadd <= 1:
        return
    rdim = xiw * c.rad
    n = c.nadd
    if c.radd[0] <= rdim and rdim <= c.radd[n - 1]:
        wa[0] = seval_c(rdim, c.uadd, c.uaddr, c.radd, n) / c.vel
        wt[0] = seval_c(rdim, c.vadd, c.vaddr, c.radd, n) / c.vel


# ------------------------------------------------------------------ CSCALC
cdef struct Vel:
    double utot, wa, wt, vt, vt_adw, va, va_adw, vd, vd_adw, ci, ci_adv, ci_vt
    double si, si_va, w, w_adv, w_vt, w_va, phi, p_adv, p_vt, p_va


cdef Vel cscalc_c(Ctx* c, int i) noexcept nogil:
    """wake.cscalc (SUBROUTINE CSCALC)."""
    cdef Vel v
    cdef double uduct, vaduct_va, vd_va, xi, wsq
    v.vt = c.vind[2 * c.ix + i]
    v.vt_adw = c.vind_adw[2 * c.ix + i]
    v.va = c.vind[i]
    v.va_adw = c.vind_adw[i]
    uduct = 0.0
    vaduct_va = 1.0
    if c.duct:
        uduct = c.urduct - 1.0
        vaduct_va = 2.0 * c.urduct
    v.vd = v.va * (vaduct_va - 1.0)
    vd_va = vaduct_va - 1.0
    v.vd_adw = vd_va * v.va_adw
    xi = c.xi[i]
    v.utot = 1.0 + uduct + c.ubody[i]
    uvadd_c(c, xi, &v.wa, &v.wt)
    v.ci = xi / c.adv - v.wt - v.vt
    v.ci_adv = -(xi / (c.adv * c.adv))
    v.ci_vt = -1.0
    v.si = v.utot + v.wa + v.va + v.vd
    v.si_va = 1.0 + vd_va
    wsq = v.ci * v.ci + v.si * v.si
    v.w = f_sqrt(wsq)
    v.w_adv = v.ci * v.ci_adv / v.w
    v.w_vt = v.ci * v.ci_vt / v.w
    v.w_va = v.si * v.si_va / v.w
    v.phi = c_atan2(v.si, v.ci)
    v.p_adv = -(v.si * v.ci_adv) / wsq
    v.p_vt = -(v.si * v.ci_vt) / wsq
    v.p_va = v.ci * v.si_va / wsq
    if c.sweep:
        c.coslam[i] = sweep_cos_c(c, i, v.ci, v.si, v.w)
    return v


cdef double sweep_cos_c(Ctx* c, int i, double ci, double si, double w) noexcept nogil:
    """cos of the effective sweep at station i (pyxrotor.sweep): the quarter-chord line's
    tangent T against the relative velocity (0, CI, SI)."""
    cdef int ix = c.ix
    cdef double sv = c.sw[i], sp = c.sw[ix + i], nv = c.sw[2 * ix + i]
    cdef double np_ = c.sw[3 * ix + i], bp = c.sw[4 * ix + i]
    cdef double b = c.beta[i], cb = c_cos(b), sb = c_sin(b), xi = c.xi[i]
    cdef double y = sv * cb - nv * sb
    cdef double yp = sp * cb - sv * sb * bp - np_ * sb - nv * cb * bp
    cdef double zp = sp * sb + sv * cb * bp + np_ * cb - nv * sb * bp
    cdef double tt = yp - y / xi, tn
    cdef double sinl
    tn = c_sqrt(1.0 + tt * tt + zp * zp)
    sinl = (tt * ci + zp * si) / (tn * w)
    if not (sinl * sinl < 1.0):                    # non-finite or degenerate: no sweep
        return 1.0
    return c_sqrt(1.0 - sinl * sinl)


cdef SecAero swept_c(SecAero r, const double* p, double cl, double dan, double rey) noexcept nogil:
    """Simple sweep (pyxrotor.sweep): r is CLCDCM in the normal plane (at the normal angle
    and Mach number); returns streamwise coefficients CL = CLn cos^2, CD = CDf + (CDn - CDf)
    cos^3 (CDf = the friction part CDMIN (Re/REREF)^REXP), derivatives with respect to the
    streamwise angle (dan = d an / d a) and W."""
    cdef double c2 = cl * cl, c3 = c2 * cl, rc, rc_rey, cdf, cdf_rey
    cdef SecAero o = r
    if rey <= 0:
        rc = 1.0
        rc_rey = 0.0
    else:
        rc = c_pow(rey / p[R_REREF], p[R_REXP])
        rc_rey = p[R_REXP] / rey
    cdf = p[R_CDMIN] * rc
    cdf_rey = cdf * rc_rey
    o.clift = r.clift * c2
    o.cl_alf = r.cl_alf * dan * c2
    o.cl_w = r.cl_w * cl * c2
    o.clmax = r.clmax * c2
    o.clmin = r.clmin * c2
    o.cdrag = cdf + (r.cdrag - cdf) * c3
    o.cd_alf = r.cd_alf * dan * c3
    o.cd_w = r.cd_w * cl * c3
    o.cd_rey = cdf_rey + (r.cd_rey - cdf_rey) * c3
    o.cmom = r.cmom * c2
    o.cm_al = r.cm_al * dan * c2
    o.cm_w = r.cm_w * cl * c2
    return o


def cscalc(s, int i):
    """wake.cscalc, returning the Velocities fields as a tuple."""
    cdef Ctx c
    load_ctx(&c, s)
    cdef Vel v = cscalc_c(&c, i)
    return (v.utot, v.wa, v.wt, v.vt, v.vt_adw, v.va, v.va_adw, v.vd, v.vd_adw, v.ci,
            v.ci_adv, v.ci_vt, v.si, v.si_va, v.w, v.w_adv, v.w_vt, v.w_va, v.phi, v.p_adv,
            v.p_vt, v.p_va)


# ------------------------------------------------------------------ GETCLCDCM, GETALF
cdef struct StaAero:
    double clift, cl_alf, cl_w, clmax, clmin, dcl_stall
    bint stallf
    double cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w


cdef inline void section_params(Ctx* c, int n, double* p) noexcept nogil:
    cdef int k
    for k in range(13):
        p[k] = c.aerodata[k * c.nax + n]


cdef StaAero getclcdcm_c(Ctx* c, int is_, double alf, double w, double rey) noexcept nogil:
    """aero.getclcdcm (SUBROUTINE GETCLCDCM), including 7.55's section fallback (and the
    Fixes.section_fallback correction); writes IAERO like the Python."""
    cdef int naero = c.naero, n, k
    cdef bint jumped
    cdef double xis = c.xi[is_], xisect1, frac
    cdef double p1[13]
    cdef double p2[13]
    cdef SecAero r1, r2
    cdef StaAero o
    n = <int>c.iaero[is_]
    if (n < 0 or n >= naero) and c.section_fallback:
        n = 0
        for k in range(naero):
            if c.xiaero[k] <= xis:
                n = k
        c.iaero[is_] = n
    if n < 0 or n >= naero:
        jumped = False
        if naero > 1:
            for k in range(naero):
                n = k
                if c.xiaero[k] <= xis:
                    c.iaero[is_] = k
                else:
                    jumped = True
                    break
        if not jumped:
            n = 0
            c.iaero[is_] = n
    # swept blade (pyxrotor.sweep): the sections are evaluated in the normal plane
    cdef double cl = 1.0, dan = 1.0, alfn = alf, wn = w, ca, sa
    if c.sweep:
        cl = c.coslam[is_]
        ca = c_cos(alf)
        sa = c_sin(alf)
        alfn = c_atan2(sa, ca * cl)
        dan = cl / (cl * cl * ca * ca + sa * sa)
        wn = w * cl
    section_params(c, n, p1)
    xisect1 = c.xiaero[n]
    r1 = clcdcm_c(alfn, wn, rey, c.vel, c.vso, p1)
    if c.sweep:
        r1 = swept_c(r1, p1, cl, dan, rey)
    o.clift, o.cl_alf, o.cl_w, o.stallf = r1.clift, r1.cl_alf, r1.cl_w, r1.stallf
    o.cdrag, o.cd_alf, o.cd_w, o.cd_rey = r1.cdrag, r1.cd_alf, r1.cd_w, r1.cd_rey
    o.cmom, o.cm_al, o.cm_w = r1.cmom, r1.cm_al, r1.cm_w
    o.clmax, o.clmin, o.dcl_stall = r1.clmax, r1.clmin, p1[R_DCL_STALL] * (cl * cl)
    if n < naero - 1:
        section_params(c, n + 1, p2)
        frac = (xis - xisect1) / (c.xiaero[n + 1] - xisect1)
        r2 = clcdcm_c(alfn, wn, rey, c.vel, c.vso, p2)
        if c.sweep:
            r2 = swept_c(r2, p2, cl, dan, rey)
        o.stallf = o.stallf or r2.stallf
        o.clift = (1.0 - frac) * o.clift + frac * r2.clift
        o.cl_alf = (1.0 - frac) * o.cl_alf + frac * r2.cl_alf
        o.cl_w = (1.0 - frac) * o.cl_w + frac * r2.cl_w
        o.clmax = (1.0 - frac) * o.clmax + frac * r2.clmax
        o.clmin = (1.0 - frac) * o.clmin + frac * r2.clmin
        o.dcl_stall = (1.0 - frac) * o.dcl_stall + frac * (p2[R_DCL_STALL] * (cl * cl))
        o.cmom = (1.0 - frac) * o.cmom + frac * r2.cmom
        o.cm_al = (1.0 - frac) * o.cm_al + frac * r2.cm_al
        o.cm_w = (1.0 - frac) * o.cm_w + frac * r2.cm_w
        o.cdrag = (1.0 - frac) * o.cdrag + frac * r2.cdrag
        o.cd_alf = (1.0 - frac) * o.cd_alf + frac * r2.cd_alf
        o.cd_w = (1.0 - frac) * o.cd_w + frac * r2.cd_w
        o.cd_rey = (1.0 - frac) * o.cd_rey + frac * r2.cd_rey
    return o


cdef inline StaAero station_aero(Ctx* c, int i, double alfa, double w, double rey) noexcept nogil:
    """oper._station_aero: GETCLCDCM, storing CL, CD, CM, STALL at station i."""
    cdef StaAero r = getclcdcm_c(c, i, alfa, w, rey)
    c.cl[i], c.cd[i], c.cm[i] = r.clift, r.cdrag, r.cmom
    c.stall[i] = 1 if r.stallf else 0
    return r


cdef struct Alpha:
    double alf, alf_cl, alf_w
    bint stallf, converged


cdef Alpha getalf_c(Ctx* c, int is_, double clift, double w) noexcept nogil:
    """aero.getalf (SUBROUTINE GETALF): 10 Newton steps from section 1's zero-lift angle."""
    cdef Alpha a
    cdef StaAero r
    cdef double dalf
    cdef int it
    a.stallf = False
    a.alf = c.aerodata[R_A0 * c.nax + 0]
    a.alf_cl = 0.0
    a.alf_w = 0.0
    a.converged = False
    for it in range(10):
        r = getclcdcm_c(c, is_, a.alf, w, 0.0)
        a.stallf = r.stallf
        dalf = -((r.clift - clift) / r.cl_alf)
        a.alf = a.alf + dalf
        a.alf_cl = 1.0 / r.cl_alf
        a.alf_w = -(r.cl_w / r.cl_alf)
        if fabs(dalf) < 1.0e-5:
            a.converged = True
            return a
    return a


def getclcdcm(s, int is_, double alf, double w, double rey):
    """aero.getclcdcm on the State's sections; returns the StationAero fields."""
    cdef Ctx c
    load_ctx(&c, s)
    cdef StaAero r = getclcdcm_c(&c, is_, alf, w, rey)
    return (r.clift, r.cl_alf, r.cl_w, r.clmax, r.clmin, r.dcl_stall, r.stallf, r.cdrag,
            r.cd_alf, r.cd_w, r.cd_rey, r.cmom, r.cm_al, r.cm_w)


def getalf(s, int is_, double clift, double w):
    cdef Ctx c
    load_ctx(&c, s)
    cdef Alpha a = getalf_c(&c, is_, clift, w)
    return (a.alf, a.alf_cl, a.alf_w, a.stallf, a.converged)


# ------------------------------------------------------------------ TPQ
def tpq(s, int itype):
    """oper.tpq (SUBROUTINE TPQ(ITYPE)): thrust, power, torque and sensitivities."""
    cdef Ctx c
    load_ctx(&c, s)
    cdef int ii = c.ii, i, j
    cdef double[::1] ti_gam = s.ti_gam, pi_gam = s.pi_gam, tw_gam = s.tw_gam
    cdef double[::1] pw_gam = s.pw_gam, tv_gam = s.tv_gam, pv_gam = s.pv_gam
    cdef double[::1] effp = s.effp, vwak = s.vwak, xw = s.xw, dxw = s.dxw
    cdef double[::1] vw_adv = s.vw_adv, vw_adw = s.vw_adw, dxw_adv = s.dxw_adv, dxw_adw = s.dxw_adw
    cdef double[:, ::1] vw_gam = s.vw_gam, dxw_gam = s.dxw_gam
    cdef double[::1] dtii = s.dtii, dpii = s.dpii, dtwi = s.dtwi, dpwi = s.dpwi
    cdef double[::1] dtvi = s.dtvi, dpvi = s.dpvi
    cdef double tinv = 0.0, pinv = 0.0, twak = 0.0, pwak = 0.0, tvis = 0.0, pvis = 0.0
    cdef double ti_adv = 0.0, pi_adv = 0.0, ti_adw = 0.0, pi_adw = 0.0
    cdef double tw_adv = 0.0, pw_adv = 0.0, tw_adw = 0.0, pw_adw = 0.0
    cdef double tv_adv = 0.0, pv_adv = 0.0, tv_adw = 0.0, pv_adw = 0.0
    cdef double tv_dbe = 0.0, pv_dbe = 0.0
    cdef double cosr, blds, adv, adw, xi, bdx, xx, xx_adv, w, ci, si, phi, alfa
    cdef double al_gi, al_dbe, al_adv, al_vt, al_va, rey, cl_al, cl_w
    cdef double cl_gi, cl_dbe, cl_adv, cl_vt, cl_va, ch_gi, ch_adv, ch_vt, ch_va
    cdef double gam, cli, chnew, chw, ch, clw, re_w, re_ch, re_gi, re_adv, re_vt, re_va
    cdef double cd_al, cd_w, cd_re, cd_gi, cd_adv, cd_vt, cd_va, cd_dbe, cdi
    cdef double eff, effi, hwc, hwc_w, hwc_ch, cosrv, cosri
    cdef double dtv, dtv_cd, dtv_si, dtv_w, dtv_ch, dtv_gi, dtv_dbe, dtv_adv, dtv_vt, dtv_va
    cdef double dpv, dpv_cd, dpv_ci, dpv_w, dpv_ch, dpv_gi, dpv_dbe, dpv_adv, dpv_vt, dpv_va
    cdef double dti, dti_ci, dti_gi, dti_vt, dti_adv, dti_adw
    cdef double dpi, dpi_si, dpi_gi, dpi_xx, dpi_va, dpi_adv, dpi_adw
    cdef double vw, wa, wt, xwi, cw, cw_adv, cw_vw, dxwi
    cdef double dtw, dtw_gi, dtw_cw, dtw_dxw, dtw_vw, dtw_adv, dtw_adw
    cdef double dpw, dpw_si, dpw_gi, dpw_adv, dpw_va, dpw_adw
    cdef Vel v
    cdef StaAero r
    cdef Alpha a

    for j in range(ii):
        ti_gam[j] = 0.0
        pi_gam[j] = 0.0
        tw_gam[j] = 0.0
        pw_gam[j] = 0.0
        tv_gam[j] = 0.0
        pv_gam[j] = 0.0

    if itype < 1 or itype > 3:
        # the Python zeroes the totals, then raises at the first station
        s.tinv = s.pinv = s.twak = s.pwak = s.tvis = s.pvis = 0.0
        s.ti_adv = s.pi_adv = s.ti_adw = s.pi_adw = 0.0
        s.tw_adv = s.pw_adv = s.tw_adw = s.pw_adw = 0.0
        s.tv_adv = s.pv_adv = s.tv_adw = s.pv_adw = 0.0
        s.tv_dbe = s.pv_dbe = 0.0
        if ii > 0:
            raise ValueError(f"TPQ: ITYPE must be 1, 2 or 3, not {itype}")

    cosr = c_cos(c.rake)
    blds = <double>c.nblds
    adv, adw = c.adv, c.adw

    with nogil:
        for i in range(ii):
            xi = c.xi[i]
            bdx = blds * c.dxi[i]
            xx = xi / adv
            xx_adv = -(xx / adv)

            v = cscalc_c(&c, i)
            w, ci, si, phi = v.w, v.ci, v.si, v.phi
            alfa = c.beta[i] - phi

            if itype == 1:
                alfa = c.beta[i] - phi
                al_gi = 0.0
                al_dbe = 1.0
                al_adv = -v.p_adv
                al_vt = -v.p_vt
                al_va = -v.p_va
                rey = c.ch[i] * fabs(w) * c.rho * c.vel * c.rad / c.rmu
                r = station_aero(&c, i, alfa, w, rey)
                cl_al, cl_w = r.cl_alf, r.cl_w
                cl_gi = cl_al * al_gi
                cl_dbe = cl_al * al_dbe
                cl_adv = cl_al * al_adv + cl_w * v.w_adv
                cl_vt = cl_al * al_vt + cl_w * v.w_vt
                cl_va = cl_al * al_va + cl_w * v.w_va
                ch_gi = 0.0
                ch_adv = 0.0
                ch_vt = 0.0
                ch_va = 0.0
            elif itype == 2:
                a = getalf_c(&c, i, c.cl[i], w)
                alfa = a.alf
                c.stall[i] = 1 if a.stallf else 0
                al_gi = 0.0
                al_dbe = 0.0
                al_adv = a.alf_w * v.w_adv
                al_vt = a.alf_w * v.w_vt
                al_va = a.alf_w * v.w_va
                cl_gi = 0.0
                cl_dbe = 0.0
                cl_adv = 0.0
                cl_vt = 0.0
                cl_va = 0.0
                gam = c.gam[i]
                cli = c.cl[i]
                chnew = 2.0 * gam / (w * cli)
                if chnew <= 0.0:
                    if i == 0:
                        c.ch[i] = c.ch[i + 1]
                    elif i == ii - 1:
                        c.ch[i] = c.ch[i - 1]
                    else:
                        c.ch[i] = 0.5 * (c.ch[i - 1] + c.ch[i + 1])
                    ch_gi = 0.0
                    ch_adv = 0.0
                    ch_vt = 0.0
                    ch_va = 0.0
                else:
                    c.ch[i] = 2.0 * gam / (w * cli)
                    ch_gi = 2.0 / (w * cli)
                    chw = -(c.ch[i] / w)
                    ch_adv = chw * v.w_adv
                    ch_vt = chw * v.w_vt
                    ch_va = chw * v.w_va
                c.beta[i] = alfa + phi
                c.beta0[i] = c.beta[i]
            else:
                ch = c.ch[i]
                c.cl[i] = 2.0 * c.gam[i] / (w * ch)
                cl_gi = 2.0 / (w * ch)
                cl_dbe = 0.0
                clw = -(c.cl[i] / w)
                cl_adv = clw * v.w_adv
                cl_vt = clw * v.w_vt
                cl_va = clw * v.w_va
                a = getalf_c(&c, i, c.cl[i], w)
                alfa = a.alf
                c.stall[i] = 1 if a.stallf else 0
                al_gi = a.alf_cl * cl_gi
                al_dbe = a.alf_cl * cl_dbe
                al_adv = a.alf_cl * cl_adv + a.alf_w * v.w_adv
                al_vt = a.alf_cl * cl_vt + a.alf_w * v.w_vt
                al_va = a.alf_cl * cl_va + a.alf_w * v.w_va
                ch_gi = 0.0
                ch_adv = 0.0
                ch_vt = 0.0
                ch_va = 0.0
                c.beta[i] = alfa + phi
                c.beta0[i] = c.beta[i]

            ch = c.ch[i]
            c.re[i] = ch * fabs(w) * c.rho * c.vel * c.rad / c.rmu
            re_w = ch * c.rho * c.vel * c.rad / c.rmu
            re_ch = fabs(w) * c.rho * c.vel * c.rad / c.rmu
            re_gi = re_ch * ch_gi
            re_adv = re_ch * ch_adv + re_w * v.w_adv
            re_vt = re_ch * ch_vt + re_w * v.w_vt
            re_va = re_ch * ch_va + re_w * v.w_va

            r = station_aero(&c, i, alfa, w, c.re[i])
            cd_al, cd_w, cd_re = r.cd_alf, r.cd_w, r.cd_rey
            cd_gi = cd_al * al_gi + cd_re * re_gi
            cd_adv = cd_al * al_adv + cd_re * re_adv + cd_w * v.w_adv
            cd_vt = cd_al * al_vt + cd_re * re_vt + cd_w * v.w_vt
            cd_va = cd_al * al_va + cd_re * re_va + cd_w * v.w_va
            cd_dbe = cd_al * al_dbe

            cli, cdi = c.cl[i], c.cd[i]
            eff = ((cli * ci - cdi * si) / (cdi * ci + cli * si)) / xx
            eff = eff * cosr
            effi = ci / (si * xx)
            effi = effi * cosr
            effp[i] = eff / effi

            hwc = 0.5 * w * ch
            hwc_w = 0.5 * ch
            hwc_ch = 0.5 * w

            cosrv = 1.0
            dtv = -hwc * cdi * si * bdx * cosrv
            dtv_cd = -hwc * si * bdx * cosrv
            dtv_si = -hwc * cdi * bdx * cosrv
            dtv_w = -hwc_w * cdi * si * bdx * cosrv
            dtv_ch = -hwc_ch * cdi * si * bdx * cosrv
            dtv_gi = dtv_cd * cd_gi + dtv_ch * ch_gi
            dtv_dbe = dtv_cd * cd_dbe
            dtv_adv = dtv_cd * cd_adv + dtv_ch * ch_adv + dtv_w * v.w_adv
            dtv_vt = dtv_cd * cd_vt + dtv_ch * ch_vt + dtv_w * v.w_vt
            dtv_va = dtv_cd * cd_va + dtv_ch * ch_va + dtv_si * v.si_va + dtv_w * v.w_va

            tvis = tvis + dtv
            tv_adv = tv_adv + dtv_adv
            tv_adw = dtv_vt * v.vt_adw + dtv_va * v.va_adw        # 7.55: not accumulated
            tv_dbe = tv_dbe + dtv_dbe
            tv_gam[i] = tv_gam[i] + dtv_gi
            for j in range(ii):
                tv_gam[j] = tv_gam[j] + dtv_vt * VG(&c, 2, i, j) + dtv_va * VG(&c, 0, i, j)

            dpv = hwc * cdi * ci * bdx * xx
            dpv_cd = hwc * ci * bdx * xx
            dpv_ci = hwc * cdi * bdx * xx
            dpv_w = hwc_w * cdi * ci * bdx * xx
            dpv_ch = hwc_ch * cdi * ci * bdx * xx
            dpv_gi = dpv_cd * cd_gi + dpv_ch * ch_gi
            dpv_dbe = dpv_cd * cd_dbe
            dpv_adv = (dpv_cd * cd_adv + dpv_ch * ch_adv + dpv_ci * v.ci_adv + dpv_w * v.w_adv
                       + hwc * cdi * ci * bdx * xx_adv)
            dpv_vt = dpv_cd * cd_vt + dpv_ch * ch_vt + dpv_ci * v.ci_vt + dpv_w * v.w_vt
            dpv_va = dpv_cd * cd_va + dpv_ch * ch_va + dpv_w * v.w_va

            pvis = pvis + dpv
            pv_adv = pv_adv + dpv_adv
            pv_adw = dpv_vt * v.vt_adw + dpv_va * v.va_adw        # 7.55: not accumulated
            pv_dbe = pv_dbe + dpv_dbe
            pv_gam[i] = pv_gam[i] + dpv_gi
            for j in range(ii):
                pv_gam[j] = pv_gam[j] + dpv_vt * VG(&c, 2, i, j) + dpv_va * VG(&c, 0, i, j)

            gam = c.gam[i]
            cosri = 1.0
            dti = gam * ci * bdx * cosri
            dti_ci = gam * bdx * cosri
            dti_gi = ci * bdx * cosri
            dti_vt = dti_ci * v.ci_vt
            dti_adv = dti_ci * v.ci_adv
            dti_adw = dti_vt * v.vt_adw
            tinv = tinv + dti
            ti_adv = ti_adv + dti_adv
            ti_adw = ti_adw + dti_adw
            ti_gam[i] = ti_gam[i] + dti_gi
            for j in range(ii):
                ti_gam[j] = ti_gam[j] + dti_vt * VG(&c, 2, i, j)

            dpi = gam * si * bdx * xx
            dpi_si = gam * bdx * xx
            dpi_gi = si * bdx * xx
            dpi_xx = gam * si * bdx
            dpi_va = dpi_si * v.si_va
            dpi_adv = dpi_xx * xx_adv
            dpi_adw = dpi_va * v.va_adw
            pinv = pinv + dpi
            pi_adv = pi_adv + dpi_adv
            pi_adw = pi_adw + dpi_adw
            pi_gam[i] = pi_gam[i] + dpi_gi
            for j in range(ii):
                pi_gam[j] = pi_gam[j] + dpi_va * VG(&c, 0, i, j)

            vw = vwak[i]
            uvadd_c(&c, xi, &wa, &wt)
            xwi = xw[i]
            cw = xwi / adv - wt - vw
            cw_adv = -(xwi / (adv * adv))
            cw_vw = -1.0

            dxwi = dxw[i]
            dtw = gam * cw * blds * dxwi
            dtw_gi = cw * blds * dxwi
            dtw_cw = gam * blds * dxwi
            dtw_dxw = gam * cw * blds
            dtw_vw = dtw_cw * cw_vw
            dtw_adv = dtw_cw * cw_adv + dtw_vw * vw_adv[i] + dtw_dxw * dxw_adv[i]
            dtw_adw = dtw_vw * vw_adw[i] + dtw_dxw * dxw_adw[i]
            twak = twak + dtw
            tw_adv = tw_adv + dtw_adv
            tw_adw = tw_adw + dtw_adw
            tw_gam[i] = tw_gam[i] + dtw_gi
            for j in range(ii):
                tw_gam[j] = tw_gam[j] + dtw_vw * vw_gam[i, j] + dtw_dxw * dxw_gam[i, j]

            dpw = gam * si * bdx * xi / adv
            dpw_si = gam * bdx * xi / adv
            dpw_gi = si * bdx * xi / adv
            dpw_adv = -(dpw / adv)
            dpw_va = dpw_si * v.si_va
            dpw_adw = dpw_va * v.va_adw
            pwak = pwak + dpw
            pw_adv = pw_adv + dpw_adv
            pw_adw = pw_adw + dpw_adw
            pw_gam[i] = pw_gam[i] + dpw_gi
            for j in range(ii):
                pw_gam[j] = pw_gam[j] + dpw_va * VG(&c, 0, i, j)

            dtii[i] = dti / bdx
            dpii[i] = dpi / bdx
            dtwi[i] = dtw / bdx
            dpwi[i] = dpw / bdx
            dtvi[i] = dtv / bdx
            dpvi[i] = dpv / bdx

    s.tinv, s.pinv, s.twak, s.pwak, s.tvis, s.pvis = tinv, pinv, twak, pwak, tvis, pvis
    s.ti_adv, s.pi_adv, s.ti_adw, s.pi_adw = ti_adv, pi_adv, ti_adw, pi_adw
    s.tw_adv, s.pw_adv, s.tw_adw, s.pw_adw = tw_adv, pw_adv, tw_adw, pw_adw
    s.tv_adv, s.pv_adv, s.tv_adw, s.pv_adw = tv_adv, pv_adv, tv_adw, pv_adw
    s.tv_dbe, s.pv_dbe = tv_dbe, pv_dbe
    s.ttot = twak + tvis
    s.ptot = pwak + pvis
    s.qtot = s.ptot * s.adv


# ------------------------------------------------------------------ VCALC, SETXW
def vcalc(s):
    """induced.vcalc (SUBROUTINE VCALC): VIND(k,I) = sum over J of VIND_GAM(k,I,J)*GAM(J),
    accumulated in J order from 0."""
    cdef Ctx c
    load_ctx(&c, s)
    cdef int ii = c.ii, k, i, j
    cdef double acc
    with nogil:
        for k in range(3):
            for i in range(ii):
                acc = 0.0
                for j in range(ii):
                    acc = acc + VG(&c, k, i, j) * c.gam[j]
                c.vind[k * c.ix + i] = acc


def setxw(s):
    """wake.setxw (SUBROUTINE SETXW): equivalent-prop radii and swirl with sensitivities."""
    cdef Ctx c
    load_ctx(&c, s)
    cdef int ii = c.ii, i, j, it
    cdef double[::1] xw_adv = s.xw_adv, xw_adw = s.xw_adw, vwak = s.vwak, vw_adv = s.vw_adv
    cdef double[::1] vw_adw = s.vw_adw, xw = s.xw, dxw = s.dxw, dxw_adv = s.dxw_adv
    cdef double[::1] dxw_adw = s.dxw_adw
    cdef double[:, ::1] xw_gam = s.xw_gam, vw_gam = s.vw_gam, dxw_gam = s.dxw_gam
    cdef double[::1] xwm_gam = np.zeros(ii)
    cdef double[::1] z_gam = np.zeros(ii)
    cdef double adv = c.adv, adw = c.adw, urduct = c.urduct
    cdef double xwm = float(s.xw0), xwm_adv = 0.0, xwm_adw = 0.0
    cdef double xi, xdx, vt, wa, si, dxwo, xwo, vw, vw_xwo, utotw, vaw, vaw_vw, vaw_xwo
    cdef double sw, sw_xwo, rez, rez_xwo, delxwo, rlx = float(s.rlx), vw_vt, vaw_adw
    cdef double sw_vt, sw_adw, z_xwo, z_xwm, z_sw, z_si, z_vt, z_va, z_adw, z_adv
    cdef Vel v
    with nogil:
        for i in range(ii):
            xi = c.xi[i]
            xdx = xi * c.dxi[i]
            v = cscalc_c(&c, i)
            vt, wa, si = v.vt, v.wa, v.si

            dxwo = f_sqrt(xwm * xwm + 2.0 * v.utot * xdx) - xwm
            xwo = xwm + 0.5 * dxwo

            for it in range(30):
                vw = vt * xi / xwo
                vw_xwo = -(vw / xwo)
                utotw = urduct
                vaw = vw * xwo / adw
                vaw_vw = xwo / adw
                vaw_xwo = vw / adw + vaw_vw * vw_xwo
                sw = utotw + wa + vaw
                sw_xwo = vaw_xwo
                rez = sw * xwo * 2.0 * (xwo - xwm) - si * xdx
                rez_xwo = sw * 2.0 * (2.0 * xwo - xwm) + sw_xwo * xwo * 2.0 * (xwo - xwm)
                delxwo = -(rez / rez_xwo)
                rlx = 1.0
                if fabs(delxwo) > 0.2 * (xwo - xwm):
                    rlx = 0.2 * (xwo - xwm) / fabs(delxwo)
                xwo = xwo + rlx * delxwo
                if fabs(delxwo) < 1.0e-6:
                    break

            vw = vt * xi / xwo
            vw_xwo = -(vw / xwo)
            vw_vt = xi / xwo

            utotw = urduct
            vaw = vw * xwo / adw
            vaw_vw = xwo / adw
            vaw_xwo = vw / adw
            vaw_adw = -(vaw / adw)

            sw = utotw + wa + vaw
            sw_xwo = vaw_xwo + vaw_vw * vw_xwo
            sw_vt = vaw_vw * vw_vt
            sw_adw = vaw_adw

            z_xwo = sw * 2.0 * (2.0 * xwo - xwm)
            z_xwm = -(sw * xwo * 2.0)
            z_sw = xwo * 2.0 * (xwo - xwm)
            z_si = -xdx

            z_xwo = z_sw * sw_xwo + z_xwo
            z_vt = z_sw * sw_vt
            z_va = z_si * v.si_va
            z_adw = z_sw * sw_adw

            z_adv = z_xwm * xwm_adv
            z_adw = z_xwm * xwm_adw + z_vt * v.vt_adw + z_va * v.va_adw + z_adw
            for j in range(ii):
                z_gam[j] = z_xwm * xwm_gam[j] + z_vt * VG(&c, 2, i, j) + z_va * VG(&c, 0, i, j)

            xw_adv[i] = -(z_adv / z_xwo)
            xw_adw[i] = -(z_adw / z_xwo)
            for j in range(ii):
                xw_gam[i, j] = -(z_gam[j] / z_xwo)

            vwak[i] = vw
            vw_adv[i] = vw_xwo * xw_adv[i]
            vw_adw[i] = vw_xwo * xw_adw[i] + vw_vt * v.vt_adw
            for j in range(ii):
                vw_gam[i, j] = vw_xwo * xw_gam[i, j] + vw_vt * VG(&c, 2, i, j)

            xw[i] = xwo

            dxw[i] = 2.0 * (xwo - xwm)
            dxw_adv[i] = 2.0 * (xw_adv[i] - xwm_adv)
            dxw_adw[i] = 2.0 * (xw_adw[i] - xwm_adw)
            for j in range(ii):
                dxw_gam[i, j] = 2.0 * (xw_gam[i, j] - xwm_gam[j])

            xwm = 2.0 * xwo - xwm
            xwm_adv = 2.0 * xw_adv[i] - xwm_adv
            xwm_adw = 2.0 * xw_adw[i] - xwm_adw
            for j in range(ii):
                xwm_gam[j] = 2.0 * xw_gam[i, j] - xwm_gam[j]
    if ii > 0:
        s.rlx = rlx
    s.xwtip = xwm


# ------------------------------------------------------------------ HELICO (whole)
cdef double PI_EXACT = 4.0 * 0.7853981633974483  # 4*ATAN(1), set exactly below
import math as _math
PI_EXACT = 4.0 * _math.atan(1.0)


def helico(s, double[::1] gam):
    """induced.helico (SUBROUTINE HELICO). gam is induced._gam_extended(s): GAM and the
    COMMON arrays after it, which HELICO reads past the blade stations."""
    cdef Ctx c
    load_ctx(&c, s)
    cdef int idim = 150, nn = 128, ii = c.ii, nblds = c.nblds, i, j, k, iimax
    cdef double adw = c.adw, xi0, xitip, xinf, xfac, dx, adwinv, qbsq, xm, x0, xp, bfac, vsum, xi
    cdef double[::1] x = np.zeros(idim + 1)
    cdef double[:, ::1] sys_ = np.zeros((4, idim + 1))
    cdef double[:, ::1] vg3 = np.zeros((ii + 1, ii + 1))
    cdef double[::1] va3 = np.zeros(ii + 1)
    cdef double[::1] xv = s.xv
    for i in range(ii):
        x[1 + i] = c.xi[i]
    xi0 = xv[0]
    xitip = xv[ii]
    x[0] = 2.0 * xi0 - c.xi[0]
    x[ii + 1] = 2.0 * xitip - c.xi[ii - 1]
    if c.duct:
        iimax = ii + 1
    else:
        xinf = 4.0 * xitip
        x[ii + 2] = 2.0 * xitip - c.xi[ii - 2]
        x[ii + 3] = 2.0 * xitip - c.xi[ii - 3]
        x[ii + 4] = 2.0 * xitip - c.xi[ii - 4]
        xfac = (x[ii + 4] - x[ii + 3]) / (x[ii + 3] - x[ii + 2])
        dx = (x[ii + 4] - x[ii + 3]) * xfac
        iimax = idim                       # if the loop ends without reaching XINF
        for i in range(ii + 5, idim):
            x[i] = x[i - 1] + dx
            if x[i] >= xinf:
                iimax = i
                break
            dx = dx * xfac

    for k in range(3):
        for i in range(ii):
            for j in range(ii):
                c.vind_gam[(k * c.ix + i) * c.ix + j] = 0.0
            c.vind_adw[k * c.ix + i] = 0.0

    adwinv = 1.0 / (adw * adw)
    qbsq = 0.25 * <double>(nblds * nblds)
    for i in range(1, iimax):
        xm, x0, xp = x[i - 1], x[i], x[i + 1]
        sys_[0, i] = (x0 + xm) / (x0 - xm)
        sys_[2, i] = (xp + x0) / (xp - x0)
        sys_[1, i] = qbsq * (1.0 / x0 + x0 * adwinv) * (xp - xm)
        sys_[3, i] = qbsq * (-2.0 * x0 * adwinv / adw) * (xp - xm)
    sys_[0, iimax] = -1.0 if c.duct else 1.0

    helico_harmonics(sys_, gam, ii, iimax, nn, c.duct, vg3, va3)

    for i in range(1, ii + 1):
        bfac = <double>nblds / (2.0 * PI_EXACT * x[i])
        vg3[i, i] = vg3[i, i] + 0.5
        for j in range(1, ii + 1):
            vg3[i, j] = vg3[i, j] * bfac
        va3[i] = va3[i] * bfac

    for i in range(ii):
        for j in range(ii):
            c.vind_gam[(2 * c.ix + i) * c.ix + j] = vg3[i + 1, j + 1]
        c.vind_adw[2 * c.ix + i] = va3[i + 1]

    for i in range(ii):
        xi = c.xi[i]
        for j in range(ii):
            c.vind_gam[(0 * c.ix + i) * c.ix + j] = c.vind_gam[(2 * c.ix + i) * c.ix + j] * xi / adw
            c.vind_gam[(1 * c.ix + i) * c.ix + j] = 0.0
        vsum = 0.0
        for j in range(ii):
            vsum = vsum + c.gam[j] * c.vind_gam[(2 * c.ix + i) * c.ix + j]
        c.vind_adw[i] = c.vind_adw[2 * c.ix + i] * xi / adw - vsum * xi / (adw * adw)
        c.vind_adw[c.ix + i] = 0.0


# ------------------------------------------------------------------ APINIT
def apinit(s, gradmo, vcalc_fn, filter_fn):
    """oper.apinit (SUBROUTINE APINIT). gradmo, vcalc_fn and filter_fn are the (dispatching)
    Python routines it calls between its station loops."""
    cdef Ctx c
    cdef int ii = s.ii, i, it
    cdef double blds = <double>s.nblds, uduct, vaduct_va, tsum, xi, utot, wa, wt, si, ci, w
    cdef double phi, alfa, rey, gam, vt_gam, vt, vt_adw, va_gam, va, va_adw, ci_vt, si_va
    cdef double wsq, w_vt, w_va, p_vt, p_va, al_vt, al_va, ch, rez, z_cl, z_w, z_g, z_al
    cdef double z_vt, z_va, z_adw, delg, dcl, rlx, dclmax, rlxmin, g_adw, cosr, t_g, t_vt
    cdef double t_adw, vhsq, vhsq_t, root, dadw, adv
    cdef StaAero r
    s.dbeta = 0.0
    uduct = 0.0
    vaduct_va = 1.0
    if s.duct:
        uduct = s.urduct - 1.0
        vaduct_va = 2.0 * s.urduct
    s.adw = s.adv * (1.0 + uduct)

    load_ctx(&c, s)
    adv = c.adv
    tsum = 0.0
    for i in range(ii):
        xi = c.xi[i]
        utot = c.urduct + c.ubody[i]
        uvadd_c(&c, xi, &wa, &wt)
        si = utot + wa
        ci = xi / adv - wt
        w = f_sqrt(ci * ci + si * si)
        phi = c_atan2(si, ci)
        alfa = c.beta[i] - phi
        rey = c.ch[i] * fabs(w) * c.rho * c.vel * c.rad / c.rmu
        station_aero(&c, i, alfa, w, rey)
        c.gam[i] = 0.5 * c.cl[i] * w * c.ch[i]
        tsum = tsum + blds * c.gam[i] * ci * c.dxi[i]

    vhsq = 0.5 * tsum / PI_X
    vhsq = f_fmax2(vhsq, -0.25)
    s.adw = s.adv * 0.5 * (1.0 + f_sqrt(1.0 + 4.0 * vhsq))

    vcalc_fn(s)

    for it in range(10):
        gradmo(s)
        load_ctx(&c, s)                    # ADW changed
        tsum = 0.0
        t_adw = 0.0
        dclmax = 0.0
        rlxmin = 1.0
        rlx = float(s.rlx)
        for i in range(ii):
            xi = c.xi[i]
            gam = c.gam[i]
            vt_gam = VG(&c, 2, i, i)
            vt = vt_gam * gam
            vt_adw = c.vind_adw[2 * c.ix + i]
            va_gam = VG(&c, 0, i, i)
            va = va_gam * gam
            va_adw = c.vind_adw[i]

            uduct = 0.0
            vaduct_va = 1.0
            if c.duct:
                uduct = c.urduct - 1.0
                vaduct_va = 2.0 * c.urduct
            utot = 1.0 + uduct + c.ubody[i]
            uvadd_c(&c, xi, &wa, &wt)

            ci = xi / c.adv - wt - vt
            ci_vt = -1.0
            si = utot + wa + va * vaduct_va
            si_va = vaduct_va
            wsq = ci * ci + si * si
            w = f_sqrt(wsq)
            w_vt = ci * ci_vt / w
            w_va = si * si_va / w
            phi = c_atan2(si, ci)
            p_vt = -(si * ci_vt) / wsq
            p_va = ci * si_va / wsq

            alfa = c.beta[i] - phi
            al_vt = -p_vt
            al_va = -p_va

            ch = c.ch[i]
            rey = ch * fabs(w) * c.rho * c.vel * c.rad / c.rmu
            r = station_aero(&c, i, alfa, w, rey)

            rez = ch * c.cl[i] * w - 2.0 * gam
            z_cl = ch * w
            z_w = ch * c.cl[i]
            z_g = -2.0
            z_al = z_cl * r.cl_alf
            z_w = z_cl * r.cl_w + z_w
            z_vt = z_w * w_vt + z_al * al_vt
            z_va = z_w * w_va + z_al * al_va
            z_adw = z_vt * vt_adw + z_va * va_adw
            z_g = z_vt * vt_gam + z_va * va_gam + z_g

            delg = -(rez / z_g)
            dcl = 2.0 * delg / (ch * w)

            rlx = 1.0
            if rlx * fabs(dcl) > 0.2:
                if dcl != 0.0:
                    rlx = f_fmin(rlx, 0.2 / fabs(dcl))
            if fabs(dcl) > fabs(dclmax):
                dclmax = dcl
            if fabs(rlx) < rlxmin:
                rlxmin = rlx

            c.gam[i] = gam + rlx * delg
            g_adw = -(z_adw / z_g)

            cosr = 1.0
            gam = c.gam[i]
            tsum = tsum + blds * gam * ci * c.dxi[i] * cosr
            t_g = blds * ci * c.dxi[i] * cosr
            t_vt = blds * gam * ci_vt * c.dxi[i] * cosr
            t_adw = t_adw + (t_g + t_vt * vt_gam) * g_adw + t_vt * vt_adw
        if ii > 0:
            s.rlx = rlx

        vhsq = 0.5 * tsum / PI_X
        vhsq = f_fmax2(vhsq, -0.2499)
        vhsq_t = 0.5 / PI_X
        root = f_sqrt(1.0 + 4.0 * vhsq)
        adv = float(s.adv)
        rez = s.adw - adv * 0.5 * (1.0 + root)
        z_adw = 1.0 - adv / root * vhsq_t * t_adw
        dadw = -(rez / z_adw)
        dadw = f_fmin2(dadw, 10.0 * s.adw)
        dadw = f_fmax2(dadw, -0.9 * s.adw)
        s.dadw = dadw
        s.adw = s.adw + dadw

        if rlxmin < 0.2:
            filter_fn(s.gam, 0.2 * ii, ii)
        if fabs(dclmax) < 0.001:
            break


# ------------------------------------------------------------------ APITER
cdef void gauss1(double[:, ::1] z, double* r, int nn) noexcept nogil:
    """GAUSS with one right-hand side (the same operations as gauss())."""
    cdef int np_, np1, n, nx, k, l
    cdef double d, pivot, ztmp, temp
    for np_ in range(nn - 1):
        np1 = np_ + 1
        nx = np_
        for n in range(np1, nn):
            d = fabs(z[n, np_]) - fabs(z[nx, np_])
            if not (d < 0.0 or d == 0.0):
                nx = n
        pivot = 1.0 / z[nx, np_]
        z[nx, np_] = z[np_, np_]
        for l in range(np1, nn):
            temp = z[nx, l] * pivot
            z[nx, l] = z[np_, l]
            z[np_, l] = temp
        temp = r[nx] * pivot
        r[nx] = r[np_]
        r[np_] = temp
        for k in range(np1, nn):
            ztmp = z[k, np_]
            for l in range(np1, nn):
                z[k, l] = z[k, l] - ztmp * z[np_, l]
            r[k] = r[k] - ztmp * r[np_]
    r[nn - 1] = r[nn - 1] / z[nn - 1, nn - 1]
    for np_ in range(nn - 2, -1, -1):
        for k in range(np_ + 1, nn):
            r[np_] = r[np_] - z[np_, k] * r[k]


def apiter(s, int ispec, int icon, int iq, induced_model, vcalc_fn, setxw_fn, tpq_fn,
           filter_fn, seval_fn, deval_fn):
    """oper.apiter (SUBROUTINE APITER(ISPEC,ICON)). The callables are the (dispatching)
    Python routines it calls once per Newton iteration."""
    if ispec < 1 or ispec > 5:
        pass                               # raised where the Python raises, below
    cdef Ctx c
    cdef int ii = s.ii, k1 = ii, k2 = ii + 1, k3 = ii + 2, it, i, j, nitera = s.nitera, n
    cdef double[::1] clmax = np.zeros(ii), clmin = np.zeros(ii), dclstall = np.zeros(ii)
    cdef double[:, ::1] q = np.zeros((iq, iq))
    cdef double[::1] dq = s.dq, dgam_a = s.dgam, dgamold = s.dgamold
    cdef double[::1] tw_gam, pw_gam, tv_gam, pv_gam
    cdef double eps = 1.0e-07, advfact = 1.0
    cdef double d_rho = s.rho, d_vel = s.vel, d_rad = s.rad
    cdef double rhov2r2 = d_rho * (d_vel * d_vel) * (d_rad * d_rad)
    cdef double w, alfa, al_dbe, al_p, ch, rey, cli, z_cl, z_w, z_gi, z_vt, z_va, z_adv, z_dbe
    cdef double z_tw, z_pw, gresmx, fresmx, aresmx, adv, adw, rlx, dgam, dcl, dclmin, dclmax
    cdef double dcllim, dbet, dadv, dadw, rms, t_spec, q_spec, p_spec, p_spec_adv, rhov3r2
    cdef double rpm, rpm_adv, pspec_rpm, pspec_adv
    cdef bint free
    cdef double d_adwfctr, d_twak, d_pwak, d_tw_adv, d_pw_adv, d_tw_adw, d_pw_adw, d_tvis, d_pvis, d_tv_adv, d_pv_adv, d_tv_adw, d_pv_adw, d_tv_dbe, d_pv_dbe, d_tspec, d_qspec
    cdef Vel v
    cdef StaAero r
    tw_gam, pw_gam, tv_gam, pv_gam = s.tw_gam, s.pw_gam, s.tv_gam, s.pv_gam

    for it in range(1, max(nitera, 1) + 1):
        free = s.free
        if free or it == 1:
            induced_model(s)
        vcalc_fn(s)
        setxw_fn(s)
        tpq_fn(s, 1)
        load_ctx(&c, s)
        adv, adw = c.adv, c.adw
        d_adwfctr, d_twak, d_pwak = s.adwfctr, s.twak, s.pwak
        d_tw_adv, d_pw_adv, d_tw_adw, d_pw_adw = s.tw_adv, s.pw_adv, s.tw_adw, s.pw_adw
        d_tvis, d_pvis, d_tv_adv, d_pv_adv = s.tvis, s.pvis, s.tv_adv, s.pv_adv
        d_tv_adw, d_pv_adw, d_tv_dbe, d_pv_dbe = s.tv_adw, s.pv_adw, s.tv_dbe, s.pv_dbe
        d_tspec, d_qspec = s.tspec, s.qspec

        gresmx = 0.0
        fresmx = 0.0
        aresmx = 0.0
        for j in range(k3 + 1):
            q[k2, j] = 0.0

        if free:
            dq[k2] = d_adwfctr * adw * d_twak / d_pwak - adv * advfact
            z_tw = d_adwfctr * adw / d_pwak
            z_pw = -(d_adwfctr * adw * d_twak / (d_pwak * d_pwak))
            for j in range(ii):
                q[k2, j] = z_tw * tw_gam[j] + z_pw * pw_gam[j]
            q[k2, k1] = z_tw * d_tw_adv + z_pw * d_pw_adv - advfact
            q[k2, k2] = z_tw * d_tw_adw + z_pw * d_pw_adw + d_adwfctr * d_twak / d_pwak
            aresmx = f_fmax(aresmx, fabs(dq[k2] / adv))
        else:
            dq[k2] = 0.0
            q[k2, k2] = 1.0

        w = 0.0
        with nogil:
            for i in range(ii):
                v = cscalc_c(&c, i)
                w = v.w
                alfa = c.beta[i] - v.phi
                al_dbe = 1.0
                al_p = -1.0
                ch = c.ch[i]
                rey = ch * fabs(w) * c.rho * c.vel * c.rad / c.rmu
                r = station_aero(&c, i, alfa, w, rey)
                clmax[i], clmin[i], dclstall[i] = r.clmax, r.clmin, r.dcl_stall
                cli = c.cl[i]

                dq[i] = ch * cli * w - 2.0 * c.gam[i]
                z_cl = ch * w
                z_w = ch * cli
                z_gi = -2.0
                z_vt = z_cl * (r.cl_alf * al_p * v.p_vt + r.cl_w * v.w_vt) + z_w * v.w_vt
                z_va = z_cl * (r.cl_alf * al_p * v.p_va + r.cl_w * v.w_va) + z_w * v.w_va
                z_adv = z_cl * (r.cl_alf * al_p * v.p_adv + r.cl_w * v.w_adv) + z_w * v.w_adv
                z_dbe = z_cl * (r.cl_alf * al_dbe)

                for j in range(ii):
                    q[i, j] = z_vt * VG(&c, 2, i, j) + z_va * VG(&c, 0, i, j)
                q[i, i] = q[i, i] + z_gi
                q[i, k1] = z_adv
                q[i, k2] = z_vt * v.vt_adw + z_va * v.va_adw
                q[i, k3] = z_dbe
                gresmx = f_fmax(gresmx, fabs(dq[i] / (0.1 * w)))

        if ispec == 1:
            t_spec = d_tspec / rhov2r2
            dq[k1] = d_twak + d_tvis - t_spec
            for j in range(ii):
                q[k1, j] = tw_gam[j] + tv_gam[j]
            q[k1, k1] = d_tw_adv + d_tv_adv
            q[k1, k2] = d_tw_adw + d_tv_adw
            q[k1, k3] = d_tv_dbe
            fresmx = f_fmax(fresmx, fabs(dq[k1]))
        elif ispec == 2:
            q_spec = d_qspec / (d_rho * (d_vel * d_vel) * (d_rad * d_rad * d_rad))
            dq[k1] = (d_pwak + d_pvis) * adv - q_spec
            for j in range(ii):
                q[k1, j] = (pw_gam[j] + pv_gam[j]) * adv
            q[k1, k1] = (d_pw_adv + d_pv_adv) * adv + d_pwak + d_pvis
            q[k1, k2] = (d_pw_adw + d_pv_adw) * adv
            q[k1, k3] = d_pv_dbe * adv
            fresmx = f_fmax(fresmx, fabs(dq[k1]))
        elif ispec == 3:
            p_spec = (<double>s.pspec) / (d_rho * (d_vel * d_vel * d_vel) * (d_rad * d_rad))
            dq[k1] = d_pwak + d_pvis - p_spec
            for j in range(ii):
                q[k1, j] = pw_gam[j] + pv_gam[j]
            q[k1, k1] = d_pw_adv + d_pv_adv
            q[k1, k2] = d_pw_adw + d_pv_adw
            q[k1, k3] = d_pv_dbe
            fresmx = f_fmax(fresmx, fabs(dq[k1]))
        elif ispec == 4:
            dq[k1] = 0.0
            for j in range(ii):
                q[k1, j] = 0.0
            q[k1, k1] = 1.0
            q[k1, k2] = 0.0
            q[k1, k3] = 0.0
        elif ispec == 5:
            rhov3r2 = d_rho * (d_vel * d_vel * d_vel) * (d_rad * d_rad)
            p_spec = (<double>s.pspec) / rhov3r2
            p_spec_adv = 0.0
            if s.lpwrvar:
                rpm = d_vel / (d_rad * adv * PI_X / 30.0)
                rpm_adv = -(rpm / adv)
                n = s.npwrvar
                s.pspec = seval_fn(rpm, s.pwrvar, s.xpwrvar, s.rpmvar, n)
                pspec_rpm = deval_fn(rpm, s.pwrvar, s.xpwrvar, s.rpmvar, n)
                pspec_adv = pspec_rpm * rpm_adv
                p_spec = (<double>s.pspec) / rhov3r2
                p_spec_adv = pspec_adv / rhov3r2
            dq[k1] = d_pwak + d_pvis - p_spec
            for j in range(ii):
                q[k1, j] = pw_gam[j] + pv_gam[j]
            q[k1, k1] = d_pw_adv + d_pv_adv - p_spec_adv
            q[k1, k2] = d_pw_adw + d_pv_adw
            q[k1, k3] = d_pv_dbe
            fresmx = f_fmax(fresmx, fabs(dq[k1]))
        else:
            s.gresmx, s.fresmx, s.aresmx = gresmx, fresmx, aresmx
            raise ValueError(f"APITER: ISPEC must be 1-5, not {ispec}")
        s.gresmx, s.fresmx, s.aresmx = gresmx, fresmx, aresmx

        dq[k3] = 0.0
        for j in range(k3 + 1):
            q[k3, j] = 0.0
        if icon == 1:
            q[k3, k1] = 1.0
        if icon == 2:
            q[k3, k3] = 1.0

        gauss1(q, &dq[0], k3 + 1)

        rlx = 1.0
        if it <= 2:
            rlx = 0.2
        for i in range(ii):
            dgam_a[i] = -dq[i]
            dgam = dgam_a[i]
            cli = c.cl[i]
            dcl = 2.0 * dgam / (c.ch[i] * w)          # 7.55: the last station's W
            dclmin = f_fmax(1.5 * dclstall[i], fabs(cli - clmin[i]))
            dclmax = f_fmax(1.5 * dclstall[i], fabs(clmax[i] - cli))
            dcllim = f_fmin2(f_fmin2(dclmax, dclmin), 0.5)
            dcllim = f_fmax2(dcllim, 0.01)
            if rlx * fabs(dcl) > dcllim:
                rlx = f_fmin2(rlx, dcllim / fabs(dcl))
            if dgam * dgamold[i] < 0.0:
                if fabs(dgam) > 0.2 * fabs(dgamold[i]):
                    rlx = f_fmin2(rlx, 0.2)

        dadv = -dq[k1]
        dadw = -dq[k2]
        dbet = -dq[k3]
        s.dadv, s.dadw, s.dbet = dadv, dadw, dbet

        if nitera == 0:
            rlx = 0.0

        if rlx * dbet > 0.05:
            rlx = f_fmin(rlx, 0.05 / dbet)
        if rlx * dbet < -0.05:
            rlx = f_fmin(rlx, -0.05 / dbet)
        if rlx * dadv > 0.5 * adv:
            rlx = f_fmin(rlx, 0.5 * adv / dadv)
        if rlx * dadv < -0.3 * adv:
            rlx = f_fmin(rlx, -0.3 * adv / dadv)
        if rlx * dadw > 0.5 * adw:
            rlx = f_fmin(rlx, 0.5 * adw / dadw)
        if rlx * dadw < -0.3 * adw:
            rlx = f_fmin(rlx, -0.3 * adw / dadw)
        s.rlx = rlx

        rms = 0.0
        for i in range(ii):
            dgam = dgam_a[i]
            c.gam[i] = c.gam[i] + rlx * dgam
            c.beta[i] = c.beta[i] + rlx * dbet
            c.beta0[i] = c.beta0[i] + rlx * dbet
            rms = rms + dgam * dgam / (1.0 + 1.0 / (adv * adv))
            dgamold[i] = dgam
        s.dbeta = s.dbeta + rlx * dbet
        s.adv = adv + rlx * dadv
        s.adw = adw + rlx * dadw
        s.rms = f_sqrt(rms / <double>ii)

        if rlx < 0.2:
            filter_fn(s.gam, 0.2 * ii, ii)

        if s.rms <= eps:
            if s.free:
                induced_model(s)
            vcalc_fn(s)
            tpq_fn(s, 1)
            s.conv = True
            return
