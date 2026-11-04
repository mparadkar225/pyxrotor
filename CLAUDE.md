# pyxrotor — notes for agents working in this repo

Python port of XROTOR 7.55. The compiled 7.55 Fortran in `reference/` is the oracle; the
port is only as good as its differential tests against it.

## Commands

```bash
git submodule update --init      # reference sources
reference/build.sh               # build both reference libraries (~1 min)
pytest                           # everything, one worker per CPU (pytest-xdist)
pytest -n0 tests/... --capture=no < /dev/null   # one process, to see Fortran output
pytest tests/differential -q     # just the port-vs-Fortran tests
HYPOTHESIS_PROFILE=deep pytest tests/differential -W error::RuntimeWarning   # 20k cases/test
python tools/common_writes.py TPQ APITER    # COMMON variables a Fortran routine assigns
python tools/minmax_sites.py xoper --routine APITER   # how each MIN/MAX treats NaN
python tools/gen_layout.py                  # regenerate src/pyxrotor/_layout.py from XROTOR.INC
```

Run the deep profile after porting or changing a routine; the default 400 cases have
missed real mismatches (an overflow in a diverging GETALF only showed up on other seeds).
Tests that run whole solves use `conftest.EXPENSIVE` (1/8 of the profile's examples).

## Porting rules

1. **XROTOR 7.55 is the source of truth.** Port from `reference/xrotor-7.55/src/*.f`.
   xrotor-python is not an oracle (it fails its own test case and uses a different
   section model); use it only as a structural example.
2. **One routine at a time, each with a differential test** in `tests/differential/`.
   The pattern: build a realistic `State` (`rotors.perturbed_states()` perturbs a converged
   Rotol-propeller solution), copy it into the reference with `statesync.state_to_ref`,
   call the Fortran (`ref.call("NAME", ...)`) and the port on the same State, then
   `statesync.assert_states_match(s, ref)`: every COMMON variable must agree, not just the
   outputs you expected. Only TEMP_X, plotting and I/O-unit blocks are excluded.
3. **Before porting, run `tools/common_writes.py ROUTINE`.** Routines that include
   XROTOR.INC share every COMMON name, so "locals" such as RLX in SETXW are solver state.
   Every listed variable must be written by the port (or the omission explained).
4. **Tolerances**: everything ported so far matches with `rtol=0`, including the full
   Newton solver. Anything looser than bitwise needs a written reason.
5. **Match the arithmetic, not just the math.** Things that have already mattered:
   * `PI` is XROTOR's `3.141592654` (`pyxrotor.constants.PI`) — except in routines that do
     not include XROTOR.INC: GRADMO, HELICO, VRTXCO compute `4*ATAN(1)` (exact), and
     VORSEGVEL uses `PI4 = 12.56637062`. Check which PI a routine sees.
   * Use `math` (or `pyxrotor.fmath`) for scalar transcendentals: same libm as the Fortran.
     NumPy's vectorized sin/cos can differ by an ulp; NumPy +,-,*,/,sqrt are exact.
   * Write `x**2` as `x * x`, `x**3` as `x * x * x`; REAL exponents use `fmath.pow`.
   * Fortran evaluates `a*b*c`, `a/b*c`, `a+b+c` left to right; keep the grouping, and note
     `RHO*VEL**2*RAD**2` is `(rho*(v*v))*(r*r)`, not a flat chain.
   * Element-wise NumPy over independent elements is bit-identical (row operations in
     GAUSS, VIND_GAM updates over J). Sums are not: use a sequential loop or
     `np.add.accumulate` over `[0.0, *terms]` (the leading 0.0 gives Fortran's +0.0 for
     all-(-0.0) sums). Never `np.sum`, which sums pairwise.
   * Python raises where Fortran returns inf/NaN. Use `pyxrotor.fmath` (`div`, `pow`,
     `sqrt`, `exp`, `log`, `sin`, `cos`) wherever a zero divisor, overflow, or invalid
     argument is possible. Diverging Newton iterations reach these paths and tests cover them.
     `fmath.pow` follows C99 pow (tested against libm in `test_fmath.py`): `math.pow(0, y<0)`
     raises, but the Fortran gets the pole, +-inf (CLCDCM when REY/REREF underflows to 0).
   * MIN/MAX with NaN depend on the call site. Each compiles to MINSD/MAXSD, which return the
     *source* operand when either operand is NaN. A constant or memory operand is always the
     source; otherwise register allocation decides. `python tools/minmax_sites.py xoper
     --routine APINIT --context 8` lists each site with its source line and operands. Port
     NaN→first-argument sites with `fmath.fmin/fmax` and NaN→second-argument sites with
     `fmin2/fmax2`. For example, APINIT's `MAX(VHSQ,-0.25)` returns -0.25, while CLCDCM's
     `MAX(0.0,x)` returns 0.0. Pin each site with a non-finite test (`test_apinit_non_finite`,
     `test_clcdcm_non_finite`).
   * IFIX/INT of NaN, ±inf or anything outside the INTEGER*4 range gives -2**31 (`fmath.ifix`).
     INTEGER arithmetic wraps (`fmath.iwrap`). VRTXCO reaches this path when ADW is NaN or 0.
   * Convert inputs to `float` at public entry points: NumPy scalars arriving from arrays
     (e.g. `xi[i]`) silently change overflow and divide-by-zero behavior.
   * Indices are 0-based in Python. Fortran IAERO = N is `iaero = N - 1`; "unset" (0) is -1.
     `State.INDEX_VARS` lists such variables; statesync converts them.
6. **Known 7.55 bugs are reproduced first, then fixed.** Port the faithful behavior (mark it
   `# 7.55:` in the code), prove it matches, then add a fix behind an option, with a test for each.
   Options that act inside the solvers live in `State.fixes` (Python-only, not COMMON, all
   False by default). `Rotor(fix_all=True)` turns them all on. Test each fix against a
   patched shim or prepared reference inputs (`test_fixes.py`).
7. **Keep Fortran names in comments** next to their Python equivalents and cite the
   source routine in each docstring. Code with no effect on state (unused locals, prints,
   plotting) may be omitted; say so in the docstring.

## The compiled core

`src/pyxrotor/_speedups.pyx` (Cython, built by `setup.py` with `-ffp-contract=off`, no
fast-math; a C compiler is required) is the only implementation of the OPER solver loop:
APITER, APINIT, TPQ, SETXW, VCALC, HELICO, with CSCALC, UVADD, SEVAL, GETCLCDCM, GETALF and
CLCDCM as C helpers, and GAUSS. The Python functions of those names (oper.py, wake.py,
induced.py, aero.clcdcm, numerics.gauss) are thin wrappers. GRADMO, VRTXCO, FILTER, DESI,
the I/O and the new models are Python; the compiled APITER/APINIT call back into them.
* Edit those routines in the .pyx and rebuild: `python setup.py build_ext --inplace` (an
  editable install picks up the in-place build). The differential tests are the check.
* The porting rules apply unchanged in C: Fortran operation order, IEEE division
  (`cdivision=True` = `fmath.div`), fmath's MIN/MAX NaN conventions per call site, and
  sqrt/log of a negative giving +NaN as in fmath.
* Convert State scalars to C doubles before any arithmetic: Cython does Python float
  arithmetic on object attributes (`s.pwak`), which raises ZeroDivisionError where C gives
  inf. Arrays go through `load_ctx`, which requires the layouts State() creates
  (C-contiguous float64; stall bool; iaero int64) and raises TypeError otherwise.
* Until 092d52e there were also pure-Python versions, compared bit for bit with the
  compiled ones (that comparison found the `fmath.pow` pole bug); see git history.
* Rotol OPER RPM solve: 5.5 ms (the pure Python was 44 ms; the Fortran takes 1.8 ms).

## Reference harness facts

* `Ref755()` loads a private copy of the library, so instances never share state.
* `ref.call(...)` passes everything by reference and returns the converted arguments, so
  output arguments can be read back. Python bool → LOGICAL, int → INTEGER*4, float → REAL*8.
  COMMON arrays can be passed as arguments directly: `ref.call("HELICO", ..., ref["GAM"], ...)`.
* A Fortran `STOP` kills the Python process: run anything that might STOP in a subprocess.
  Under xdist (the default) a STOP crashes one worker and the test is reported failed. In
  a single process (`-n0`) STOP exits with status 0 before pytest prints its summary or
  writes a report, so a truncated run looks like a pass: a run without the final
  `N passed` line was cut short. CI also checks that `--junitxml` wrote its report
  (`test_stop_guard.py` pins both behaviors; a C exit handler calling back into Python
  segfaults after interpreter shutdown, so don't).
  Strategies must generate inputs the shims accept:
  `blades()` once drew stations from `hub`, below REFARBI's XI0 = (hub*R)/R by an ulp.
* Interactive routines (ASKR prompts) can't be called; transcribe them into
  `reference/shim/ref755_shim.f` statement for statement, with prompts turned into arguments,
  and keep fixed-form lines within 72 columns. OPER commands are there as REFVELO, REFRPM, ...
* `reference/build.sh` must be rerun after editing the shim. It replaces the library atomically,
  but a long test run started before a rebuild mixes old and new shims: restart it.
* Fortran output and list-directed input go through `pyxrotor.fortio`, which is tested against
  gfortran (`REFFMT`, `REFRDL`). gfortran rounds ties to even on the exact binary value,
  writes `0.10000+101` for 3-digit exponents, and doesn't write trailing `X` spacing at the
  end of a record. List-directed REAL*8 output (`fortio.write_list`, tested with `REFLDW`)
  gives each value a 26-column block with 17 significant digits: F form in 21 columns
  plus 5 blanks, otherwise 1P E form with a 3-digit exponent.
* Don't import a test module inside a `@given` test: its own `@given` decorators then run
  mid-test, and Hypothesis fails it as nested. Import at module top.
* A Fortran runtime error (a READ without `ERR=`/`END=`, a prompt hitting EOF) kills the test
  process while `Ref755` has stdout redirected, so pytest prints nothing. Rerun the test with
  `--capture=no < /dev/null` to see the message, and test such paths in a subprocess.

## Known issues in XROTOR 7.55 (all reproduced faithfully unless noted)

* **Stale aero-section indices after DESI.** If aero sections are defined before the first
  design, `IAERO` is computed against the pre-design radial grid (all zeros at start-up), and
  DESGEN never refreshes it. Outer stations then extrapolate past the last section and get a
  negative Cd. In the P-51D test (`test_p51d_stale_sections`), 7.55 reports 85.27%. The same
  blade re-analysed with the correct sections gives 79.80% and draws 5% more power; a
  design made with the correct sections gives 80.26%. ARBI and LOAD call SETIAERO, so they
  are not affected. The fix is `Rotor(refresh_design_sections=True)`, which runs SETIAERO
  after SETX in DESGEN. The shim's `REFDESGEN(OK, LFIXIA)` is the oracle for the fixed path.
* **DESI SETCLD uses the current (pre-DESGEN) XI grid**, so a linear CL (`EDIT CL`) is
  evaluated on the previous design's stations. On a first design every station is at
  XI = 0, so the whole blade gets the root CL. DESGEN's stall check then tests those
  values. Fix: `Rotor(fix_linear_cl=True)` (`State.fixes.linear_cl`) re-evaluates the last
  SETCLD (`State.cldes_linear`) on the design grid and checks stall there. The oracle is
  `REFDESGENX(OK, LFIXIA, LFIXCL, CL1, CL2)`. On a 3-blade test design with CL 0.7→0.4,
  7.55 gives 86.4% (CL 0.7 everywhere) and the fix gives 83.6%.
* **DESI `ITER n` (with an inline argument) sets NITERA**, the OPER limit. Only the prompted
  form sets the design limit NITERD.
* **DEITER's ducted branch reads GAMDES**, a local that is never initialized, so a ducted
  design is undefined in 7.55. pyxrotor raises NotImplementedError. DESINI's ZETA is likewise
  undefined if neither DEST nor DESP is set; pyxrotor raises ValueError.
* **GETCLCDCM's fallback picks the wrong section.** When a station's IAERO is invalid (as
  INIT leaves it, 0 everywhere, until SETIAERO runs, or after sections are deleted), the
  recovery search exits with N at the first section *outboard* of the station instead of
  the inner bounding one. It stores the inner index, so only that first lookup is wrong. A
  station outboard of the last section finds no bounding pair and gets section 1, the root
  airfoil, permanently. Fix: `Rotor(fix_section_fallback=True)`
  (`State.fixes.section_fallback`) picks the SETIAERO index. The oracle is the unmodified
  routine with those IAERO entries pre-filled (`test_fixes.py`).
* **GETALF starts every station from section 1's zero-lift angle**, not the station's own.
* **TPQ: TV_ADW and PV_ADW are overwritten at each station, not accumulated**, so the
  viscous thrust/power sensitivity to ADW in APITER's Jacobian comes from the tip station only.
* **APITER's CL-change limiter uses the tip station's W for every station** (the leftover
  value of W from the residual loop).
* **HELICO reads GAM beyond the II blade stations** (up to IIMAX-1, the far-field grid),
  i.e. stale values; past IX it reads into the next COMMON arrays (DTII, ...). The port
  reproduces this via `induced._gam_extended`.
* **SETXW overwrites the COMMON relaxation factor RLX** with its Newton step factor.
* ADWCORR is dead code (never called); APITER's 0.75R ADVFACT is computed then set to 1.0.
* **SETCAS's increments are COMMON.** The local names DADV and DBET are the solver's Newton
  steps, so ASEQ leaves its advance-ratio step in DADV and BSEQ its blade-angle step (in
  degrees) in DBET. The raw step is stored before the `first == last` early return, and the
  signed one after it.
* **CASPLT solves each non-engine case twice** (APER(4,2) first, then the case type's own
  solve). Afterwards it restores the original operating point and solves it again.
* **ADDC and GETCAS don't set CASPAR(0)**, the case-type code, so those cases keep whatever
  code was stored there before.
* **CGET loses results.** GETCAS reads records into 128 characters, but a SHOCAS line with
  results is 159 characters long, so power, thrust, torque and efficiency are dropped
  (fix: `read_cases(full_lines=True)`).
* **GETFLT over-counts a comma after a blank** (`1,2 ,3`), and the READ that follows has no
  `END=`, so XROTOR stops with "Fortran runtime error: End of file". The port raises EOFError
  (`test_getflt_aborts`).
* **SAVE writes LVDUCT where LOAD reads FREE**, so after LOAD the free-wake flag means "ADW
  differs from ADV*URDUCT" (fix: `load(keep_free_wake=True)`). SAVE stores BETA0, not BETA.
* **LOAD stops on a line reading `END` in place of the URDuct record**: it then reads the
  URDUCT value as the integer NADD (Fortran runtime error). The port raises LoadError.
* **FNAME is CHARACTER*80**: XROTOR silently truncates longer file paths (pyxrotor doesn't).
  SAVE asks on stdin before overwriting, and LOAD asks on stdin when the file's station
  count differs from II. In tests, use short, fresh paths and feed prompts with
  `test_restart.fortran_stdin`.
* **NACELL's sources are shifted one interval downstream.** The "HHY FIX" sets source
  strengths from area changes, but SOURCE(i) = A(i) - A(i-1) is placed at the midpoint of
  interval i. So the first interval gets no source and the last interval's change is
  dropped. On a P-51-like spinner this makes the root UBODY -0.199 instead of -0.146 (fix:
  `nacelle(..., fix_source_positions=True)`; the shim's `REFNACELL(U, LFIXSRC)` is the
  oracle). NACE CONST also sets XW0 = 0. A bad number in the file zeroes UBODY but leaves
  VBODY and XW0 unchanged. An empty file stops XROTOR (READ without END=).
* **OPER ADW sets the factor ADWFCTR**, not ADW: the free-wake condition is
  ADWFCTR·ADW·Twake/Pwake = ADV. It has no effect with a rigid wake, and a free-wake DESI
  design resets it to 1.
* **VPUT/VGET slipstream profiles.** SAVVEL writes the far-downstream averaged induced
  velocities (2x the disk values, from Va = Vt*r/ADW) at the *disk* radii, not the
  contracted wake radii. VGET's weights (axial 0 to 1, swirl 0 or +-1) are the user's job;
  -1 is for a counter-rotating rotor. UVADD imposes nothing outside the profile's radius
  range.
* **LOAD doesn't spline the slipstream it reads.** UADDR/VADDR stay stale (zero in a fresh
  session), so a profile from a restart file is interpolated with zero slopes. Fix:
  `Rotor(fix_slipstream_splines=True)` (`State.fixes.slipstream_splines`).
* **DESI fails** (GETALF inversion, NaN wake) when the helical tip Mach number reaches about 1.
* The generic compressibility model (dCD = K (M − Mcrit)^3, local Mach capped at 0.99) is
  very sensitive to Mcrit for transonic tips.

## New models (beyond XROTOR)

Code that isn't a port lives in its own modules (`field.py`, ...), says so in its
docstring, and is validated in `tests/model/` by exact limits and independent methods,
never only by its own derivation. Keep XROTOR's conventions and quantities where they meet
(e.g. `InducedField` reproduces VPUT's disk averages exactly).

* `field.InducedField`: circumferentially averaged induced velocity of XROTOR's rigid
  helical wake anywhere in space. Swirl is exact by Stokes' theorem (zero ahead of the
  disk). Axial and radial velocities come from semi-infinite vortex cylinders (Branlard &
  Gaunaa 2015). Elliptic integrals use Carlson's forms in NumPy; the closed forms of R_C
  lose precision as x -> y, which broke R_J until R_C used duplication. No wake
  contraction yet.
* Measured with it on the Rotol blade: 0.2 R ahead of a disk the axial induced velocity is
  0.59-0.84 of the disk value (least at mid-span), well below the on-axis actuator-disk
  estimate. 0.2 R behind, only 0.58-0.71 of the far-wake value has developed, so VGET with
  weight 1 overstates a close rear rotor's inflow.

* `pair.RotorPair`: each rotor is a full XROTOR solve. It sees the other's `InducedField`,
  imposed through UADD/VADD with the nodes *at its own station radii*: the front wake makes
  the axial velocity jump at each trailing-vortex radius, which a spline between other
  nodes would ring across. XROTOR's free-wake condition (equivalent prop: SETXW, TPQ's
  TWAK/PWAK) uses WA/WT, so each wake pitch responds to the other rotor. Iterate with
  under-relaxation of the imposed inflow (about 0.35x per iteration at w = 0.7). Limits in
  `tests/model/test_pair.py`:
  * At 50 R spacing the front is the isolated rotor, and the rear is XROTOR's
    VPUT/VGET(1, -1).
  * A near-unloaded rear leaves the front isolated.
  * Coincident co-rotating rows are 2-3% below a rotor with twice the blades (each row's
    tip loss is the 2-blade one).
  * Equal-torque contra-rotation leaves no net far-wake angular momentum.
* `RotorPair(interaction="weights")` is a port, not a new model: CROTOR's CR iteration
  (crotor.f v0.85, CROTOR 755-ES1.3, GPL v2+, github.com/mparadkar225/crotor; its XROTOR is
  7.55 with the solver unchanged). The oracle is the shim's transcription (REFCRLOAD,
  REFCRVCLR, REFCRRPM; CROTOR's own COMMON is /CRSHIM/), `tests/differential/test_crotor.py`.
  CROTOR shares one COMMON between the rotors, and two leaks matter bitwise: LOADCR copies
  only NADD profile entries (the rest stay from the rotor solved before, `_shared_tail`),
  and PUTVEL sets NADD = II before its UVADD calls (WA/WT cancel, so rounding only). Its
  default weights are the field model's zero-spacing limit, to first order in the spacing
  (`test_zero_spacing_is_crotors_weights`). Sample rotors: `tests/data/crotor/`.
* NACA Report 747 (`validation/naca747`, `tests/measurement/test_naca747.py`): the tested
  "single" 4-blade propeller is two 2-blade hubs co-rotating 0.16 R apart, so model it as a
  co-rotating RotorPair. Dual-minus-single peak efficiency matches (+0.022/+0.029 at 45/50°,
  model +0.023/+0.029); the model's dual is ~0.6° of pitch more loaded than the test while
  single matches, and figure 12 (front/rear CP) puts that mostly on the front rotor
  (+3-23% vs rear -1 to +13%). Slipstream contraction is negligible there: an
  actuator-disk estimate (`validation/naca747/contraction_estimate.py`) shrinks the
  slipstream 0.1-1.5% at 0.16 R and changes the rear's power by < 0.15% (< 0.8% even fully
  developed). Digitizing: a crop whose lowest grid line is not the axis shifts every
  value by one step; check CT·J/CP against the report's efficiency curves. Solid vs
  dashed lines: separate by connected-component length (`validation/naca747/digitizing`).
* `cylinder_velocity` on a sheet: form 1 - k0^2 as ((a-r)/(a+r))^2 (`_scaled_pi`).
  Subtraction underflows within ~1e-8 of r = a and gave NaN.

## Licensing

GPL v3 or later (XROTOR is GPL v2 or later; xrotor-python is GPL v3).
