# pyxrotor — notes for agents working in this repo

Python port of XROTOR 7.55. The compiled 7.55 Fortran in `reference/` is the oracle; the
port is only as good as its differential tests against it.

## Commands

```bash
git submodule update --init      # reference sources
reference/build.sh               # build both reference libraries (~1 min)
pytest                           # everything (~1.5 min)
pytest tests/differential -q     # just the port-vs-Fortran tests
HYPOTHESIS_PROFILE=deep pytest tests/differential -W error::RuntimeWarning   # 20k cases/test
python tools/common_writes.py TPQ APITER    # COMMON variables a Fortran routine assigns
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
   * MIN/MAX with NaN: in the compiled reference every site tested returns the *first*
     argument (`fmath.fmin/fmax`, same as Python's built-ins). gfortran documents this as
     unspecified and a standalone probe behaved differently, so pin new sites with a
     non-finite test against the compiled routine (see `test_clcdcm_non_finite`).
   * Convert inputs to `float` at public entry points: NumPy scalars arriving from arrays
     (e.g. `xi[i]`) silently change overflow and divide-by-zero behavior.
   * Indices are 0-based in Python. Fortran IAERO = N is `iaero = N - 1`; "unset" (0) is -1.
     `State.INDEX_VARS` lists such variables; statesync converts them.
6. **Known 7.55 bugs are reproduced first, then fixed.** Port the faithful behavior (mark it
   `# 7.55:` in the code), prove it matches, then add a fix behind an option, with a test for each.
7. **Keep Fortran names in comments** next to their Python equivalents and cite the
   source routine in each docstring. Code with no effect on state (unused locals, prints,
   plotting) may be omitted; say so in the docstring.

## Reference harness facts

* `Ref755()` loads a private copy of the library, so instances never share state.
* `ref.call(...)` passes everything by reference and returns the converted arguments, so
  output arguments can be read back. Python bool → LOGICAL, int → INTEGER*4, float → REAL*8.
  COMMON arrays can be passed as arguments directly: `ref.call("HELICO", ..., ref["GAM"], ...)`.
* A Fortran `STOP` kills the Python process: run anything that might STOP in a subprocess.
* Interactive routines (ASKR prompts) can't be called; transcribe them into
  `reference/shim/ref755_shim.f` statement for statement, with prompts turned into arguments,
  and keep fixed-form lines within 72 columns. OPER commands are there as REFVELO, REFRPM, ...
* `reference/build.sh` must be rerun after editing the shim.

## Known issues in XROTOR 7.55 (all reproduced faithfully unless noted)

* **Stale aero-section indices after DESI.** If aero sections are defined before the first
  design, `IAERO` is computed against the pre-design radial grid and DESI never refreshes
  it, so outer stations extrapolate past the last section (negative Cd). In a P-51D test this
  overstated design efficiency (85.3% vs 79.8%). ARBI and LOAD call SETIAERO and are fine.
  Fix: call SETIAERO after DESI. (DESI not yet ported.)
* **GETCLCDCM's fallback picks the wrong section.** When a station's IAERO is invalid (as
  INIT leaves it, 0 everywhere, until SETIAERO runs), the recovery search exits with N at the
  first section *outboard* of the station instead of the inner bounding one.
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
* **DESI fails** (GETALF inversion, NaN wake) when the helical tip Mach number reaches about 1.
* The generic compressibility model (dCD = K (M − Mcrit)^3, local Mach capped at 0.99) is
  very sensitive to Mcrit for transonic tips.

## Licensing

GPL v3 or later (XROTOR is GPL v2 or later; xrotor-python is GPL v3).
