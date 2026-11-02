# pyxrotor — notes for agents working in this repo

Python port of XROTOR 7.55. The compiled 7.55 Fortran in `reference/` is the oracle; the
port is only as good as its differential tests against it.

## Commands

```bash
git submodule update --init      # reference sources
reference/build.sh               # build both reference libraries (~1 min)
pytest                           # everything
pytest tests/differential -q     # just the port-vs-Fortran tests
HYPOTHESIS_PROFILE=deep pytest tests/differential -W error::RuntimeWarning   # 20k cases/test, ~6 min
```

Run the deep profile after porting or changing a routine; the default 400 cases have
missed real mismatches (an overflow in a diverging GETALF only showed up on other seeds).

## Porting rules

1. **XROTOR 7.55 is the source of truth.** Port from `reference/xrotor-7.55/src/*.f`.
   xrotor-python is not an oracle (it fails its own test case and uses a different
   section model); use it only as a structural example.
2. **One routine at a time, each with a differential test** in `tests/differential/`,
   using Hypothesis to drive inputs through both the port and `Ref755`. Set inputs in
   COMMON with `ref["NAME"] = ...`, call the Fortran with `ref.call("NAME", ...)`, and
   compare outputs and any COMMON state the routine writes.
3. **Tolerances**: `rtol=0` for straight-line arithmetic (SETX matches bit for bit);
   `harness.RTOL` (1e-12) for iterative solvers. Anything looser needs a written reason.
4. **Match the arithmetic, not just the math.** Things that have already mattered:
   * `PI` is XROTOR's `3.141592654` (`pyxrotor.constants.PI`), never `math.pi`.
   * Use `math.sin/cos/sqrt/exp` for scalars inside loops: same libm as the Fortran.
     NumPy's vectorized transcendentals can differ by an ulp.
   * Write `x**2` from Fortran as `x * x` (gfortran's form).
   * Keep loop order and accumulation order in iterative solvers; vectorize only where a
     differential test proves it changes nothing.
   * Fortran evaluates `a*b*c` and `a/b*c` left to right; keep the same grouping.
     `x**3` with an integer exponent is `(x*x)*x`; with a REAL exponent it is libm `pow`.
   * Python raises where Fortran quietly produces inf/NaN: use `aero._div` for divisions
     that can meet zero and `aero._pow` for powers that can overflow. Differential tests
     deliberately drive iterations into divergence, so these paths do get exercised.
   * Convert inputs to `float` at public entry points: NumPy scalars arriving from arrays
     (e.g. `xi[i]`) silently change overflow and divide-by-zero behavior.
   * Indices are 0-based in Python. Fortran IAERO = N is `iaero = N - 1`; "unset" (0) is -1.
5. **Known 7.55 bugs are reproduced first, then fixed.** Port the faithful behavior behind
   a compatibility flag, prove it matches, then make the fix the default, with a test for each.
6. **Keep Fortran names in comments** next to their Python equivalents
   (`dxi  # DXI`), and cite the source routine in each docstring.

## Reference harness facts

* `Ref755()` loads a private copy of the library, so instances never share state.
* `ref.call(...)` passes everything by reference and returns the converted arguments, so
  output arguments can be read back. Python bool → LOGICAL, int → INTEGER*4, float → REAL*8.
* A Fortran `STOP` kills the Python process: run anything that might STOP in a subprocess.
* Interactive routines (ASKR prompts) can't be called; transcribe them into
  `reference/shim/ref755_shim.f` statement for statement, with prompts turned into arguments,
  and keep fixed-form lines within 72 columns.

## Known issues in XROTOR 7.55

* **Stale aero-section indices after DESI.** If aero sections are defined before the first
  design, `IAERO` is computed against the pre-design radial grid and DESI never refreshes
  it, so outer stations extrapolate past the last section (negative Cd). In a P-51D test this
  overstated design efficiency (85.3% vs 79.8%). ARBI and LOAD call SETIAERO and are fine.
  Fix: call SETIAERO after DESI.
* **GETCLCDCM's fallback picks the wrong section.** When a station's IAERO is invalid (as
  INIT leaves it, 0 everywhere, until SETIAERO runs), the recovery search exits with N at the
  first section *outboard* of the station instead of the inner bounding one. Ported
  faithfully and tested, including what it writes back to IAERO.
* **GETALF starts every station from section 1's zero-lift angle**, not the station's own.
* **DESI fails** (GETALF inversion, NaN wake) when the helical tip Mach number reaches about 1.
* The generic compressibility model (dCD = K (M − Mcrit)^3, local Mach capped at 0.99) is
  very sensitive to Mcrit for transonic tips.

## Licensing

GPL v3 or later (XROTOR is GPL v2 or later; xrotor-python is GPL v3).
