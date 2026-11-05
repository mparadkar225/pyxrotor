# Reference implementations

pyxrotor is checked routine by routine against compiled Fortran. Both sources are git
submodules pinned to exact commits; `build.sh` compiles them into `reference/build/`.

| Reference | Source | Role |
|---|---|---|
| XROTOR 7.55 | `xrotor-7.55/` (MIT source, mirrored at mparadkar225/xrotor) | **The oracle.** Every ported routine is tested against it. |
| xrotor-python | `xrotor-python/` (DARcorporation, 2019, GPL v3) | Structural reference only (see below). |
| CR-176883 data banks | `cr176883/` (transcribed from NASA CR-176883, 1986) | Oracle for `pyxrotor.databank`. |

## XROTOR 7.55 (`libxrotor755.so`, `ref755.py`)

The whole program is built as a shared library. The main program becomes a subroutine,
and `shim/ref755_shim.f` adds non-interactive versions of the start-up code and ARBI,
transcribed statement for statement with prompts replaced by arguments.

`commons.py` parses `XROTOR.INC` into the exact memory layout of all 27 COMMON blocks
(a test checks every block size against the compiled symbol), so Python can read and
write any XROTOR state variable by name:

```python
from reference.ref755 import Ref755
ref = Ref755()                         # private copy of the library: isolated state
ref["II"], ref["XI0"] = 40, 0.2
ref.call("SETX")                       # any subroutine; arguments by reference
xi = ref.get("XI")[:40]
_, vso, rho, rmu = ref.call("ATMO", 7.62, 0.0, 0.0, 0.0)   # outputs come back
```

Limitations: a Fortran `STOP` ends the Python process, and CHARACTER arguments are not
supported in `call` (set CHARACTER COMMON variables directly instead).

## xrotor-python (`build/xrotor-python/`, `refxp.py`)

A 2019 conversion of XROTOR to Fortran 90 modules, with DESI, MODI, plotting, and the
parametric section model removed and replaced by linear interpolation of tabulated polars.
The build switches its C interface from float32 to float64.

It currently fails its own bundled test case (Newton's iteration diverges at every
optimization level, with or without the patch), and its section model is not 7.55's, so
nothing is validated against it. It is kept for its design: all COMMON state moved into
one explicit context object, and a small library API (`set_case`, `operate`, ...), which
is a useful pattern for pyxrotor's own structure.

## CR-176883 airfoil data banks (`libcr176883.so`, `crdb.py`)

The Clark Y, NACA 4-digit and NACA 16-series data banks of NASA CR-176883, transcribed
from the report's listings (see `cr176883/README.md`). Built with the same settings as
XROTOR, so `pyxrotor.databank` is compared with them bit for bit:

```python
from reference.crdb import CRDB
out = CRDB().call("CLCD16", 0.5, 4.0, 0.0, 0.0, 1, 1, 1.0, 0.06, 0.0, 0.3, 0.75,
                  wa, vis, den, 0.0, 0.0, 0.0, 0.0)    # arguments come back, as in ref.call
```

## Build settings

`-O1 -fdefault-real-8 -fdefault-double-8 -ffp-contract=off`: all reals are IEEE double
(matching Python floats), and fused multiply-add is disabled so results are the same on
x86-64 and ARM. `BUILDINFO.json` records the compiler and submodule commits used.

Note that the stand-alone interactive build used during exploration (`-fdefault-real-8`
only) keeps DOUBLE PRECISION at 16 bytes, which touches a few `DBLE`/`DEXP` calls in the
stall model; its printed results agree with this build to print precision, not bitwise.
