# pyxrotor

A Python + NumPy port of [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/) 7.55, Mark
Drela and Harold Youngren's propeller and windmill design and analysis code, written to be
used as a library and verified routine by routine against the original Fortran.

## Status

| Area | Status |
|---|---|
| Reference harness (XROTOR 7.55 as a callable library, full COMMON-block access) | Done |
| `SETX` radial grid | Ported, bit-identical to 7.55 |
| Parametric section aero (`CLCDCM`, `GETCLCDCM`, `GETALF`) | Next |
| OPER (operating-point solver, wake models, sweeps) | To do |
| DESI (minimum-induced-loss design), MODI | To do |
| BEND, NOIS, file formats | To do |

## Quick start

```bash
git clone --recursive https://github.com/mparadkar225/pyxrotor
cd pyxrotor
reference/build.sh            # needs gfortran, gcc, cmake, libX11 headers
pip install -e ".[test]"
pytest
```

## How it is checked

* **Against the Fortran** (`tests/differential/`): each ported routine is run beside the
  compiled 7.55 routine on randomized inputs (Hypothesis), comparing outputs and COMMON
  state. Straight-line routines must match exactly; iterative solvers to about 1e-12.
* **Against measurement** (`tests/measurement/`, `validation/`): wind-tunnel data, starting
  with ARC R&M 2357 (Rotol wooden Spitfire IX propeller), at engineering tolerances.
* **The references themselves** (`tests/reference/`): the library build reproduces the
  interactive program's results.

See [`reference/README.md`](reference/README.md) for how the references are built.

## License

GPL v3 or later. XROTOR is GPL (v2 or later) and pyxrotor is a derivative work.
