# pyxrotor

A Python + NumPy port of [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/) 7.55, Mark
Drela and Harold Youngren's propeller and windmill design and analysis code, written to be
used as a library and verified routine by routine against the original Fortran.

## Status

| Area | Status |
|---|---|
| Reference harness (XROTOR 7.55 as a callable library, full COMMON-block access) | Done |
| `State`: all 27 COMMON blocks, generated from XROTOR.INC | Done |
| Start-up (`INIT`, `SETDEF`), `ATMO`, `ARBI` blade input | Ported, bit-identical |
| Section aero (`CLCDCM`, `GETCLCDCM`, `GETALF`, `SETIAERO`, `SORTAR`) | Ported, bit-identical |
| Induced velocity: graded momentum, potential (Goldstein), discrete vortex wake | Ported, bit-identical |
| OPER solver (`TPQ`, `APINIT`, `APITER`, `APER`) and wake geometry | Ported, bit-identical |
| OPER commands: `VELO`, `ANGL`, `ADVA`, `RPM`, `THRU`, `TORQ`, `POWE` (fix pitch or rpm), engine power line | Ported, bit-identical |
| Sweeps (`ASEQ`, `RSEQ`, `BSEQ`, `VSEQ`: `SETCAS` plus `CASPLT`'s solves), case table (`ADDC`, `CLRC`, `CASE`, `CPUT`, `CGET`) | Ported, bit-identical; files byte-identical |
| Slipstream profiles `VPUT`/`VGET` (`SAVVEL`, `GETVEL`), for coupling rotors |  Ported, bit-identical; files byte-identical |
| Nacelle/spinner body (`NACE`: `NACELL`, `SRCVEL`), OPER `ADW` factor | Ported, bit-identical (source-position bug reproduced; fix optional) |
| `SAVE`/`LOAD` restart files (7.55 format) | Ported, bit-identical; files byte-identical (old pre-6.9 format not ported) |
| DESI design (`SETCLD`, `DESINI`, `DEITER`, `DESGEN`) and commands `INPU`, `EDIT`, `FORM`, `ITER`, `N` | Ported, bit-identical (stale-section bug reproduced; fix optional) |
| MODI | To do |
| **Beyond XROTOR:** induced-velocity field off the disk (`Rotor.induced_field()`), the base for rotor-rotor interaction | New model, validated against exact limits and an independent discrete vortex wake (`tests/model`) |
| BEND, NOIS, xrotor.def, PVAR files | To do |

```python
from pyxrotor import Rotor, AeroSection

r = Rotor()                                   # XROTOR start-up state
r.atmosphere(0.0)                             # km
r.set_sections([AeroSection(xisect=0.0, cdmin=0.010), AeroSection(xisect=0.7, cdmin=0.008)])
r.arbi(nblds=2, vel=50.0, tip_radius=1.64, hub_radius=0.33,
       r_over_R=[0.25, 0.5, 0.75, 1.0], c_over_R=[0.14, 0.16, 0.14, 0.05],
       beta_deg=[45.0, 30.0, 22.0, 17.0])
r.velocity(40.0)
r.solve_rpm(1200.0)                           # or solve_thrust(800.0, fix="rpm"), ...
print(r.performance())                        # thrust, power, efficiency, CT, CP, ...

d = Rotor(fix_all=True)                       # 7.55 bug fixes on (default: 7.55 exactly)
d.atmosphere(7.62)
d.set_sections([...])
d.design(nblds=4, tip_radius=1.702, hub_radius=0.20, vel=112.0, rpm=1437.0,
         power=1118550.0, cl=0.5)             # DESI INPU: minimum-induced-loss blade
d.redesign(nblds=3, cl_root_tip=(0.6, 0.4))   # DESI EDIT

d.sweep_velocity(60.0, 140.0, 10.0, constant_speed=True)   # VSEQ: list of case dicts
d.sweep_advance_ratio(0.3, 1.2, 0.1)          # ASEQ (appends, like XROTOR; clear_cases())
d.write_cases("p51.cas")                      # CPUT
d.save("p51.xr")                              # SAVE: readable by XROTOR 7.55
r2 = Rotor(); r2.load("p51.xr")               # LOAD: XROTOR's own .xr files too
r2.nacelle("spinner.nac")                     # NACE: body area file -> UBODY at the blades

front.solve_rpm(1500.0)                       # XROTOR's way to put a rotor behind another:
p = front.slipstream_profile()                # VPUT (or front.save_slipstream("f.vel"))
rear.set_slipstream(p["r"], p["u"], p["v"],   # VGET: swirl weight -1 = counter-rotating
                    axial_weight=1.0, swirl_weight=-1.0)
r2.nacelle(z=[0, 0.5, 1.0], area=[-0.1, -0.35, -0.45], rotor_z=0.15, wake_area=-0.3)
```

## 7.55 bugs and fixes

By default pyxrotor reproduces XROTOR 7.55 bit for bit, bugs included (they're listed in
[`CLAUDE.md`](CLAUDE.md)). The ones that change results have a fix you can turn on:

| Option | Fixes |
|---|---|
| `Rotor(refresh_design_sections=True)` | a design keeps the previous grid's aero-section indices (outer stations extrapolate past the last section) |
| `Rotor(fix_linear_cl=True)` | a linear design CL is set on the previous grid (a first design gets the root CL everywhere) |
| `Rotor(fix_section_fallback=True)` | a station with no valid section index uses the wrong section (past the last section: the root one) |
| `Rotor(fix_slipstream_splines=True)` | LOAD leaves the slipstream profile's spline slopes stale (zero) |
| `Rotor(fix_all=True)` | all four |
| `nacelle(..., fix_source_positions=True)` | NACE's body sources are shifted one interval downstream |
| `load(..., keep_free_wake=True)` | LOAD sets the free-wake flag from SAVE's duct flag |
| `read_cases(..., full_lines=True)` | CGET drops power, thrust, torque and efficiency |

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
  state. Everything ported so far matches bit for bit, including the Newton solvers.
* **Against measurement** (`tests/measurement/`, `validation/`): wind-tunnel data, starting
  with ARC R&M 2357 (Rotol wooden Spitfire IX propeller), at engineering tolerances.
* **The references themselves** (`tests/reference/`): the library build reproduces the
  interactive program's results.
* **New models** (`tests/model/`), which have no Fortran to match: exact limits and
  independent methods (closed forms against brute-force Biot-Savart, the averaged field
  against a discrete helical vortex wake).

See [`reference/README.md`](reference/README.md) for how the references are built.

## License

GPL v3 or later. XROTOR is GPL (v2 or later) and pyxrotor is a derivative work.
