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
| **Beyond XROTOR:** contra-rotating and tandem pairs (`RotorPair`): coupled solve; fixed pitch, power or thrust with a torque split, power at fixed pitch; gear ratio; shared nacelle | New model, validated by limits (`tests/model/test_pair.py`) and NACA Report 747 tunnel data |
| **Beyond XROTOR:** section estimates (`pyxrotor.sections`: Clark Y, NACA 16-series) with a critical Mach number calibrated on transonic propeller tests | New model, checked against six flight and tunnel datasets (`tests/measurement/test_calibrated_sections.py`) |
| **Beyond XROTOR:** round-shank drag (`pyxrotor.shank`, TR-852's measured cylinder drag on a no-lift placeholder section) | New model, exact against XROTOR's viscous totals (`tests/model/test_shank.py`); checked on the YP-39 and RM L6L27a (`tests/measurement/test_shank_drag.py`) |
| **Beyond XROTOR:** NASA CR-176883's airfoil data banks (`pyxrotor.databank`: Clark Y, NACA 4-digit, NACA 16-series CL and CD from thickness, camber, angle and Mach), and `fit_section` / `clark_y_section` / `naca16_section` to turn one into XROTOR sections (the transonic comparisons rerun with them still need the calibrated shift: `validation/cr176883_sections`) | Port of the report's Fortran, transcribed in `reference/cr176883`: bit for bit against it (`tests/differential/test_databank.py`), and the transcription reproduces the report's sample output (`tests/model/test_databank_sample.py`) |
| **Beyond XROTOR:** swept (scimitar) blades in the compiled core (`Rotor.set_sweep`, `pyxrotor.sweep`): the quarter-chord line's geometry, the effective sweep of the local flow at every iteration, sections in the normal plane, and the curved blade's own induction (swept discrete vortex wake and curved bound vortex) | Exact limits (zero sweep, rigid rotation and axial shift give XROTOR's VRTXCO), an independent 3-D geometry check, convergence with stations (`tests/model/test_sweep.py`); SR-3 propfan (`tests/measurement/test_propfan.py`): efficiency to Mach 0.85, but the blade angle is 5 deg low and the fixed-pitch power trend wrong in sign |
| **Beyond XROTOR:** transonic lift loss past a lift-divergence Mach number (`Rotor.set_transonic_lift`, `pyxrotor.transonic`), fitted to the Clark Y data bank's lift | Fit and exact limits (`tests/model/test_transonic.py`); P-47 flight tests (`tests/measurement/test_transonic_lift.py`, `validation/transonic_lift`) |
| **Beyond XROTOR:** transonic solves (`pyxrotor.robust`: power at fixed rpm, fixed pitch) that keep working near the local-Mach cap | New, checked against XROTOR's own converged solutions (`tests/model/test_robust.py`) |
| CROTOR's counter-rotation iteration (`RotorPair(interaction="weights")`: VPUT/VGET with fixed weights, CROTOR 755-ES1.3 `CROT` `RPM` and `POWE` with `FIX`/`SYNC`) | Ported, bit-identical to a shim transcription of crotor.f; it is the field model's zero-spacing limit |
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

### Transonic tips

```python
from pyxrotor import robust, sections

secs = [sections.naca16(t, x, 0.5) for x, t in [(0.2, 0.30), (0.5, 0.09), (1.0, 0.05)]]  # r/R, t/c
r.set_sections(secs)
sol = robust.solve_power(r, 1.5e6, rpm=1300.0)   # constant speed: blade angle for 1.5 MW
sol.usable, sol.exact, sol.method, sol.efficiency, sol.thrust, sol.beta75_deg
```

XROTOR's compressibility model with the usual thickness estimate Mcrit = 0.85 - 1.2 t loses
far too much efficiency at transonic tips (0.1-0.9 at helical tip Mach 0.95-1.17, on every
propeller tested). `pyxrotor.sections` adds 0.15 to that estimate at 5% thickness and above
and 0.20 at 4.5% and below, which fits P-47 and YP-39 flight tests, YP-47M and SR-2 tunnel
tests and RM L51L28 static tests (`validation/`); `mcrit_shift=0.0` gives the plain estimate.
Near the 0.99 local-Mach cap XROTOR's Newton iteration often fails or ends in a small limit
cycle, its power solve more often than its fixed-pitch one. `robust.solve_power` tries POWE,
then POWE continued in Mach number from the solution without compressibility, then
interpolation between fixed-pitch solves (accepting small limit cycles);
`robust.solve_fixed_pitch` is RPM with the limit-cycle acceptance. `sol.exact` says whether
the rotor holds a converged XROTOR solution of the request. Swept (scimitar) blades are not
predicted: XROTOR's straight lifting line needs 6-8 deg less blade angle than the SR-3
propfan (`validation/nasa_propfan`).

### Contra-rotating propellers

```python
from pyxrotor.pair import RotorPair

pair = RotorPair(front, rear, spacing=0.25)   # two set-up Rotors, rear 0.25 m behind
res = pair.solve_power(1.5e6, rpm=1300.0, torque_split=0.5)   # equal torques
res.thrust, res.efficiency, res.front["torque"], res.rear["torque"], res.net_torque
res.swirl_power_fraction                      # far-wake swirl loss (~0 when balanced)
pair.field(r, z)                              # the pair's induced velocity anywhere
```

Each rotor is a full XROTOR solve that sees the other's averaged induced velocity at the
actual spacing (axial only ahead of a disk; the slipstream's axial velocity and swirl
behind it). The two are iterated to convergence. Not modeled yet: slipstream contraction,
viscous wakes, blade-passing unsteadiness.

`RotorPair(front, rear, 0.0, interaction="weights")` couples the rotors the way XROTOR
itself can, as Philip Carter's CROTOR automates it: each rotor's VPUT slipstream is imposed
on the other with fixed weights (default front 0.5 axial, no swirl; rear 0.5 axial, swirl
-1), alternating solves until each thrust settles to 0.1%. It reproduces CROTOR bit for bit
and is what the field model tends to as the spacing goes to zero (error about 1e-4 per
1e-4 R); at 0.1-0.2 R the field model puts 5-9% more of the thrust on the front rotor.

`pair.solve_power_split(power, rpm, power_split=0.5, fix="rpm"|"pitch", locked=True)` is
CROTOR's POWE: the power split between the rotors, at fixed rpm (blade angles set), at fixed
pitch with each rotor at its own rpm (`locked=False`), or at fixed pitch geared together at
`gear_ratio` (`locked=True`, CROTOR's SYNC; weights interaction only, and the total power
is only near the target: use `solve_power_fixed_pitch` for the exact one).

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

`pip install` compiles the operating-point solver loop with Cython, so it needs a C
compiler. A solve takes about 5.5 ms (the Fortran: 1.8 ms) for the Rotol propeller.

## How it is checked

* **Against the Fortran** (`tests/differential/`): each ported routine is run beside the
  compiled 7.55 routine on randomized inputs (Hypothesis), comparing outputs and COMMON
  state. Everything ported so far matches bit for bit, including the Newton solvers.
* **Against measurement** (`tests/measurement/`, `validation/`): wind-tunnel data at
  engineering tolerances: ARC R&M 2357 (Rotol wooden Spitfire IX propeller) and NACA
  Report 747 (four-blade single- and dual-rotating propellers, for `RotorPair`). On 747 the
  counter-rotation peak-efficiency gain matches the test to 0.001 at 45-50° and is about
  0.01 optimistic at low pitch; see [`validation/naca747`](validation/naca747/README.md).
  Compressibility: P-47 and YP-39 flight tests, the YP-47M altitude-wind-tunnel series, RM
  L51L28 static tests and NASA's SR-2/SR-3 propfans (`validation/naca_p47`, `naca_yp39`,
  `naca_yp47m`, `naca_l51l28`, `nasa_propfan`).
* **The references themselves** (`tests/reference/`): the library build reproduces the
  interactive program's results.
* **New models** (`tests/model/`), which have no Fortran to match: exact limits and
  independent methods (closed forms against brute-force Biot-Savart, the averaged field
  against a discrete helical vortex wake).

See [`reference/README.md`](reference/README.md) for how the references are built.

## License

GPL v3 or later. XROTOR is GPL (v2 or later) and pyxrotor is a derivative work.
