# NASA CR-176883 airfoil data banks

Korkan, K. D., Camba, J., III, and Morris, P. M., "Aerodynamic data banks for Clark-Y,
NACA 4-digit and NACA 16-series airfoil families", Texas A&M University, NASA CR-176883,
1986 ([NTRS 19860021221](https://ntrs.nasa.gov/citations/19860021221)).

The report is three Fortran programs that give section CL, CD (and CM) from thickness,
camber, angle of attack and Mach number. Their data are tables in DATA statements and
curve fits:

| bank | source data | file | report pages (PDF) |
|---|---|---|---|
| Clark Y | Haines & Monaghan, R&M 2036 (1945): t/c 0.04-0.35, Mach to ~0.9 | `clarky.f` | 40-62 (48-69) |
| NACA 4-digit, camber at 0.5c | McCormick et al., "A study of torpedo propellers", 1956 | `naca4.f` | 102-111 (109-118) |
| NACA 16-series | curve fits (t/c 0.03-0.4, Mach to supersonic) | `naca16.f` | 141-158 (148-165) |
| standard atmosphere | ATMCON, shared | `atmcon.f` | |

## Transcription
Every line was read from the 300 dpi page images (OCR was too noisy: minus signs, short
numbers). Kept statement for statement, with the original constants and operation order.
Omitted: the limitation warnings (WRITE/FORMAT only), the main programs' output, and the
moment and coordinate routines (CMFF, CYCOOR, QCCM, FDCOOR, FDCR, AREA, the 16-series
CMnFF), which don't affect CL or CD. The main programs' CL/CD arithmetic is described in
each file's header (the 4-digit one is the subroutine `NACA4`).

Renamed so that all files link into one library: the 16-series CLCD is `CLCD16`, the Clark Y
CLCD is `CLCDCY`, the 4-digit ATMCON (other arguments) is `ATMCN4`; the Clark Y ATMCON is
the same as `atmcon.f` and is dropped.

## Checks against the report
`sample_*.f` reproduce the report's sample cases from the transcription compiled as written
(default REAL*4); `sample_*.txt` hold the printed values, read from the output pages.
`tests/model/test_databank_sample.py` compiles and runs them:

| bank | sample case | printed rows | largest difference |
|---|---|---|---|
| Clark Y | t/c 0.117 at Mach 0-0.4 and t/c 0.30 at Mach 0.6-0.85, 10,000 ft, chord 0.75 ft | 42 | 2e-6 in CL, 1e-6 in CD |
| 4-digit | CLD 0.2, t/c 0.12, Mach 0.5, 10,000 ft, alpha -10 to 8 deg | 10 | 1e-6 (CL), CD exact |
| 16-series | t/c 0.04, CLD 0.22815 (ALPHALO -2 deg), Mach 0.45 and 0.70 | 6 | 1e-6 (CD) |

The last-digit differences are single-precision rounding on the 1986 machine. The 16-series
sample output is not the listing's DATA case: that version of the program took ALPHALO as
input (CLD = 0.1096 sqrt((1+t)/(1-t)) x 2 = 0.228150).

The samples don't reach every table (Clark Y: CLSF, CDSF and most CDnnF Mach tables; 4-digit:
camber above 0.04; 16-series: supersonic and t/c > 0.25). Those were checked by reading only.

## What the banks say about pyxrotor.sections
`databank.fit_section` fits XROTOR's section parameters to a bank at Mach 0.3, Re 2e6
(`tests/model/test_databank_sample.py` records these):

| | 16-series (t 4-9%) | Clark Y (t 6-12%) |
|---|---|---|
| Mcrit (CD rise 0.002 at CLDMIN) vs 0.85 - 1.2 t | +0.02 to +0.08 | -0.02 to +0.08 |
| DCDCL2 vs the estimates' 0.007 | ~0.06 (8x) | 0.006-0.008 |
| CDMIN vs the estimates | 15-50% lower | within 0.001, CL at minimum 0.19 (est. 0.35) |

So the calibrated +0.15/+0.20 in `pyxrotor.sections` is not a section property: the 2-D data
put the drag rise near the thickness estimate, and the propeller tests want much more
(tip relief, the cubic drag law's shape, the compressible stall limit).

## Listing oddities (kept as printed)
- **4-digit CDL, X8.** The continuation line reads `01,.2,.3,.4,.5` with the `0` in column 6.
  A zero there is not a continuation mark, so it must be the letter O (a valid one), and
  X8(13) = 1. as the program compiles. The table's steps say .1 was meant; `NACA4(...,
  LFIXX8=.TRUE.)` (`databank.naca4(fix_x8=True)`) uses .1. It changes the drag due to lift
  for camber above 0.04 between X = 0 and 1.
- **4-digit CLDVSZ.** Its two tables (camber from CLD, CLD from camber) are only roughly
  inverse: CLD 0.4 at 12% thickness comes back as 0.445.
- **Clark Y:** see the header of `clarky.f` (HMLF's dead drag-critical table, CDSF's
  interval lookup, AOF's slope table, CD23F's `.2`, CD30F's non-monotonic Mach row).
- **16-series AOAT** (a warning check, omitted): `DIMENSION ALL(5)` with six DATA values.
