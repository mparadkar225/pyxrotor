# Propeller sizing for a 500 mph, 2900 hp fighter

A design exploration with pyxrotor (not a validation). `study.py` designs a
minimum-induced-loss blade for each layout at the maximum-speed point, over diameters of
11.2–14.1 ft, helical tip Mach numbers of 0.95–1.15 at 500 mph and design CL 0.4 and 0.5. It
then checks each design off-design. All results are in `results.json`.

## Requirement and assumptions
| | |
|---|---|
| maximum speed | 500 mph at 30,000 ft with 2900 hp (WEP held to critical altitude); also 25,000 and 35,000 ft at the same power |
| high-speed cruise | 350 mph at 25,000 ft, 1300 shp, 80% rpm (constant speed) |
| take-off | sea level, 2900 hp, maximum rpm, judged at 130 mph (XROTOR's static power is unreliable with stalled roots) |
| sections | NACA 16-series, t/c 0.30 at 0.2R (faired, no round shank) to 0.055 at the tip; calibrated Mcrit |
| layout | spinner 0.22R; contra rows 0.12 D apart, equal torques, same rpm |

The cruise power comes from scaling the max-speed thrust power by ρV³ for parasite drag,
with about 5% induced drag.

## Candidates (design CL 0.5; CL 0.4 is 0.01–0.02 worse everywhere)
| layout | D | prop rpm (gear at 2800) | Mt at 500 mph | AF per blade | β0.75 design | η 500 mph at 25k / 30k / 35k ft | η cruise | η climb 130 mph |
|---|---|---|---|---|---|---|---|---|
| 4-blade | 13.1 ft | 1083 (0.387) | 1.05 | 126 | 53° | 0.827 / 0.829 / 0.785 | 0.870 | 0.703 |
| 4-blade | 14.1 ft | 1007 (0.360) | 1.05 | 109 | 53° | 0.837 / 0.843 / 0.800 | 0.881 | 0.725 |
| 5-blade | 12.1 ft | 1171 (0.418) | 1.05 | 120 | 53° | 0.824 / 0.827 / 0.781 | 0.865 | 0.687 |
| 5-blade | 13.1 ft | 1083 (0.387) | 1.05 | 103 | 53° | 0.835 / 0.841 / 0.797 | 0.878 | 0.711 |
| 6 contra (3+3) | 12.1 ft | 938 (0.335) | 0.95 | 131 / 114 | 62° / 60° | 0.878 / 0.872 / 0.846 | 0.898 | 0.721 |
| 6 contra (3+3) | 13.1 ft | 979 (0.350) | 1.00 | 96 / 84 | 58° / 57° | 0.884 / 0.870 / 0.834 | 0.912 | 0.749 |
| 8 contra (4+4) | 11.2 ft | 1021 (0.365) | 0.95 | 119 / 102 | 62° / 60° | 0.879 / 0.874 / 0.848 | 0.895 | 0.705 |
| 8 contra (4+4) | 12.1 ft | 938 (0.335) | 0.95 | 101 / 87 | 61° / 60° | 0.886 / 0.883 / 0.860 | 0.902 | 0.728 |

Contra AF and β are front / rear.

## Findings
- **Single rotation peaks at tip Mach 1.05–1.10 at 500 mph.** Lower rpm loads the blade so
  heavily that swirl losses dominate: at tip Mach 0.95 the induced efficiency (including
  swirl) is only 0.85. Higher rpm loses more to compressibility.
- **Contra-rotation recovers the swirl.** It gains 0.03–0.05 at maximum speed and 0.02–0.04
  at cruise, and its best rpm is lower (tip Mach 0.95–1.0): with no swirl to recover, only
  compressibility is left.
- **A 5th blade buys 0.01 over 4** at the same diameter, or lets the diameter shrink by
  about 1 ft for the same efficiency.
- **8 contra vs 6 contra.** At the same diameter the 8-blade gains 0.005–0.01. Its use is a
  smaller diameter (11–12 ft) at the same efficiency, with narrower blades.
- **Altitude.** At 35,000 ft the same 2900 hp loads the blades more, and every layout loses
  0.03–0.05 (the contras the least).

## Limits
- Absolute efficiencies are probably 0.03–0.06 optimistic (the calibrated model's low-speed
  bias against the wartime tests); compare layouts with each other.
- No shank or spinner drag, blade weight, noise, or installation effects.
- The contra rows were designed without compressibility, each in the other's field. The
  single-rotation blades were designed with it, which is worth about 0.01 at tip Mach 1.05.
- Some contra solves near the local-Mach cap fail; they are recorded as not ok.
