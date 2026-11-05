# Digitizing the YP-47M reports

The scripts are as they were run on 300 dpi page images from naca-archive's debug-figs
workflow, from a folder holding `<NTRS id>/pages/p<N>-<NN>.png`. The images aren't kept here.

**Blade form (figure 1 of each report):** `bladeform.py` → `bladeform_traces.json`
- The curves are followed column by column (`dig2.follow`), skipping grid columns and runs
  along grid rows. Stretches lying along a grid line are bridged with `dig2.extend`.
- Calibration comes from the detected grid rows and columns. Each curve's scale is value =
  a + b·k, where k is the row count above the bottom frame line, read from the axis labels
  (`FIGS`).
- The design CL of the Aeroproducts, 838 and 836 blades is read by eye (flat outboard,
  `manual`).
- `make_geometry.py` → `../geometry.json`. It also replaces the Curtiss 836 b/D inboard of
  0.55R with values read by eye (the trace followed the blade-angle curve where they cross),
  and adds rounded tips.
- **Check:** integrated from 0.15R, the traces give activity factors of 107.5 (Aeroproducts,
  quoted 107) and 150.2 (Curtiss 732, quoted 150). The Curtiss 838 and 836 traces give 102
  and 81 against the quoted 120 and 98. Their b/D axes read correctly (the 838's peak is on the
  ".83" label), and the two ratios are nearly equal (1.18 and 1.21). So Curtiss probably quoted
  these two on another basis. The geometry here is as drawn.

**Efficiency against CP (9 cross-plot figures):** `symread.py` → `syms_raw.json` and
numbered overlays.
- Symbols are found with naca-archive's detector (`nacaarchive.symbols`). It also picks up
  letters in the legend and labels, and it often confuses the four shapes.
- So each data symbol's curve, and hence its J, was assigned on the overlays and zoomed crops
  by eye. Symbols the detector missed were read by eye against the grid. The fig. 6 panel of
  the Aeroproducts report was read entirely by eye.
- The values are in `../efficiency_vs_cp.csv`. Positions are good to about ±0.005 in CP and
  ±0.5 points in efficiency.
- Points on the dashed curves ("extrapolated from cross plot") are left out, except one
  symbol on each of the Aeroproducts J 3.40 and 3.50 curves.
- The Curtiss 836 figure 5 legend reads J 1.00, 1.00, 1.60. The second (□) is taken as 1.30:
  the text discusses "1.32 to 1.00".
