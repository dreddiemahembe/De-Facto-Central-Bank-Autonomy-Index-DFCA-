# Changes

## 1.1 (October 2026)

The technical note describes version 1.0. These are the differences.

- **The geometric mean is now taken over autonomy scores (1 − capture), not over capture scores.**
  The note's reason for the geometric mean is that a bank should not offset heavy monetisation
  with a clean record elsewhere. Taken over capture scores, a clean dimension pulls the capture
  average down, which is the offsetting the note rules out. When version 1.0 was run on the
  simulated panel, a bank captured on three of four dimensions scored DFCA 0.90 and the median
  country-year scored 0.94. Version 1.1 scores the first at 0.09 and the median at 0.56.
- **Geometric aggregation shifts scores into [0.1, 1] and back.** Version 1.0 floored scores at
  0.0001, so one dimension at its worst value pulled the whole index towards zero. The shift keeps
  0 at 0 and 1 at 1. The floor is the argument `geo_floor`; 0.1 is a choice, and the README shows
  how much it matters.
- **`dims` is an argument.** The same code builds the reduced core and other indices.
- **Rank bands report how many countries were ranked.** In the simulated test 16 of 20 countries
  have a 2015 value, so the intervals are positions out of 16, not out of 20.
- **Added:** input checks, `coverage_audit()`, `coverage_matrix()`, `impute_within_country()`,
  `level_bands()`, a within-country variance share in `diagnostics()`, 25 pass/fail checks in the
  test harness, `demo_figures.R`.

Numbers in the technical note that do not change: the dimension correlations (0.28 to 0.49),
PC1 share (53%), winsorised maximum (90.5 to 8.8) and the listwise loss (511 of 780). Numbers that
do: recovery of the latent factor (0.64 becomes 0.70), the rank interval width (3.5 out of 20
becomes 4 out of 16) and the validation regression.
