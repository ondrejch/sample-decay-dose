# Paper support scripts

These scripts regenerate the numbers and figures that the leaky-box journal paper quotes. The paper source
(`../journal_paper/`) is a separate Overleaf clone that this repository ignores. Both scripts write into the run
directories they read, which git also ignores.

| Script | Reads | Produces | Paper items |
|---|---|---|---|
| `xe_convergence.py [run\|analyse\|all]` | `run_2026-10-23/leaky_boxes_xe13*_spd*.json` | `convergence.csv`, `conservation.csv`, `figs/leaky_boxes_xe135_convergence.png` | Xe-135 step sweep, Xe-136 conservation and its trapezoid excess (Section 5, Appendix A) |
| `psar_numbers.py [section ...]` | `run_2026-10-22/` (EOL and MHA runs, `msrr_psar_eol.f71`) | `paper_numbers_psar.json`, `paper_mha_attribution.csv` | Part II numbers listed below |

- `xe_convergence.py run` repeats the ORIGEN runs (SCALE, about an hour). The default `analyse` phase works from the
  JSON snapshots the runs leave, so it needs no SCALE.
- `psar_numbers.py` runs no ORIGEN case; it reads the F71 files of the runs with `obiwan` (a few minutes). Sections:
  - `totals`: the MHA and EOL doses, integrated from t = 0.
  - `attribution`: per-nuclide MHA doses at the cell exit (Box B) and in the systems-pit release (Box C) for
    FGR-11 + FGR-12 and FGR-11 alone, and the restriction to the PSAR mass chains.
  - `two-hour`: the 2-hour EOL doses, integrated over the internal positions of the Box B and Box C sub-runs. F71
    times count from the start of each sub-run, so the script anchors them to the end time in the box json.
  - `positions`: dose-driver activities at F71 positions 2 and 6 and the share-weighted dose change.
  - `pit-rate`: the systems-pit leak rate that reproduces the published 48.66 mrem. Each contributor runs as a
    parent-free three-box chain with its own half-life (short-lived daughters with their parent), scaled to its
    ORIGEN pit dose at 1 %/d. The estimate reproduces the ORIGEN cell-exit dose as the pit rate grows without
    bound, which checks it at the other end of the range.
  - `daughters`: an upper estimate of the Cs, Rb and Sr daughters of the liberated inventory that PSAR Eq. 13.1-1
    transports and this model keeps in the salt.
  - `hand-check`: the analytic three-box I-131 dose, with I-131 bred from Te-131m and the fast A = 131 precursors.
  - `iodine-vapour`: EOL and MHA doses with the ICRP-72 Table A.3 elemental iodine vapour coefficients.
