# Examples

This directory contains runnable example scripts and scenario folders for the
sample decay and dose workflow.

**What’s here**
- Example scripts for F71-based dose calculations.
- Example scripts for F33-based irradiation workflows.
- Scenario folders with input decks and helper data.

**Example cases**
- `Nash_flibe` — FLiBe irradiation example.
- `aluminum` — Aluminum irradiation example.
- `cobalt_steel` — Cobalt-in-steel irradiation example.
- `co_stellite` — Stellite / steel irradiation and shielding variants.
- `flibe_salt` — FLiBe salt irradiation example.
- `hotcell` — Hotcell lead-shield dose example.
- `irradiator` — Wastewater irradiator model example.
- `pipe_gas` — Offgas pipe / MHA nuclides example.
- `radiator` — Radiator decay-time scan example.
- `shielded_sample` — Shielded salt sample dose example.

Each case folder includes its own `README.md`.

**How to run**
- Top-level scripts expect to be started from inside this `examples/` directory
  (they read their SCALE inputs via `'../<file>'`), e.g.
  `cd examples && python calc_f71_doses.py`.
- Some scripts expect SCALE/ORIGEN and MAVRIC outputs to exist.

**Top-level scripts and required inputs**

Most top-level scripts read a user-supplied SCALE output file placed in the
repository root (`SCALE_FILE.f71` from a TRITON sequence, or
`SCALE_FILE.mix0007.f33` from an ORIGEN irradiation) while being launched from
`examples/`. These input files are not distributed with the repository —
generate them with your own SCALE model first (`*.f71` / `*.f33` are
gitignored). This includes `examples/msrr.f71`, used by the
`OrigenFromTritonMHA` workflow (default `'../msrr.f71'` relative to the run
directory): it is a user-supplied file, not part of the repository.

| Script | Input required (in repo root, launched from `examples/`) |
| --- | --- |
| `calc_f71_doses.py`, `calc_f71_doses_mass.py` | `../SCALE_FILE.f71` |
| `calc_f71_doses_decaytime.py` | `../SCALE_FILE.f71`, `../SCALE_FILE_60days.f71` |
| `calc_f71_storage_tank_doses.py`, `calc_f71_tank_doses.py` | `../SCALE_FILE.f71` |
| `calc_irradiation_doses.py`, `calc_irradiation_doses_decaytime.py` | `../SCALE_FILE.mix0007.f33` (F33 + atom densities) |
| `fuelsalt_doseplot.py`, `fuelsalt_doses_decaytime.py` | previously generated dose CSVs / F71 files |
| `storage_tank_doses_scan.py`, `storage_tank_doses_scan_parallel*.py` | `../SCALE_FILE.f71` |
| `decay_atoms_cvs.py` | nuclide/atoms CSV file |
| `decay_atoms_noshield.py` | none (atom densities defined in script) |
| `plot_doses.py`, `storage_tank_plots.py` | previously generated dose JSON/CSV outputs |

Scenario folders (`Nash_flibe`, `irradiator`, ...) are self-contained; see
their individual READMEs.
