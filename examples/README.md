# Examples

This directory contains runnable example scripts and scenario folders for the
sample decay and dose workflow.

**What’s here**
- Example scripts for F71-based dose calculations.
- Example scripts for F33-based irradiation workflows.
- Scenario folders with input decks and helper data.

**Example cases**
- `Nash_flibe`: FLiBe irradiation example.
- `aluminum`: Aluminum irradiation example.
- `cobalt_steel`: Cobalt-in-steel irradiation example.
- `co_stellite`: Stellite / steel irradiation and shielding variants.
- `flibe_salt`: FLiBe salt irradiation example.
- `hotcell`: Hotcell lead-shield dose example.
- `irradiator`: Wastewater irradiator model example.
- `pipe_gas`: Offgas pipe / MHA nuclides example.
- `radiator`: Radiator decay-time scan example.
- `shielded_sample`: Shielded salt sample dose example.

Each case folder includes its own `README.md`.

**How to run**
- Top-level scripts expect to be started from inside this `examples/` directory
  (they read their SCALE inputs via `'../<file>'`), e.g.
  `cd examples && python calc_f71_doses.py`.
- Some scripts expect SCALE/ORIGEN and MAVRIC outputs to exist.
- ORIGEN and MAVRIC cases run in subdirectories of the launch directory
  (`run_*/`, `_decaybox_*/`). These directories are gitignored.

**Extra dependencies**
- `scipy` is used by `co_stellite/pb_steel_irr_plot_hcdoses_gamma.py`. Install it with `pip install '.[examples]'`.
- `MSRRpy` builds the SS-316 plus cobalt compositions in the `co_steel_*` scripts of
  `co_stellite/` and `cobalt_steel/`. Install it separately.

**Detector positions**
- `DoseEstimator` places the point detector `det_x` = 30 cm from the sample centre.
- The tank estimators (`DoseEstimatorSquareTank`, `DoseEstimatorStorageTank`,
  `DoseEstimatorGenericTank`) place it `det_standoff_distance` from the outer
  surface of the last shielding layer (1 cm by default). They recompute `det_x`
  in `mavric_deck()`, so set `det_standoff_distance` to move the detector.
- `HandlingContactDoseEstimatorGenericTank` and `MHATank` use a contact detector
  0.1 cm and a handling detector 30 cm from that surface. Their `contact_dose`
  and `handling_dose` properties hold the two totals.

**Top-level scripts, inputs and outputs**

Most top-level scripts read a user-supplied SCALE output file placed in the
repository root (`SCALE_FILE.f71` from a TRITON sequence, or
`SCALE_FILE.mix0007.f33` from an ORIGEN irradiation) while being launched from
`examples/`. Generate these input files with your own SCALE model first. The
repository does not distribute them, and `*.f71` / `*.f33` are gitignored.
The MSRR inventory conventionally kept at `examples/msrr.f71` is also
user-supplied. The `OrigenFromTritonMHA` workflow reads it through its default
path `'../msrr.f71'`, which resolves against the launch directory. For example,
launching `pipe_gas/pipe_doses.py` from `examples/pipe_gas/` reads
`examples/msrr.f71`.

Each top-level script writes its results to its own file, so runs of different
scripts in `examples/` keep their results side by side. The generated
`responses*.json` and `doses*.json` files are gitignored.

| Script | Reads (launched from `examples/`) | Writes |
| --- | --- | --- |
| `calc_f71_doses.py` | `../SCALE_FILE.f71` | printed doses |
| `calc_f71_doses_mass.py` | `../SCALE_FILE.f71` | `responses_f71_mass.json`, `doses_f71_mass.json` |
| `calc_f71_doses_decaytime.py` | `../SCALE_FILE.f71` | `responses_f71_decaytime.json`, `doses_f71_decaytime.json` |
| `calc_f71_storage_tank_doses.py`, `calc_f71_tank_doses.py` | `../SCALE_FILE.f71` | printed doses |
| `calc_irradiation_doses.py` | `../SCALE_FILE.mix0007.f33` | printed doses |
| `calc_irradiation_doses_decaytime.py` | `../SCALE_FILE.mix0007.f33` | `responses_irr_decaytime.json`, `doses_irr_decaytime.json` |
| `fuelsalt_doses_decaytime.py` | `../SCALE_FILE_60days.f71` (60 EFPD burn) | `responses_fuelsalt_decaytime.json`, `doses_fuelsalt_decaytime.json` |
| `fuelsalt_doseplot.py` | `responses_fuelsalt_decaytime.json` | `dose_g_fuelsalt_1MW*.png` |
| `storage_tank_doses_scan.py` | `../SCALE_FILE.f71` | `doses_storage_tank_scan.json` |
| `storage_tank_doses_scan_parallel.py` | `../SCALE_FILE.f71` | `doses_storage_tank_scan_parallel.json` |
| `storage_tank_doses_scan_parallel_noshield.py` | `../SCALE_FILE.f71` with `--run-analysis`, then `doses_storage_tank_noshield.json` | `doses_storage_tank_noshield.json` |
| `storage_tank_plots.py` | `<dir>/doses_storage_tank_scan_parallel.json` for each result folder in `my_data` | `dose_g_storage_tank-*.png`, `storage_dose.xlsx` |
| `decay_atoms_cvs.py` | nuclide CSV (`--atoms-csv`, default `/tmp/atoms.csv`) with `--run-analysis`. Without options it plots `<dir>/doses_decay_atoms_cvs.json` for each result folder in `plot()` | `doses_decay_atoms_cvs.json`, or `dose_g_offgas_tank-*.png` and `offgas_dose.xlsx` |
| `decay_atoms_noshield.py` | nuclide CSV (`--atoms-csv`, default `/tmp/atoms.csv`) with `--run-analysis`, then `doses_decay_atoms_noshield.json` | `doses_decay_atoms_noshield.json` |
| `plot_doses.py` | `responses_f71_mass.json`, `responses_f71_decaytime.json` or `responses_irr_decaytime.json`, selected by `LABEL` | `dose_mass*.png` or `dose_decay_time*.png` |

Scenario folders (`Nash_flibe`, `irradiator`, ...) are self-contained; see
their individual READMEs.
