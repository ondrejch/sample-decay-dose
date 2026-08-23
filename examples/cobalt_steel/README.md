# Cobalt Steel

Example case under `examples/`.

**Scripts**
- `co_steel_irr_doses.py` — Example use case of SampleDose - irradiation of SS-316 with a specified wt% of cobalt
- `co_steel_irr_plot_gamma.py` — Plotting script for calc_doses_mass

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `SCALE_FILE.mix0007.f33` in `examples/` (read as `'../SCALE_FILE.mix0007.f33'`).
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
