# Co Stellite

Example case under `examples/`.

**Scripts**
- `co_steel_irr_doses.py` — Example use case of SampleDose - irradiation of SS-316 with a specified wt% of cobalt
- `co_steel_irr_doses_shield.py` — Irradiation of SS-316 with a specified wt% of cobalt in a steel pipe
- `co_steel_irr_hcdoses_shield.py` — Irradiation of SS-316 with a specified wt% of cobalt in a steel pipe, handling and contact doses
- `co_steel_irr_hcdoses_shield_pb.py` — Irradiation of SS-316 with a specified wt% of cobalt in a steel pipe, handling and contact doses
- `pb_steel_irr_plot_hcdoses_gamma.py` — Finds minimum lead thickness for shielding
- `steel_irr_plot_gamma.py` — Plotting script for calc_doses_mass
- `steel_irr_plot_hcdoses_gamma.py` — Plotting script for calc_doses_mass
- `stell_irr_doses.py` — Example use case of SampleDose - irradiation of stellite
- `stell_irr_plot_gamma.py` — Plotting script for calc_doses_mass

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `steel.f33` in `examples/` (read as `'../steel.f33'`).
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
