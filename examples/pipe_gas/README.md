# Pipe Gas

Example case under `examples/`.

**Scripts**
- `dose_plot.py` — Plotting script for pipe_doses.py
- `immediate_doses.py` — Prints dose at t=0
- `pipe_doses.py` — Irradiation of SS-316 with a specified wt% of cobalt in a steel pipe

**Data / inputs**
- `msrr.f71` — input or output data file
- `responses.json` — input or output data file

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `msrr.f71` in `examples/` (read as `'../msrr.f71'`).
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
