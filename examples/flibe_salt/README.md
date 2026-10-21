# Flibe Salt

Example case under `examples/`.

**Scripts**
- `flibe_doseplot.py` — Plotting gamma doses
- `flibe_irr_dose.py` — Example use case of SampleDose - irradiation of natural FLiBe (`my_flibe_nat`);
  a Li-7 enriched composition `my_flibe_li7` is defined for comparison

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `SCALE_FILE.mix0001.f33` in `examples/` (read as `'../SCALE_FILE.mix0001.f33'`).
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
