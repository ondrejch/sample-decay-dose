# Irradiator

Example case under `examples/`.

**Scripts**
- `irradiator.py` — Wastewater irradiator using SampleDose - irradiation of SS316

**Data / inputs**
- `run_*` — run directories with SCALE/MAVRIC outputs
- `irradiator/` — case data or helper files
- `EIRENE.mix0002.f33` — input or output data file

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `EIRENE.mix0002.f33` next to the script.
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
