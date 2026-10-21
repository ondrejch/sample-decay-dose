# Radiator

Example case under `examples/`.

**Scripts**
- `radiator_doses.py` — contact and handling doses from the salt in a tube-bundle radiator, read from one
  F71 position (case 2, end of the burn) without a decay step

**Data / inputs**
- `origen.f71` — input or output data file

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `origen.f71` next to the script.
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
