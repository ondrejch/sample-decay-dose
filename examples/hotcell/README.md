# Hotcell

Example case under `examples/`.

**Scripts**
- `leadcell-plots.py` — Plotting helper for the hotcell lead-shield example.
- `leadcell.py` — Example hotcell is 70 x 70 x 61 cm, close enough to 70 x 70 x 70 cm

**Run**
- From this directory: `python <script>.py`
- Needs a user-supplied `msrr.f71` next to the script.
- SCALE inventory files (`*.f71` / `*.f33`) are user-supplied and gitignored.
- `leadcell.py` computes the MAVRIC adjoint flux for every case (`reuse_adjoint_flux = False`).
  Enable the reuse only when an adjoint flux file from a run with the same geometry and detectors exists.
- Outputs: `responses.json` with the contact and handling responses per decay day, read by `leadcell-plots.py`.
