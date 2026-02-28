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
- Run from the repo root, e.g. `python examples/calc_f71_doses.py`.
- Some scripts expect SCALE/ORIGEN and MAVRIC outputs to exist.
