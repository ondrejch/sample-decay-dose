# Irradiator

Example case under `examples/`.

**Scripts**
- `irradiator.py`: wastewater irradiator using SampleDose, with an irradiated SS-316 plate as the source.

**Inputs**
- `EIRENE.mix0002.f33`: user-supplied ORIGEN library next to the script.
- SCALE inventory files (`*.f71` / `*.f33`) are gitignored.

**Outputs**
- `run_irradiator_1m/`: ORIGEN irradiation and decay of the plate, created by the run.
- `run_irradiator_1m_MAVRIC/`: MAVRIC dose calculation, created by the run.

**Run**
- From this directory: `python irradiator.py`
