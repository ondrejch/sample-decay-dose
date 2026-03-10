# AGENS.md

Project agent guidance for `sample-decay-dose`.

## Purpose

This repository computes dose rates from decaying samples using SCALE/ORIGEN/MAVRIC workflows and post-processing utilities.

## Repo Map

- `sample_decay_dose/`: core Python library (`SampleDose.py`, helpers, constants, data).
- `examples/`: runnable scripts for common workflows and plotting.
- `leaky_box_origen/`: leaky-box ORIGEN utilities and related scripts.
- `test/`: unit tests.
- Root `*.json5`, `*.csv`, `*.png`, `*.xlsx`: generated artifacts from sample runs.

## Environment

- Python `>=3.10`
- Install deps:
  - `pip install -r requirements.txt`
  - or `pip install ./`
- SCALE binaries must be available via `SCALE_BIN`, for example:
  - `export SCALE_BIN=/opt/scale6.3.2-mpi/bin`

## Validation Commands

- Preferred: `pytest -q`
- Alternative: `PYTHONPATH=. python -m unittest discover -s test -p 'test_*.py'`

## Agent Working Rules

- Keep changes focused and minimal; avoid touching generated artifacts unless explicitly requested.
- Prefer editing source under `sample_decay_dose/`, `examples/`, `leaky_box_origen/`, and `test/`.
- If a change affects physics assumptions or units, document the assumption in code comments and/or README.
- Add or update tests when behavior changes.
- When running heavy SCALE-dependent scripts, confirm intent before long runs.

## Typical Outputs

- Leaky-box runs: `boxA*.json5`, `boxB*.json5`, `boxC*.json5`, `leaky_boxes*.json/.xlsx/.png`
- Dose post-processing: `leaky_boxes_dose*.csv/.png`

