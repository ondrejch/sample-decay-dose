# AGENTS.md

Project agent guidance for `sample-decay-dose`.

## Purpose

This repository computes dose rates from decaying samples using SCALE/ORIGEN/MAVRIC workflows and post-processing utilities.

## Repo Map

- `sample_decay_dose/`: core Python library (`SampleDose.py`, helpers, constants, data).
- `examples/`: runnable scripts for common workflows and plotting.
- `leaky_box_origen/`: leaky-box ORIGEN utilities and related scripts.
- `concrete_irrad/`: concrete-lid activation + MAVRIC dose how-to and the concrete/rebar mixed-layer homogenizer (`concrete_rebar_mixer.py`).
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

## Scientific Writing Style

Applies to the paper drafts, summaries, and any prose written
for readers outside a single session.

1. **Sentences over dashes.** Write parenthetical material as its own sentence
   instead of an em-dash aside. En-dashes are allowed for numeric ranges only
   (`5--10\%`). If a sentence needs three or more commas to hold together,
   split it.
2. **Positive first.** Say what a quantity, model, or result *is* and what it
   *does*, with the concrete object named (observable, number, mechanism).
   Put any boundary in a short follow-up sentence. A stack of "is not / does
   not / rather than" clauses is a signal to restructure.
3. **No contrast slogans.** Avoid "X is this. Not that." shapes and their
   inline cousins ("A, not B"). Use one declarative sentence, or two
   sentences: claim first, limitation second.
4. **Negation has a job.** Keep "not/never" for genuine exclusions and scope
   limits: falsifier boundaries, applied cuts ("p+A is excluded from the
   quantitative tables"), refused claims. Do not use negation as the default
   way to define or praise something.
5. **Concrete beats abstract.** Name the observable, the value, and the
   mechanism before the interpretation. When editing, prefer the version a
   reader can picture.
6. **Professional register.** No colloquial competitive language ("win",
   "lose", "beat", "champion"). State the metric outcome directly: "has the
   lowest figure of merrit in five of six runs", "is ranked first by log-pull RMS",
   "agreement improves under the corrected convention". Model comparisons are
   reported as measurements, not contests.

