# SKILLS.md

Project-specific working skills for `sample-decay-dose`.

## 1) Core Library Changes

Use when modifying dose/decay logic in `sample_decay_dose/`.

- Read affected module and nearby tests first.
- Preserve public API behavior unless change is requested.
- Add/adjust unit tests in `test/` for all logic changes.
- Run: `pytest -q`

## 2) Example Workflow Changes

Use when updating scripts in `examples/`.

- Keep scripts runnable as standalone examples.
- Prefer clear input/output naming; avoid overwriting unrelated artifacts.
- Validate quickly with a reduced or existing example input when possible.

## 3) Leaky-Box Utility Changes

Use when editing `leaky_box_origen/`.

- Confirm whether SCALE runtime is required before executing long jobs.
- Separate code edits from expensive data-generation runs.
- If OCR/PDF extraction is involved, note extra tools: `pdftoppm`, `tesseract`.

## 4) Test and Regression Skill

Use before finalizing non-trivial changes.

- Run targeted tests first, then full suite:
  - `pytest -q test`
- If outputs are numerical, check for unit consistency and order-of-magnitude sanity.

## 5) Documentation Skill

Use when behavior, setup, or workflow changes.

- Update relevant README sections with exact commands.
- Keep examples synchronized with current filenames and expected outputs.

