# Leaky Box (ORIGEN)

Leaky-box utilities built around SCALE/ORIGEN. This is where leakage tests and
F71-based leakage simulations live.

**Primary entrypoints**
- `LeakyBox.py` is the maintained implementation.
- `BoxF71.py` is a legacy variant kept in sync for reference.
- `AnalysisRerun.py` rebuilds spreadsheets from saved JSON outputs.

**Outputs**
- Runs write all artifacts under `leaky_box_origen/run_YYYY-MM-DD/`.
- ORIGEN working directories (`_box_*`), JSON, Excel, CSV, and plots are all written to that run directory.

**Notes**
- Defaults (F71 path, case/time selection, leak elements) are set in `LeakyBox.py`.
- Override F71 defaults with environment variables:
  - `LEAKYBOX_F71_PATH`
  - `LEAKYBOX_F71_CASE`
- Source PDFs used for DCF extraction live under `PDF/`.
- Script data CSV files (DCF tables) live under `leaky_box_origen/data/`.
