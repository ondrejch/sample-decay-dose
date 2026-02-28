# Leaky Box (ORIGEN)

Leaky-box utilities built around SCALE/ORIGEN. This is where leakage tests and
F71-based leakage simulations live.

**Primary entrypoints**
- `LeakyBox.py` is the maintained implementation.
- `BoxF71.py` is a legacy variant kept in sync for reference.
- `AnalysisRerun.py` rebuilds spreadsheets from saved JSON outputs.

**Outputs**
- Runs create `_box_*` working directories under this folder.
- JSON, Excel, and plots are written to the current working directory.

**Notes**
- Defaults (F71 path, case/time selection, leak elements) are set in `LeakyBox.py`.
