# Pipe Gas

Contact and handling doses from MHA off-gas nuclides inside a 1/16" OD SS-316H tube, for 0 to 1 day of decay.
The gas is the Te, I, Xe, Br, Kr and H retained from the MSRR fuel inventory (`OrigenFromTritonMHA`).
`MHATank` scales it to about 1.0 mg and fills the bore of a 10 ft tube.

**Scripts**
- `pipe_doses.py`: runs the ORIGEN decay and the MAVRIC contact and handling dose calculations.
- `dose_plot.py`: plots the gamma contact and handling doses of `pipe_doses.py` versus decay time.
- `immediate_doses.py`: prints the contact and handling doses of `pipe_doses.py` as LaTeX table rows.

**Inputs**
- `msrr.f71`: user-supplied MSRR TRITON inventory in `examples/`, read as `'../msrr.f71'` when the scripts run
  from this directory.
- SCALE inventory files (`*.f71` / `*.f33`) are gitignored.

**Detector positions**
- Contact detector 0.1 cm and handling detector 30 cm from the tube outer surface.

**Outputs**
- `responses.json`: all dose responses per decay time. `dose_plot.py` and `immediate_doses.py` read this file.
- `doses.json`: `{'contact': ..., 'handling': ...}` dose totals per decay time.

**Run**
- From this directory: `python pipe_doses.py`, then `python dose_plot.py` or `python immediate_doses.py`.
