# Co Stellite

Irradiation and dose examples for SS-316 with a cobalt impurity and for stellite.

**Scripts**
- `co_steel_irr_doses.py`: irradiation of SS-316 with a specified wt% of cobalt, bare sample.
- `co_steel_irr_doses_shield.py`: the same steel inside a steel pipe.
- `co_steel_irr_hcdoses_shield.py`: the steel inside a steel pipe, with contact and handling doses.
- `co_steel_irr_hcdoses_shield_pb.py`: contact and handling doses for the pipe wrapped in 0.1–10 cm of lead.
- `stell_irr_doses.py`: irradiation of stellite, bare sample.
- `steel_irr_plot_gamma.py`: plots the gamma dose of `co_steel_irr_doses.py` or `co_steel_irr_doses_shield.py`.
  Set `GEOMETRY` in the script to `'bare'` or `'shield'` to match the run.
- `steel_irr_plot_hcdoses_gamma.py`: plots the gamma contact and handling doses of `co_steel_irr_hcdoses_shield.py`.
- `pb_steel_irr_plot_hcdoses_gamma.py`: fits the lead-thickness scan of `co_steel_irr_hcdoses_shield_pb.py` at the
  longest decay time and prints the lead thickness that brings the dose down to 80 mrem/h.
- `stell_irr_plot_gamma.py`: plots the gamma dose of `stell_irr_doses.py`.

**Case directory**

Every script reads the sample mass and the irradiation time from the name of the current directory.
The directory name must contain `dose-<years>year_` and `_<mass>g`, for example `dose-2year_10g`.
Run the calculation script and its plot script from the same case directory, for example:

```
mkdir -p examples/co_stellite/dose-2year_10g
cd examples/co_stellite/dose-2year_10g
export STELLITE_SCALE_OUT=/path/to/triton/msrr.out
python ../co_steel_irr_doses.py --wtpctCo 0.1
python ../steel_irr_plot_gamma.py
```

**Inputs** (user-supplied, placed in the parent of the case directory, `examples/co_stellite/` above)
- `steel.f33`: ORIGEN library for the steel, read as `'../steel.f33'` by the four `co_steel_*` scripts.
- `stellite.f33`: ORIGEN library for the stellite, read as `'../stellite.f33'` by `stell_irr_doses.py`.
- TRITON output with the irradiation flux of mixture 8140 (steel) and mixture 9000 (stellite).
  The scripts read its path from the `STELLITE_SCALE_OUT` environment variable.
  The default is `../msrr.out`.
  Use the same value for a calculation script and its plot script, since the plot titles quote this flux.
- SCALE inventory files (`*.f71` / `*.f33`) are gitignored.

**Extra dependencies**
- `MSRRpy` builds the SS-316 plus cobalt composition (`MSRRpy.mat.material_types.Solid`) in the four
  `co_steel_*` scripts. Install it separately.
- `scipy` is used by `pb_steel_irr_plot_hcdoses_gamma.py` for the fit (`pip install '.[examples]'`).

**Detector positions**
- `co_steel_irr_doses.py`, `stell_irr_doses.py`: 30 cm from the sample centre (`DoseEstimator`).
- `co_steel_irr_doses_shield.py`: 30 cm from the pipe outer surface (`det_standoff_distance = 30.0`).
- `co_steel_irr_hcdoses_shield*.py`: contact detector 0.1 cm and handling detector 30 cm from the outer surface of
  the last layer.

**Outputs** (in the case directory)
- `responses.json`: all dose responses per decay time. The plot scripts read this file.
- `doses.json`: total dose per decay time, and per lead thickness for the `_pb` script.
  For the two `hcdoses` scripts each entry holds `{'contact': ..., 'handling': ...}`, the separate contact and
  handling totals.
