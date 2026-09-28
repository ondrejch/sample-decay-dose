# Cobalt Steel

Irradiation of an SS-316 coupon (1 g, 2 years at 1e13 n/cm2/s) for several cobalt impurity levels.

**Scripts**
- `co_steel_irr_doses.py`: irradiation of SS-316 with the wt% of cobalt given by the required `--wtpctCo` option.
  The dose detector sits 30 cm from the sample centre.
- `co_steel_irr_plot_gamma.py`: plots the gamma doses of the `co_steel_irr_doses.py` runs at 0.02, 0.05, 0.10,
  0.20 and 0.40 wt% cobalt.

**Run**

Run each cobalt level in its own subdirectory named `co_<wt%>pct`, then plot from this directory:

```
cd examples/cobalt_steel
for co in 0.02 0.05 0.10 0.20 0.40; do
    mkdir -p co_${co}pct
    (cd co_${co}pct && python ../co_steel_irr_doses.py --wtpctCo ${co})
done
python co_steel_irr_plot_gamma.py
```

**Inputs**
- `SCALE_FILE.mix0007.f33`: user-supplied ORIGEN library placed in `examples/cobalt_steel/`.
  The script reads it as `'../SCALE_FILE.mix0007.f33'` from inside a `co_<wt%>pct/` subdirectory.
- SCALE inventory files (`*.f71` / `*.f33`) are gitignored.

**Extra dependencies**
- `MSRRpy` builds the SS-316 plus cobalt composition (`MSRRpy.mat.material_types.Solid`) in
  `co_steel_irr_doses.py`. Install it separately.

**Outputs**
- `co_<wt%>pct/responses.json`: all dose responses per decay time, read by `co_steel_irr_plot_gamma.py`.
- `co_<wt%>pct/doses.json`: total dose per decay time.
