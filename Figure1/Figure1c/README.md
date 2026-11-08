# Figure 1c: grid-level EDH30 trends

[Chinese details](README.zh-CN.md) | [Environment](../../README.md) | [Sources/rights](../../RIGHTS.md)

## Bundled files

- `EDH30_slope_per_decade.tif`: plotting slopes (degrees C hours per decade).
- `EDH30_HAC_BH_adjusted_pvalue.tif`: map significance input.
- `Figure1c_Source_Data.csv.gz`: pixel coordinates, slope, HAC standard error, pointwise 95% interval, raw/BH-adjusted P.
- `prepare_grid_EDH30_trends.py`: OLS/HAC/BH analysis; HAC code also used by 1d.
- `plot_Figure1c.py`, `Figure1c.pdf`: plotting entry and reference panel.

Only complete 1986-2025 records enter analysis. HAC uses three lags and a finite-sample correction; BH adjustment is across all tested complete cells. Plotting slopes are float32 and source estimates float64, so equivalent calculation orders differ in final digits; the archived slope comparison uses rtol=1e-6, atol=1e-5. Source-table columns contain more model information than the plotting TIFFs.

From the repository root, also retaining 1a's GeoPackage:

```text
python Figure1/Figure1c/plot_Figure1c.py --data-dir Figure1/Figure1c --output Figure1/Figure1c/Figure1c.pdf --overwrite
```

For optional reconstruction, first obtain 1a's external inputs:

```text
python Figure1/Figure1c/prepare_grid_EDH30_trends.py --tif-dir /path/to/monthly --area-dir /path/to/work/area --output-dir /path/to/work/grid_trends
python Figure1/Figure1c/plot_Figure1c.py --data-dir /path/to/work/grid_trends --output /path/to/work/grid_trends/Figure1c.pdf
```

Use a new output directory. The alternative `--annual-dir` requires original maize-masked annual TIFFs named `gee_era5land_annual_heat_maize_area_YEAR.tif`, containing `annual_EDH30_YEAR` and tag `maize_mask_applied=True`. Arbitrary unmasked grids are not accepted. Full upstream reconstruction was not rerun for this deposit.
