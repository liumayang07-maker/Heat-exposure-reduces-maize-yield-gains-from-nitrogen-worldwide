# Figure 1b: continental annual EDH30

[Chinese details](README.zh-CN.md) | [Environment](../../README.md) | [Sources/rights](../../RIGHTS.md)

## Bundled files

- `continent_annual_EDH30_1986_2025.csv`: 240 rows, six continents by 40 years, unsmoothed area-weighted exposure and coverage.
- `plot_Figure1b.R`: six separate GLS AR(1)/REML fits, writing `Figure1b.pdf` and `GLS_models.csv`.
- `GLS_models.csv`: slopes, standard errors, intervals, raw/BH-adjusted P and AR(1) parameters.
- `Figure1b_Source_Data.xlsx`: five sheets: `Annual_data`, `Plot_data`, `GLS_models`, `Model_parameters`, `Coefficient_covariance`.
- `prepare_continent_EDH30.py`, `maize_continent_crosswalk.csv.gz`: upstream aggregation and fixed original assignments (840,477 rows).
- `Figure1b.pdf`: archived panel reference without final manual annotations.

The plotted series is not smoothed. Shading is the pointwise 95% confidence interval for the GLS fitted mean. Blank crosswalk assignments are excluded. Area weighting uses the corresponding year's area, with 2020 reused in 2021-2025.

From the repository root:

```text
Rscript Figure1/Figure1b/plot_Figure1b.R --overwrite
python Figure1/build_source_data.py --rscript Rscript --overwrite
```

For external reconstruction, obtain monthly climate files and annual area tables as described for 1a, then run:

```text
python Figure1/Figure1b/prepare_continent_EDH30.py --tif-dir /path/to/monthly --area-dir /path/to/work/area --output /path/to/work/continent_annual_EDH30_1986_2025.csv
```

Aggregation rejects an existing output. To fit the new CSV, place it under its expected filename in a separate package copy. The crosswalk contains IDs/labels rather than original geometry or area values; this deposit does not grant a new blanket licence for it.
