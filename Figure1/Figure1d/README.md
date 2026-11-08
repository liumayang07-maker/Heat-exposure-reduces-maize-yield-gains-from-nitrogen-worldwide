# Figure 1d: regional summer EDH30

[Chinese details](README.zh-CN.md) | [Environment](../../README.md) | [Sources/rights](../../RIGHTS.md)

## Bundled files

- `region_summer_1985_2025.csv`: 615 rows, 15 regions by 41 years, with coverage and season fields.
- `regions.csv`: original region definitions and order.
- `maize_region_crosswalk.csv.gz`: 840,477 original grid assignments, including blanks; only used upstream.
- `prepare_regional_summer_EDH30.py`: seasonal area-weighted aggregation.
- `plot_Figure1d.py`, `plot_heatmap.R`: analysis/plotting entries.
- `Figure1d_heatmap.pdf`, `Figure1d_trend.pdf`: reference panels.
- `Figure1d_Source_Data.xlsx`: `Annual_data`, `Heatmap_data`, `Trend_models`, `Model_parameters`, `Coefficient_covariance`.

The heatmap uses 1985-2025 (41 years), within-region sample-SD Z-scores and Euclidean/Ward.D2 clustering. Trends use raw exposure in 1986-2025 (40 years), OLS with three-lag HAC uncertainty and BH adjustment across 15 slopes. Intervals are pointwise 95% intervals. Preserve these different year ranges.

Summer is June-August for assigned northern regions and preceding December-February for Brazil, Argentina and South Africa, labelled by ending year. Summer 1985 requires December 1984. Original whole-region season assignments remain unchanged. Exposure is area-weighted cumulative EDH30, not exposure divided by hot-day count.

From the repository root:

```text
python Figure1/Figure1d/plot_Figure1d.py --rscript Rscript --overwrite
```

Retain `Figure1c/prepare_grid_EDH30_trends.py` and `Figure1/build_source_data.py`, which provide shared functions. The entry writes two PDFs and the workbook; temporary tables are cleaned automatically.

For optional external reconstruction, obtain 1984-2025 monthly exports (1a exporter with `--start-year 1984`) and 1985-2020 area tables (extractor with `--start-year 1985`), then run:

```text
python Figure1/Figure1d/prepare_regional_summer_EDH30.py --tif-dir /path/to/monthly --area-dir /path/to/work/area --output /path/to/work/region_summer_1985_2025.csv
python Figure1/Figure1d/plot_Figure1d.py --input /path/to/work/region_summer_1985_2025.csv --output-dir /path/to/work/plots --rscript Rscript
```

The fixed regional assignment derives from GADM 4.1 countries and Natural Earth Chinese admin-1 units. No original GADM geometry is distributed. The derived crosswalk is not declared unrestricted; clarify applicable terms before open release. Do not silently substitute different boundaries. Original area HDF5 and climate exports remain external.
