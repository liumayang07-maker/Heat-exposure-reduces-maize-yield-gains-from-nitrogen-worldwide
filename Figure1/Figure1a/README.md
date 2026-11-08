# Figure 1a: mean annual EDH30, 1986-2025

[Chinese details](README.zh-CN.md) | [Environment](../../README.md) | [Sources/rights](../../RIGHTS.md)

## Bundled files and plotting

`mean_annual_EDH30_1986_2025.tif` is the direct single-band plotting input in EPSG:4326. `world_boundaries.gpkg` contains Natural Earth countries and is also used by 1c. `plot_Figure1a.R` writes the reference `Figure1a.pdf`.

EDH30 is `sum(max(hourly_temperature_C - 30, 0))`, in degrees C hours. Annual totals use UTC calendar years. The map averages available maize-masked years; valid zeros remain distinct from missing cells. This is not growing-season or hot-day-normalized exposure. Positive harvested area defines annual masks; 2021-2025 reuse 2020 area. No additional area weighting is applied to the mean map.

From the repository root:

```text
Rscript Figure1/Figure1a/plot_Figure1a.R --overwrite
```

## Optional external-input reconstruction

The three numbered Python scripts export climate data, extract area and prepare the mean. Install `requirements-upstream.txt`; obtain `Harvested_area_1961-2020.h5` (Maize group) from the provider. Authenticate Earth Engine and supply an enabled Google Cloud project:

```text
earthengine authenticate
python Figure1/Figure1a/01_export_era5land_heat_annual_tifs.py --project YOUR_PROJECT_ID
```

This submits 40 Google Drive export tasks (1986-2025), each with 12 monthly EDH30 bands; it does not download files. Download manually as `era5land_heat_monthly_YEAR.tif`. `--dry-run` still requires authentication. Do not resubmit existing tasks. `--start-year 1984` exports 42 years for 1d's preceding-December requirement.

Replace example paths and use fresh output locations:

```text
python Figure1/Figure1a/02_prepare_maize_area.py --area-h5 /path/to/Harvested_area_1961-2020.h5 --output-dir /path/to/work/area --start-year 1985
python Figure1/Figure1a/03_prepare_mean_EDH30.py --tif-dir /path/to/monthly --area-dir /path/to/work/area --output /path/to/work/mean_annual_EDH30_1986_2025.tif
```

Area extraction defaults to 1986-2020; the explicit 1985 start also supports 1d. To plot a reconstructed raster, use its expected filename in a separate package copy. Climate exports and the area HDF5 are not bundled. Current provider access and area-data redistribution permission remain unresolved. External reconstruction was not rerun for this deposit.
