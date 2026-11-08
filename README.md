<a id="overview"></a>

# Heat exposure reduces maize yield gains from nitrogen worldwide

Data and analysis code accompanying the manuscript with this working title. This deposit contains **Figures 1–4 and Extended Data Figure 3**, including their deposited supplementary analyses. The repository is private; editor/reviewer access must be arranged separately.

Start with [Figure 1](#figure-1). See [sources and rights](#sources-and-rights) for attribution and applicable access conditions. No publication or archive DOI has been assigned here.

[中文说明 / Chinese README](README.zh-CN.md)

## Contents

- [Environment](#environment)
- [Figure 1 overview and commands](#figure-1)
- [Figure 1a](#figure-1a)
- [Figure 1b](#figure-1b)
- [Figure 1c](#figure-1c)
- [Figure 1d](#figure-1d)
- [Figure 2 files and commands](#figure-2)
- [Figure 3 files and commands](#figure-3)
- [Figure 4 files and commands](#figure-4)
- [Extended Data Figure 3](#extended-data-figure-3)
- [Data sources and rights](#sources-and-rights)

<a id="environment"></a>

## Environment

Tested runtime: Python 3.12.5 and R 4.5.3 on Windows. The tested direct Python dependencies and optional upstream dependencies are listed below. Use a separate environment:

```text
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install numpy==2.4.4 pandas==2.3.3 scipy==1.17.1 statsmodels==0.14.6 matplotlib==3.10.8 rasterio==1.5.0 geopandas==1.1.3
```

For optional reconstruction from external climate and area inputs, also install:

```text
python -m pip install h5py==3.16.0 earthengine-api==1.7.42
```

In R, install required packages if unavailable:

```r
install.packages(c("terra", "sf", "ggplot2", "data.table", "nlme", "pheatmap"))
```

Verified R packages: terra 1.9.1, sf 1.1.0, ggplot2 4.0.2, data.table 1.18.2.1, nlme 3.1.168, pheatmap 1.0.13. These are recorded versions, not an R lockfile. R needs Cairo PDF support; sf/terra need their usual system geospatial libraries. Python plots use Arial; missing fonts can change appearance. Put `Rscript` on PATH or supply its path to entries accepting `--rscript`.

Scripts do not install packages automatically. Fresh environment installation and other operating systems have not been validated in this deposit.

## Reproduction levels

1. **Bundled inputs:** regenerate five panel PDFs and six GLS model rows; source workbooks are supplied outputs, without downloads or the original working directory.
2. **External inputs:** obtain the area HDF5 and monthly climate exports, then run preparation scripts. This requires external access and potentially substantial compute; it was not rerun for this deposit.
3. **Final layout:** the assembled manuscript figure was manually edited. Reproducing individual panels does not recreate that layout.

No open-source or blanket open-data licence has been selected. Third-party terms remain applicable; see [Sources and rights](#sources-and-rights).

<a id="figure-1"></a>

## Figure 1: heat exposure across global maize systems

Scripts, archived inputs, numerical results and reference PDFs for panels 1a-d. Preserve the subdirectory structure: panels share code and geometry.

| Panel | Entry | Essential inputs | Outputs |
| --- | --- | --- | --- |
| [1a](#figure-1a) | `Figure1a/plot_Figure1a.R` | mean TIFF and Natural Earth GeoPackage | PDF |
| [1b](#figure-1b) | `Figure1b/plot_Figure1b.R` | continental annual CSV | PDF and GLS model CSV |
| [1c](#figure-1c) | `Figure1c/plot_Figure1c.py` | slope/P TIFFs and 1a geometry | PDF |
| [1d](#figure-1d) | `Figure1d/plot_Figure1d.py` | regional annual CSV and region definitions | two PDFs |

### Run from the repository root

Archived outputs already exist. These commands explicitly replace them in your checkout. Run in a separate copy to preserve the references.

```text
Rscript Figure1/Figure1a/plot_Figure1a.R --overwrite
Rscript Figure1/Figure1b/plot_Figure1b.R --overwrite
python Figure1/Figure1c/plot_Figure1c.py --data-dir Figure1/Figure1c --output Figure1/Figure1c/Figure1c.pdf --overwrite
python Figure1/Figure1d/plot_Figure1d.py --rscript Rscript --overwrite
```

Source-data workbooks are supplied as archived numerical outputs; a separate workbook-building wrapper is not included. Panel 1d imports HAC functions from `Figure1c/prepare_grid_EDH30_trends.py`; 1c shares upstream functions with `Figure1a/03_prepare_mean_EDH30.py`. The 1c map reads `Figure1a/world_boundaries.gpkg`.

Crosswalks support external-input reconstruction, preserving the actual spatial assignments. The 1c compressed source table preserves pixel-level model results beyond the two plotting TIFFs.

See [environment](#environment) and [sources/rights](#sources-and-rights). The Chinese README retains detailed field definitions; they do not establish that upstream reconstruction was repeated for this deposit. Individual PDFs do not recreate the manually assembled final figure.

<a id="figure-1a"></a>

## Figure 1a: mean annual EDH30, 1986-2025


### Bundled files and plotting

`mean_annual_EDH30_1986_2025.tif` is the direct single-band plotting input in EPSG:4326. `world_boundaries.gpkg` contains Natural Earth countries and is also used by 1c. `plot_Figure1a.R` writes the reference `Figure1a.pdf`.

EDH30 is `sum(max(hourly_temperature_C - 30, 0))`, in degrees C hours. Annual totals use UTC calendar years. The map averages available maize-masked years; valid zeros remain distinct from missing cells. This is not growing-season or hot-day-normalized exposure. Positive harvested area defines annual masks; 2021-2025 reuse 2020 area. No additional area weighting is applied to the mean map.

From the repository root:

```text
Rscript Figure1/Figure1a/plot_Figure1a.R --overwrite
```

### Optional external-input reconstruction

The three numbered Python scripts export climate data, extract area and prepare the mean. Install the optional upstream dependencies listed in [Environment](#environment); obtain `Harvested_area_1961-2020.h5` (Maize group) from the provider. Authenticate Earth Engine and supply an enabled Google Cloud project:

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

<a id="figure-1b"></a>

## Figure 1b: continental annual EDH30


### Bundled files

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
```

For external reconstruction, obtain monthly climate files and annual area tables as described for 1a, then run:

```text
python Figure1/Figure1b/prepare_continent_EDH30.py --tif-dir /path/to/monthly --area-dir /path/to/work/area --output /path/to/work/continent_annual_EDH30_1986_2025.csv
```

Aggregation rejects an existing output. To fit the new CSV, place it under its expected filename in a separate package copy. The crosswalk contains IDs/labels rather than original geometry or area values; this deposit does not grant a new blanket licence for it.

<a id="figure-1c"></a>

## Figure 1c: grid-level EDH30 trends


### Bundled files

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

<a id="figure-1d"></a>

## Figure 1d: regional summer EDH30


### Bundled files

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

Retain `Figure1c/prepare_grid_EDH30_trends.py`, which provides shared HAC functions. The entry writes two PDFs; the supplied source workbook is an archived output and is not rewritten.

For optional external reconstruction, obtain 1984-2025 monthly exports (1a exporter with `--start-year 1984`) and 1985-2020 area tables (extractor with `--start-year 1985`), then run:

```text
python Figure1/Figure1d/prepare_regional_summer_EDH30.py --tif-dir /path/to/monthly --area-dir /path/to/work/area --output /path/to/work/region_summer_1985_2025.csv
python Figure1/Figure1d/plot_Figure1d.py --input /path/to/work/region_summer_1985_2025.csv --output-dir /path/to/work/plots --rscript Rscript
```

The fixed regional assignment derives from GADM 4.1 countries and Natural Earth Chinese admin-1 units. No original GADM geometry is distributed. The derived crosswalk is not declared unrestricted; clarify applicable terms before open release. Do not silently substitute different boundaries. Original area HDF5 and climate exports remain external.

<a id="sources-and-rights"></a>

## Data sources, access and rights

Source pages checked on 2026-11-08. This record documents provenance and access, not legal clearance for an open release. The repository remains private. No open-source licence for author code or blanket licence for author outputs has been selected; rights are reserved pending an explicit choice. Review access does not itself grant unrestricted redistribution.

| Source | Use here | Verified terms/access | Treatment |
| --- | --- | --- | --- |
| [ERA5-Land / Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY) | Hourly `temperature_2m`, converted to EDH30 | Official catalogue supplies acknowledgement and licence conditions | Derived results and code included; raw hourly/monthly files excluded |
| [Natural Earth](https://www.naturalearthdata.com/about/terms-of-use/) | Country geometry; continental and Chinese admin-1 assignments | Source states map data are public domain | Geometry-only GeoPackage and attribution included |
| [GADM](https://gadm.org/license.html) | Countries in the original regional assignment | Academic/non-commercial use permitted; redistribution of GADM data requires prior permission | No original GADM boundaries included. Derived ID-to-region table retained for private research review; unrestricted redistribution of that table is not asserted or verified |
| [Area dataset DOI](https://doi.org/10.11888/Terre.tpdc.300446), [source article](https://doi.org/10.1038/s41597-023-02526-z), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10495426/) | `Harvested_area_1961-2020.h5`, `Maize` group, mask and annual weights | DOI portal inaccessible during this check. Article CC BY 4.0 does not alone establish the HDF5 redistribution licence | HDF5 and annual area tables excluded; obtain data and applicable permission from provider |
| [Zenodo 7408852](https://zenodo.org/records/7408852) | Alternate area/fertilizer dataset record | Metadata public; files explicitly restricted on check date | Not represented as an anonymous download mirror |

ERA5-Land acknowledgement: **Generated using Copernicus Climate Change Service Information (2026).** Exposure metrics are modifications/analyses by the manuscript authors. Neither the European Commission nor ECMWF is responsible for any use made of the Copernicus information or data it contains. The year denotes this deposit, not a new climate-data version.

Made with Natural Earth. Country geometry was converted to GeoPackage without simplification or reprojection. Crosswalks contain IDs and labels, not boundary polygons or harvested-area values. They preserve actual assignments rather than substituting different boundaries.

### Before public release

Confirm area-file access and terms covering redistribution of derived spatial assignments and area-weighted results. Clarify the derived regional crosswalk's status with the applicable rights holder if needed. Select an author-approved code/output licence separately from third-party terms. Arrange editor/reviewer access while private. No provider contact or permission request has been made here.

<a id="figure-2"></a>

## Figure 2 files and commands

Figure2 contains 96 files: analysis code, its inputs, fitted models needed by the plotting code, and archived results. Source Data workbooks are supplied for inspection; workbook-assembly scripts, editable layout assets, duplicate formatted regression tables and verification reports are omitted.

| Directory | Files | Contents |
| --- | ---: | --- |
| `Figure2/Figure2a/` | 14 | Comparison input, overall REML code/results, histogram, map and spatial inputs; optional maize-area preparation code. |
| `Figure2/Figure2b/` | 17 | Nitrogen-group analysis, fit/plot/simple-slope code, comparison/model inputs, transformation parameters, fitted model, coefficient/covariance/Wald/curve results, PDF and Source Data. |
| `Figure2/Figure2c/` | 15 | Environment-type × nitrogen-rate analysis with fit/plot code and corresponding inputs, model, results, PDF and Source Data. |
| `Figure2/Figure2d/` | 15 | Three-way interaction analysis with fit/plot code and corresponding inputs, model, results, PDF and Source Data. |
| `Figure2/Extended_Data_Figure/` | 8 | Three input tables, plotting code, three PDFs and Source Data. |
| `Figure2/Supplementary_Figure_heat_classification/` | 6 | Profile input, classification/diagnostic code and results, two PDFs and Source Data. |
| `Figure2/Supplementary_Table_SD_imputation/` | 4 | Observed/current SD input, imputation code, reimputed SD table and grouped-median-CV workbook. |
| `Figure2/SD_sensitivity/` | 8 | Fitting and summarization code plus six core result CSVs for Baseline, CV05, CV10 and CV20. |
| `Figure2/Weather_preprocessing/` | 9 | Download/metric code, request and environment/R1 inputs, and four derived weather tables. |

Additional R packages used by Figure 2: `Matrix`, `metafor`, `patchwork`, `scales`, `openxlsx`, `mclust`, `cluster` and `lubridate`. Install these only as needed. Weather downloading additionally needs Python `cdsapi` and the user's own CDS credentials; credentials and hourly weather archives are not included.

Commands below run from the repository root. Scripts locate deposited inputs relative to their own location. The archived results can be inspected directly. Re-running scripts can replace derived outputs; use a separate working copy if retaining the archived versions.

```text
python Figure2/Figure2a/plot_map.py
Rscript Figure2/Figure2a/plot_histogram.R
Rscript Figure2/Figure2b/plot_Figure2b.R
Rscript Figure2/Figure2c/plot_Figure2c.R
Rscript Figure2/Figure2d/plot_Figure2d.R
Rscript Figure2/Extended_Data_Figure/reproduce_extended_data.R
Rscript Figure2/Supplementary_Figure_heat_classification/reproduce_classification.R
Rscript Figure2/Supplementary_Table_SD_imputation/reproduce_SD_imputation.R
```

Optional model fitting uses the deposited comparison tables by default:

```text
Rscript Figure2/Figure2a/fit_overall_REML.R
Rscript Figure2/Figure2b/fit_Figure2b.R
Rscript Figure2/Figure2c/fit_Figure2c.R
Rscript Figure2/Figure2d/fit_Figure2d.R
Rscript Figure2/Figure2b/calculate_simple_slopes.R
Rscript Figure2/SD_sensitivity/run_sensitivity.R
```

The sensitivity runner invokes `summarize_sensitivity.py` internally and generates a run-status CSV in addition to the six deposited result tables. Re-fitting these models can take substantial time. The three `model_REML.rds` files are retained because the panel plotting scripts read them. Source Data workbooks are archived outputs and are not regenerated by an omitted workbook builder.

External weather inputs are optional for consulting the deposited results. To rebuild them, replace the placeholders below with your own directories; actual machine paths are not stored in these examples. The downloader requires CDS access and can download a large volume of data. `--check-only` performs a local check without downloading.

```text
python Figure2/Weather_preprocessing/01_download_weather.py --output-dir <weather_directory> --check-only
python Figure2/Weather_preprocessing/01_download_weather.py --output-dir <weather_directory>
Rscript Figure2/Weather_preprocessing/02_calculate_weather_metrics.R --raw-dir <weather_directory>/raw_utc --output-dir <derived_output_directory>
python Figure2/Figure2a/prepare_maize_area.py --area-h5 <maize_area_hdf5> --output <mean_area_tif>
```

The weather request table describes ERA5-Land retrieval through CDS. Consult the existing attribution and applicable terms for external climate, crop-area and boundary data before redistribution. This deposit does not grant additional rights to those datasets. No weather downloads or Figure 2 model refitting were performed when preparing this deposit.

<a id="figure-3"></a>

## Figure 3 files and commands

Figure3 contains 33 files: analysis code, inputs, fitted model and archived numerical/figure outputs. Editable AI/TIFF layouts, duplicate formatted regression tables, table-assembly scripts and an archived run-status table are omitted.

| Directory | Files | Contents |
| --- | ---: | --- |
| `Figure3/Figure3a/` | 9 | Environment and profile inputs, region polygons and world boundaries, `reproduce_Figure3a.py`, two regional summary CSVs, PDF and Source Data. |
| `Figure3/Figure3b/` | 17 | Fit/plot code, comparison/model inputs, transformation parameters, fitted RDS, country coding, regional support, coefficients/covariance/variance/model-parameter and Wald-test tables, prediction curves, PDF and Source Data. |
| `Figure3/SD_sensitivity/` | 7 | Analysis entry and six core result CSVs: coefficients, regional heat coefficients, joint tests, variance components, model parameters and curve summaries. |

Figure3a starts with 856 environments and profile inputs for 610 non-Low environments (61,000 bins). Spatial assignment retains 596 non-Low environments in the 12 regions. The two regional summaries each contain 1,200 rows. Colored curves show the mean smoothed intensity with Q25–Q75 bands; dashed curves show mean cumulative exposure. Figure3b uses a separate region-based inclusion rule: 4,541 comparisons, 846 environments and 277 references. Its 108-coefficient REML model produces 3,636 prediction rows (12 regions × 3 nitrogen-rate scenarios × 101 exposure points). The two panels therefore use different inclusion sets.

Figure3a uses Python `numpy`, `pandas`, `geopandas`, `matplotlib` and `openpyxl`. Figure3b uses R `data.table`, `Matrix`, `metafor` and `ggplot2`; sensitivities also use the shared Python summarizer's dependencies listed for Figure 2. Arial is used for figure text. Commands below run from the repository root:

```text
python Figure3/Figure3a/reproduce_Figure3a.py --output-dir <output_directory>
Rscript Figure3/Figure3b/plot_Figure3b.R
```

The Figure3b plot reads the supplied `model_REML.rds`, `model_input.csv`, regional support/tests and transformations. Keep both comparison and model inputs: they serve fitting and plotting respectively. Country-coding and covariance outputs document the fitted model. Source Data workbooks are supplied as archived outputs; Figure3a's analysis script also writes its own workbook.

Optional re-fitting and SD sensitivity analysis:

```text
Rscript Figure3/Figure3b/fit_Figure3b.R
Rscript Figure3/SD_sensitivity/run_sensitivity.R
```

The sensitivity entry calls `Figure2/SD_sensitivity/run_sensitivity.R`, which uses the shared SD inputs under `Figure2/Supplementary_Table_SD_imputation/` and invokes `Figure2/SD_sensitivity/summarize_sensitivity.py`. These files must remain in the repository; the Figure3 branches in the shared code are required. The four scenarios are Baseline, CV05, CV10 and CV20. The six deposited CSVs are consultation outputs; temporary models and old status reports are omitted. The runner creates a new status table and refuses to replace existing result tables with different contents. Re-fitting can take substantial time. Plotting/fitting can replace derived panel outputs; use a separate working copy to preserve the deposited versions.

Figure3a reuses upstream weather/profile and heat-classification inputs documented for Figure 2; its spatial files are retained locally for independent panel use and retain the same source attribution. The numerical inputs are supplied with each panel, so default fitting/plotting does not require the original local master table. No Figure3 model refitting or weather downloading was performed while preparing this deposit.

<a id="figure-4"></a>

## Figure 4 files and commands

Figure4 contains 36 files: simulation-related inputs, losslessly compressed daily simulation outputs, post-processing/plotting code, numerical results and PDFs. Duplicate PNG/SVG/TIFF exports, old layouts, verification reports, input/table assembly helpers and formatted Word tables are omitted.

| Location | Files | Contents |
| --- | ---: | --- |
| `Figure4/` | 1 | Current combined figure, `Figure4_combined.pdf`. |
| `Figure4/APSIM_inputs/` | 4 | 81-site soil/management configuration, two recovered daily weather tables, and the ERA5-Land hourly downloader. |
| `Figure4/APSIM_outputs/` | 2 | `baseline_daily.csv.gz` and `tmax_plus2_daily.csv.gz`; complete original CSV bytes restored by gzip decompression. |
| `Figure4/N2O_analysis/` | 25 | Seven analysis/plotting scripts, 14 numerical CSVs and four PDFs. |
| `Figure4/Yield_analysis/` | 4 | Yield code, site-season and regional-annual tables, and PDF. |

The two simulation output files are about 33.20 and 31.55 MiB compressed. Python readers use pandas' automatic gzip detection. Decompression preserves source-row references used by the harvest-endpoint analysis. The recovered weather tables contain 134,198 baseline and 127,152 warming site-days; they are recovered report records with incomplete date coverage, not complete original APSIM driver files.

The archived simulations cover 81 sites and 2011–2022 sowing years. The warming scenario raises daily maximum temperature by 2°C, with minimum temperature unchanged. Total nitrogen is 240 kg N ha^-1 per season: basal application or a 50:50 basal/topdressing split at End of juvenile. The irrigated management includes irrigation at End of juvenile and Flag leaf. All treatments, including the treatment labelled Rainfed, also receive establishment irrigation; its amount has not been independently documented. `xIW70` denotes the user-confirmed initial-water setting of 70% of field capacity. These management definitions were supplied by the study author; the original APSIM project, model version and cultivar parameters remain unavailable. This is an archive of simulation outputs and post-processing, not a complete runnable APSIM project.

The paired seasonal comparison contains 2,913 baseline–warming pairs; `excluded_pairs.csv` records three excluded management seasons at APSIM_012 in 2013. Main-figure cumulative bars select complete sites independently within each scenario (81 baseline and 80 warming sites), whereas paired analyses use the paired endpoint tables. Preserve these distinct inclusion rules. The archived CSVs provide numerical source data; a separate Source Data workbook is not required for using this deposit.

Python code uses `numpy`, `pandas`, `matplotlib`, `geopandas`, `shapely`, `pyproj`, `pypdf`, `pymupdf` and `Pillow`. The R panel uses `ggplot2`, `scales`, `ragg` and `svglite`, plus base R `grid`. See the scripts for any additional package requirements. Commands run from the repository root; plotting can overwrite derived files, so use a separate working copy to preserve the deposited outputs.

```text
python Figure4/N2O_analysis/plot_Figure4_combined.py
python Figure4/N2O_analysis/plot_N2O_cumulative.py
python Figure4/N2O_analysis/plot_supplementary_N2O_annual.py
Rscript Figure4/N2O_analysis/plot_warming_management_comparison.R
python Figure4/Yield_analysis/plot_supplementary_yield.py
```

The combined figure requires `Rscript` on PATH; otherwise pass `--rscript <Rscript_executable>`. Its shared Python module also reads the deposited `Figure2/Figure2a/study_regions.geojson` and `world_boundaries.gpkg`. `Extended_Data_Figure4.pdf` shows cumulative processes; `Extended_Data_N2O_increments.pdf` shows increments. These are distinct outputs.

Optional upstream post-processing and external weather downloading:

```text
python Figure4/N2O_analysis/analyze_baseline.py --hourly-dir <hourly_weather_directory>
python Figure4/N2O_analysis/analyze_warming.py
python Figure4/APSIM_inputs/ERA5_hourly/scripts/download_era5_land_hourly.py --dry-run
python Figure4/APSIM_inputs/ERA5_hourly/scripts/download_era5_land_hourly.py --output-dir <external_weather_directory>
```

The baseline binning step requires external hourly weather files named `<xSite>_2010_2025_hourly_utc.csv`; these are not deposited. The downloader uses CDS ERA5-Land timeseries and requires the user's own credentials, accepted data terms and the Python dependencies `requests`, `xarray`, `h5netcdf`, `h5py` and `cdsapi`, in addition to numpy/pandas. Actual downloading can be large. No download or APSIM simulation was run while preparing this deposit. Existing downstream CSVs/PDFs can be inspected without these external inputs.

<a id="extended-data-figure-3"></a>

## Extended Data Figure 3

`Extended_Data_Figure3/` is a compact deposit of 229 files (approximately 194.5 MiB). It contains the complete figure PDF and code; the five-band mean EDH30 raster; fixed prediction parameters; RF training inputs, folds, tuning/validation results and the final model; SHAP values and plotting inputs; 96 annual NPZ predictor matrices; 30 period-mean yield rasters; regional numerical results; and supplementary figure PDFs. Model and numerical data TIFF/NPZ files are inputs/results, rather than duplicate figure exports.

The 98 upstream climate/soil/maize-area source rasters and 576 annual yield TIFFs are stored locally outside this deposit. PNG/SVG/figure TIFF duplicates, verification-only scripts/reports, child READMEs, a workbook builder and the root duplicate validation PDF are also omitted. The canonical validation PDF is under `RF_model/figures/`. Six JSON files remain because they provide runtime configuration, feature order, method details or required integrity records. Local absolute provenance paths have been replaced by external file identifiers.

| Location | Role |
| --- | --- |
| `RF_model/` | Exact 6,040-row training table with 17 predictors, original training-source inputs, fold/tuning/validation tables, `rf_model.joblib`, training/prediction/validation code and PDFs. The model was saved with scikit-learn 1.8.0. |
| `RF_model/SHAP/` | Complete SHAP values, importance and prediction decomposition, calculation code, RCS/histogram results, plotting code and PDF. SHAP explains the fitted model on the training set, not out-of-fold performance or causal effects. |
| `prediction/` | 96 NPZ matrices, shared parameters and masks, 30 period means, regional outputs and prediction/preparation/summarization code. |
| `Extended_Data_Figure3_a-e/` | Mean annual EDH30 raster, map code, necessary configuration and PDF. |
| `Extended_Data_Figure3_f-g/` | Regional nitrogen-gain heatmap code and PDF. |
| `Supplementary_Figures/` | Historical-yield and nitrogen-gain maps, using the deposited period means. |

The historical period is 2011–2022; future periods are 2030–2050 and 2070–2090 under SSP2-4.5/SSP5-8.5. Historical ERA5-Land hourly heat calculations and future ACCESS-CM2 daily-sine estimates retain different methods. The future inputs remain provisional; file integrity does not validate climate extrapolation. Maps use 35,267 common-support cells; regional gain tables retain each package's own support. Predictions concern the archived study-region maize cells, not all global maize area.

From the repository root:

```text
python Extended_Data_Figure3/assemble_Extended_Data_Figure3.py
python Extended_Data_Figure3/Extended_Data_Figure3_f-g/plot_yield_heatmaps.py
python Extended_Data_Figure3/Supplementary_Figures/plot_supplementary.py
python Extended_Data_Figure3/RF_model/plot_validation.py
python Extended_Data_Figure3/RF_model/SHAP/scripts/plot_extended_shap.py
python Extended_Data_Figure3/prediction/summarize_outputs.py
python Extended_Data_Figure3/prediction/run_prediction.py --package historical_2011_2022
```

Plotting can overwrite derived outputs; use a separate copy to preserve this deposit. The last two commands only check deposited file integrity and prediction inputs by default. They do not train or run RF inference. Batch inference requires `--execute --output-dir <new_external_package_directory>`; future packages additionally require `--accept-provisional-climate`. Consult `--help` for package names and options. Annual outputs can be regenerated from NPZ and the model without the omitted source climate rasters.

To reconstruct period summaries from external annual outputs, use `summarize_outputs.py --annual-root <external_root>`, where files are under `<external_root>/<package>/outputs/`. To rebuild NPZ or shared soil inputs, preparation scripts require `--source-dir <external_sources_directory>`; its layout must match `prediction/manifest.json`'s `external_source_files`. These optional upstream inputs are not deposited. `export_crop_year_climate.js` requires replacing `YOUR_MAIZE_ASSET_ID` with an accessible Earth Engine asset; original asset IDs and user credentials are not supplied.

Python dependencies include numpy, pandas, scipy, scikit-learn, joblib, matplotlib, geopandas, shapely, pyproj, rasterio, statsmodels and patsy. Full SHAP recalculation additionally needs shap (archived method: 0.52.0). Arial is used for figures. Geospatial boundaries are shared with `Figure2/Figure2a/`. No RF training, SHAP recalculation, batch inference or climate export was performed while preparing this compact deposit.
