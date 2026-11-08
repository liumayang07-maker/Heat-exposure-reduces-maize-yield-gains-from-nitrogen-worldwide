<a id="overview"></a>

# Heat exposure reduces maize yield gains from nitrogen worldwide

Data and analysis code accompanying the manuscript with this working title. This initial deposit contains **Figure 1 only**. Other analyses are not yet included. The repository is private; editor/reviewer access must be arranged separately.

Start with [Figure 1](#figure-1). See [verification](#verification) for the exact tests performed and [sources and rights](#sources-and-rights) for attribution and unresolved permissions. No publication or archive DOI has been assigned here.

[中文说明 / Chinese README](README.zh-CN.md)

## Contents

- [Environment](#environment)
- [Figure 1 overview and commands](#figure-1)
- [Figure 1a](#figure-1a)
- [Figure 1b](#figure-1b)
- [Figure 1c](#figure-1c)
- [Figure 1d](#figure-1d)
- [Data sources and rights](#sources-and-rights)
- [Verification](#verification)

<a id="environment"></a>

## Environment

Tested runtime: Python 3.12.5 and R 4.5.3 on Windows. The tested direct Python dependencies and optional upstream dependencies are listed below. Use a separate environment:

```text
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install numpy==2.4.4 pandas==2.3.3 scipy==1.17.1 statsmodels==0.14.6 matplotlib==3.10.8 rasterio==1.5.0 geopandas==1.1.3 openpyxl==3.1.5
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

1. **Bundled inputs:** regenerate five panel PDFs, six GLS model rows and two source workbooks, without downloads or the original working directory.
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
| [1d](#figure-1d) | `Figure1d/plot_Figure1d.py` | regional annual CSV and region definitions | two PDFs and source workbook |

### Run from the repository root

Archived outputs already exist. These commands explicitly replace them in your checkout. Run in a separate copy to preserve the references.

```text
Rscript Figure1/Figure1a/plot_Figure1a.R --overwrite
Rscript Figure1/Figure1b/plot_Figure1b.R --overwrite
python Figure1/build_source_data.py --rscript Rscript --overwrite
python Figure1/Figure1c/plot_Figure1c.py --data-dir Figure1/Figure1c --output Figure1/Figure1c/Figure1c.pdf --overwrite
python Figure1/Figure1d/plot_Figure1d.py --rscript Rscript --overwrite
```

`build_source_data.py` refits 1b through R and writes its source workbook. It also provides the workbook writer imported by 1d. Panel 1d imports HAC functions from `Figure1c/prepare_grid_EDH30_trends.py`; 1c shares upstream functions with `Figure1a/03_prepare_mean_EDH30.py`. The 1c map reads `Figure1a/world_boundaries.gpkg`.

Crosswalks support external-input reconstruction, preserving the actual spatial assignments. The 1c compressed source table preserves pixel-level model results beyond the two plotting TIFFs.

See [environment](#overview), [verification](#verification) and [sources/rights](#sources-and-rights). The Chinese README retains detailed field definitions; they do not establish that upstream reconstruction was repeated for this deposit. Individual PDFs do not recreate the manually assembled final figure.

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
python Figure1/build_source_data.py --rscript Rscript --overwrite
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

Retain `Figure1c/prepare_grid_EDH30_trends.py` and `Figure1/build_source_data.py`, which provide shared functions. The entry writes two PDFs and the workbook; temporary tables are cleaned automatically.

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

<a id="verification"></a>

## Figure 1 deposit verification

Verified on 2026-11-08 using Python 3.12.5 and R 4.5.3 on Windows. Exact direct dependency versions are recorded in the environment section above. Existing installed environments were used; a fresh dependency installation and other operating systems were not tested.

### Isolation and unchanged references

The 35 selected files were copied from the manuscript archive into an independent repository directory. Only documentation was revised there; scientific scripts, numerical inputs and archived outputs were not changed. A second complete Figure1 copy was used for execution, so its overwritten outputs did not replace the deposit references. The original archive is checked by SHA256 before push. The seven excluded files are the two final-layout assets and five historical audit files described in the preparation plan; no original file was deleted.

### Commands executed

The five commands in [Figure 1 commands](#figure-1) were executed with the actual Rscript path substituted. All returned exit code 0. Execution times were approximately 4.75 s (1a), 1.38 s (1b), 2.70 s (1b workbook), 4.94 s (1c) and 8.80 s (1d); these are local observations, not performance guarantees.

| Check | Result |
| --- | --- |
| 1b GLS model table | Six rows reproduce the archived table; numeric comparison rtol=1e-12, atol=1e-12 |
| 1b workbook | Five sheets reproduce: 240 annual rows, 240 plot rows, six GLS rows, six parameter rows, 24 covariance rows |
| 1d workbook | Five sheets reproduce: 615 annual rows, 615 heatmap rows, 15 trends, 15 parameter rows, 60 covariance rows |
| Workbook numerical tolerance | rtol=1e-12, atol=1e-12; sheet names/order and text values also compared |
| Five PDFs | One page each, equal page dimensions; full-page RGB rendering at 144 dpi has zero differing components against archived references |
| 1a raster | Single-band EPSG:4326; 599,246 valid, nonnegative pixels, including 28,070 zeros |
| 1c source table | 540,269 rows; pointwise intervals independently checked with t(0.975,38) and supplied HAC SE |
| 1c BH correction | Recomputed across all supplied raw P values; 308,069 adjusted P values below 0.05 |
| 1c rasters versus source table | All pixel coordinates mapped back to both rasters; valid-cell counts match. Slope tolerance rtol=1e-6, atol=1e-5; adjusted-P tolerance rtol=1e-6, atol=1e-7 |

Workbook/PDF metadata timestamps need not match byte-for-byte. Numerical and rendered output comparisons, rather than file hashes, establish reproduction.

### Limits

No Earth Engine task was submitted, no external climate or area data downloaded, and no complete upstream aggregation rerun. The 1a mean and 1c trends were not independently rebuilt from hourly/monthly inputs in this deposit check. Their archived inputs, plotting and source-table consistency were checked. The manually assembled final figure was not reconstructed. Source-access and redistribution questions remain documented in [Sources and rights](#sources-and-rights); passing numerical checks does not resolve them.

Repository presentation: only the English and Chinese READMEs and the Figure1 directory are retained at the root. Dependency commands and verification information are embedded in the READMEs; scientific files are unchanged.
