<a id="overview"></a>

# Heat exposure reduces maize yield gains from nitrogen worldwide

Data and analysis code accompanying the manuscript with this working title. This initial deposit contains **Figure 1 only**. Other analyses are not yet included. The repository is private; editor/reviewer access must be arranged separately.

Start with [Figure 1](#figure-1). See [verification](#verification) for the exact tests performed and [sources and rights](#sources-and-rights) for attribution and unresolved permissions. No publication or archive DOI has been assigned here.

中文说明：本仓库目前仅包含 Figure1，保持私有。全部使用说明集中在本页；各面板的中文方法与字段解释可展开查看。本地原始稿件目录未修改。

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

Tested runtime: Python 3.12.5 and R 4.5.3 on Windows. Direct Python dependencies are pinned in `requirements.txt`; `requirements-upstream.txt` adds optional external-input preparation dependencies. Use a separate environment:

```text
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install -r requirements.txt
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

See [environment](#overview), [verification](#verification) and [sources/rights](#sources-and-rights). Chinese panel notes retain detailed fields and historical checks; they do not establish that upstream reconstruction was repeated for this deposit. Individual PDFs do not recreate the manually assembled final figure.

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

<details>
<summary>中文方法与字段说明</summary>

#### 指标与数据范围

EDH30 为每小时超过 30°C 的温度部分之和：`sum(max(T_hour_C - 30, 0))`。年度值按 UTC 日历年累计；对各年先应用正玉米收获面积掩膜，再对像元有效年份取多年平均。它不是平均气温、生育期暴露或按热日数归一化的指标。

GeoTIFF 含 599,246 个有效像元，其中 28,070 个为零；数值范围为 0–32,135.62 °C·h。零值与缺失值不同。文件不含逐年值或每个像元的有效年份数，不能默认每个像元均有完整的 40 年。地图颜色没有再按面积加权。

</details>

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

<details>
<summary>中文方法与字段说明</summary>

#### 核心数据字段

每行代表一个洲在一个 UTC 日历年的玉米收获面积加权平均暴露。

| 字段 | 定义 |
| --- | --- |
| `continent` | North America、South America、Europe、Asia、Africa 或 Oceania。 |
| `year` | 气候年份，1986–2025。 |
| `area_source_year` | 面积来源年份；1986–2020 使用当年，2021–2025 使用 2020 年。 |
| `total_maize_area` | 该洲当年正玉米收获面积之和，沿用源表面积数值。 |
| `valid_area` | 12 个月 EDH30 均有效的玉米网格所对应的面积之和。 |
| `valid_area_share` | `valid_area / total_maize_area`。 |
| `valid_cells` | 参与面积加权的有效玉米网格记录数；不应解读为互不重复的气候栅格数。 |
| `weighted_sum` | 有效网格的 `annual_EDH30 × maize_area` 之和。 |
| `area_weighted_mean` | `weighted_sum / valid_area`，单位 °C·h，表示面积加权的年度累计暴露。 |

年度 EDH30 是全年每小时 `max(T_C - 30, 0)` 的总和；上下半年均按同一 UTC 日历年归属，包括南半球。真零值保留，缺失月份不补零。面积仅在上游空间聚合时用于加权，拟合洲际年序列时不再进行第二次面积加权。

#### 模型与图注

每个洲独立拟合 `area_weighted_mean ~ year_c`，其中 `year_c = year - mean(year)`；使用 `nlme::gls`、`corAR1(~ year)` 和 REML。每洲 40 个年度值，不使用移动平均。斜率及其标准误乘以 10，单位为 °C·h/十年。

斜率置信区间和拟合均值置信带使用 t 分布（自由度 38）。阴影是点态 95% 拟合均值区间，不是预测区间或同时置信带。六个斜率检验的双侧 P 值进行 BH 校正；原值保存在 `p_value`，校正值保存在 `p_BH_six_continents`。

最终 `Figure1.jpg` 中标为 `P` 的数字对应 BH 校正值，例如 Asia 的原始 P 约为 0.074，而图中 0.106 为校正值；正式图注需明确。当前 PDF 是用于拼图的无图例分图，洲别颜色和统计结果由脚本及模型表记录，不等同于已经完成最终拼图排版。

趋势同时受到气候变化及逐年玉米面积空间分布变化影响，不能解释为固定空间范围的纯气候趋势，也不是平均气温变化或产量效应。

</details>

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

<details>
<summary>中文方法与字段说明</summary>

#### 方法与字段

EDH30 是全年 UTC 每小时 `max(T_C - 30, 0)` 的累计，不是生育期暴露或按热日数归一化的指标。月度累计采用原地图流程的 float32；先施加逐年正玉米收获面积掩膜。2021–2025 年沿用 2020 年面积。只有连续 40 年均有效的像元进入分析，零值保留，缺失值不补零。像元内不再按面积加权。

每个像元拟合含截距的 OLS 年度趋势，斜率乘以 10。HAC 使用 Newey–West/Bartlett 核、3 阶滞后和 `40/38` 有限样本修正，双侧 t 检验自由度为 38。恒定序列按原版规则设斜率和标准误为 0、P 值为 1，并纳入 BH 校正总体。BH 不是仅在地图显示范围或显著像元内计算。

| Source Data 字段 | 定义 |
| --- | --- |
| `longitude`, `latitude` | 原气候栅格像元中心坐标，EPSG:4326，单位度。 |
| `slope_per_decade` | 原版 HAC 函数按中心化时间计算的 float64 OLS 每十年斜率，°C·h/十年。地图 TIFF 按原版 17 号脚本的充分统计量公式计算并保存为 float32；等价公式的运算顺序和存储精度导致少量末位差异，按原版容差 `rtol=1e-6, atol=1e-5` 对照，不强行替换原值。 |
| `hac_se_per_decade` | HAC 斜率标准误，°C·h/十年。 |
| `ci_low_per_decade`, `ci_high_per_decade` | 斜率 ± `t(0.975, 38) × HAC_SE`；点态 95% 区间，不是 BH 调整后的区间，不是拟合均值或预测区间。 |
| `p_value` | 双侧 HAC t 检验原始 P 值。 |
| `p_BH_all_grids` | 全部完整 40 年受检像元共同校正后的 P 值。 |

该面板不绘制年度拟合曲线，因此不输出每像元 40 行拟合值，也不强行套用 1b 的六洲 GLS 参数表。原版没有保存截距、完整系数协方差或斜率置信区间；本次仅从原版斜率和 HAC 标准误计算上述斜率区间，不虚构其他历史模型结果。

</details>

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

<details>
<summary>中文方法与字段说明</summary>

#### 两个子图的区别

热图使用 **1985–2025 年，共 41 年**；每个区域用自身的均值和样本标准差（分母 `n-1`）计算 Z-score。先对 15 个区域的 41 维 Z-score 序列计算欧氏距离，再用 R `hclust(method='ward.D2')` 得到区域排列。年份保持时间顺序；不聚类年份，不显示聚类树。

下方趋势使用 **1986–2025 年，共 40 年**，不是热图 Z-score，而是原始面积加权夏季 EDH30。各区域分别拟合含截距的 OLS，采用 Newey–West/Bartlett HAC 标准误、3 阶滞后和 `40/38` 有限样本修正。双侧 t 检验与点态 95% 斜率区间均使用自由度 38。15 个斜率 P 值共同做一次 BH 校正，星号采用 `* q<0.05`、`** q<0.01`、`*** q<0.001`。区间没有进行多重比较调整。

1985 年数据仅用于热图，不再承担原探索性分析中的百分比基准。保留原稿两个子图不同的年份范围，不擅自统一，也不改用区域 GLS。

#### 季节与区域定义

季节按 UTC 月份累计，以结束年份标记：北半球为当年 6–8 月；南半球为上年 12 月至当年 2 月。1985 年南半球夏季需要 **1984 年 12 月**。Brazil、Argentina、South Africa 的整个指定区域均按南半球季节处理，不按单个网格纬度重新划季节，沿用原分析设定。

每季先对每个玉米网格累计三个月 EDH30，只保留三个月均有效的记录，再计算 `sum(season_EDH30 * maize_area) / sum(valid_maize_area)`。面积使用季节结束年份，2021–2025 年沿用 2020 年。真实零暴露保留；缺失月份不补零。面积网格中心映射到实际气候栅格仿射坐标，不假定精确 0.1°。

原固定空间归属来自 GADM 4.1 国家范围及 Natural Earth admin-1 中国省级范围。中国四区域的原始省名如下，未按其他分区体系调整：

- `CHN_NE`：Heilongjiang、Liaoning、Jilin、Inner Mongol。
- `CHN_NCP`：Beijing、Tianjin、Hebei、Henan、Shandong。
- `CHN_SW`：Sichuan、Yunnan、Guizhou、Chongqing、Guangxi。
- `CHN_NW`：Shaanxi、Shanxi、Gansu、Ningxia、Qinghai、Xinjiang。

其余区域为 USA、IND、MEX、UKR、NGA、ETH、CAN、FRA、BRA、ARG、ZAF。原归属表使用网格中心与区域边界相交判断；349,414 行归属目标区域，491,063 行没有目标区域归属，没有多重归属记录。空白不等于温度缺失或零暴露，包含分析范围之外的网格及可能的边界缺口。为固定实际输入并避免复制约 295 MB 的边界文件，直接保留原两列归属，不把它替换成 1b 的洲别表。

指标是夏季累计超过 30°C 的温度量，单位 °C·h；**不是除以热日数的暴露强度**。最终拼图中 “on hot days” 的措辞不能据此解释为按热日归一化。趋势同时涉及气象变化和逐年玉米面积空间分布变化，不是固定权重的纯气候效应。

#### Source Data 与字段

| Sheet | 内容 |
| --- | --- |
| `Annual_data` | 615 行完整核心输入；包含 1985 年，趋势模型只取 1986–2025 年。 |
| `Heatmap_data` | 615 行实际热图输入：区域、季节年份、原始均值、Z-score 和显示顺序。 |
| `Trend_models` | 15 行斜率、HAC 标准误、95% 区间、原始与 BH 校正 P 值及显示顺序；就是柱图和星号使用的数据。 |
| `Model_parameters` | 15 行中心化基准年、截距、年度斜率、残差标准差、样本数、自由度及 HAC 阶数。 |
| `Coefficient_covariance` | 每区域 2 × 2 HAC 系数协方差，共 60 行。 |

输入表的 `season_year` 是季节结束年；`area_source_year` 是面积年份；`start_month`、`end_month`、`season_days` 记录实际季节。`total_maize_area`、`valid_area`、`missing_temperature_area` 分别是区域总玉米面积、完整季节有效面积和温度不完整面积，单位沿用面积来源的 ha；`valid_area_share = valid_area / total_maize_area`。`weighted_sum` 为加权分子，`area_weighted_mean` 为均值（°C·h）。`maize_cell_count`、`valid_cell_count`、`zero_season_cell_count` 是玉米网格记录数，不是互不重复的气候像元数。

`slope_per_decade`、`HAC_se_per_decade`、`ci95_lower`、`ci95_upper` 单位为 °C·h/十年。`p_raw` 是原始 P，`p_BH` 是 15 区域共同校正值。

参数表的 `year_center=2005.5`；`intercept` 是该年拟合值（°C·h），`slope_per_year` 是年度斜率（°C·h/年），`residual_sigma=sqrt(SSR/38)`。协方差表的 `term_row`、`term_col` 对应 `intercept` 和 `year_c`，元素单位为两系数单位之积；此处为 HAC 协方差，不是普通独立残差 OLS 协方差。截距及完整协方差由现有原版独立校验模型导出，不伪称为历史已保存结果。该面板没有拟合曲线，因此不额外保存 600 行未用于绘图的拟合值。

</details>

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

Verified on 2026-11-08 using Python 3.12.5 and R 4.5.3 on Windows. Exact direct dependency versions are recorded in `requirements.txt`, `requirements-upstream.txt` and the repository README. Existing installed environments were used; a fresh dependency installation and other operating systems were not tested.

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

Workbook/PDF metadata timestamps need not match byte-for-byte. Numerical and rendered output comparisons, rather than file hashes, establish reproduction. Input checksums in `MANIFEST.sha256` identify the deposited reference files.

### Limits

No Earth Engine task was submitted, no external climate or area data downloaded, and no complete upstream aggregation rerun. The 1a mean and 1c trends were not independently rebuilt from hourly/monthly inputs in this deposit check. Their archived inputs, plotting and source-table consistency were checked. The manually assembled final figure was not reconstructed. Source-access and redistribution questions remain documented in [Sources and rights](#sources-and-rights); passing numerical checks does not resolve them.

Documentation consolidation: thirteen Markdown documents were combined into this README after numerical verification. Scientific code, data and archived figures are unchanged; the previous numerical checks remain applicable.
