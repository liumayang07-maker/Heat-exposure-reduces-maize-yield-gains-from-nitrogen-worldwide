# Heat exposure reduces maize yield gains from nitrogen worldwide

[English README](README.md)

本仓库保存论文的数据与分析代码。目前包含 **Figure1–Figure4、Extended Data Figure3 及已归档的对应补充分析**。仓库保持私有，编辑和审稿人的访问需另行安排。

根目录仅保留两份 README 和 `Figure1/`、`Figure2/`、`Figure3/`、`Figure4/`、`Extended_Data_Figure3/` 主体目录。代码、数据与分图保持原始归档版本；本地 `manuscripts/` 原件未修改。

## 目录

- [环境与安装](#环境与安装)
- [运行步骤](#运行步骤)
- [Figure1a](#figure1a)
- [Figure1b](#figure1b)
- [Figure1c](#figure1c)
- [Figure1d](#figure1d)
- [Figure2 文件与运行命令](#figure2)
- [Figure3 文件与运行命令](#figure3)
- [Figure4 文件与运行命令](#figure4)
- [Extended Data Figure3](#extended-data-figure3)
- [数据来源与使用条件](#数据来源与使用条件)

## 环境与安装

已验证环境：Windows、Python 3.12.5、R 4.5.3。以下为实际使用的直接依赖版本，不是完整环境锁定文件。建议建立独立环境，按所用 Shell 激活后安装：

```text
python -m venv .venv
python -m pip install numpy==2.4.4 pandas==2.3.3 scipy==1.17.1 statsmodels==0.14.6 matplotlib==3.10.8 rasterio==1.5.0 geopandas==1.1.3
```

仅当需要从外部气候与面积输入重建结果时，额外安装：

```text
python -m pip install h5py==3.16.0 earthengine-api==1.7.42
```

R 中按需安装：

```r
install.packages(c("terra", "sf", "ggplot2", "data.table", "nlme", "pheatmap"))
```

已验证的 R 包版本：terra 1.9.1、sf 1.1.0、ggplot2 4.0.2、data.table 1.18.2.1、nlme 3.1.168、pheatmap 1.0.13。R 需支持 Cairo PDF，sf/terra 需具备其地理空间系统依赖。Python 图件使用 Arial；缺失字体可能改变排版。将 `Rscript` 加入 PATH，或在接受 `--rscript` 的 Python 入口中指定其完整路径。未验证全新环境安装及其他操作系统。

## 运行步骤

从仓库根目录运行。以下命令会显式替换当前副本中已有输出；若要保留参照图件，请先复制整个仓库到独立目录。无需下载外部数据即可重绘：

```text
Rscript Figure1/Figure1a/plot_Figure1a.R --overwrite
Rscript Figure1/Figure1b/plot_Figure1b.R --overwrite
python Figure1/Figure1c/plot_Figure1c.py --data-dir Figure1/Figure1c --output Figure1/Figure1c/Figure1c.pdf --overwrite
python Figure1/Figure1d/plot_Figure1d.py --rscript Rscript --overwrite
```

保持主体目录结构。1c 使用 1a 的 `world_boundaries.gpkg`；1d 使用 1c 的 HAC 函数。两份 crosswalk 用于外部输入的上游重建，不是现有分图重绘的必需数据。1c 压缩 Source Data 保留标准误、置信区间等栅格之外的模型信息。

完整外部输入重建命令见英文 README 对应面板。该过程需要自行获取面积 HDF5、认证 Earth Engine 并下载气候导出结果；不要把这些上游任务当成普通重绘命令重复提交。

## Figure1a

多年平均年度 EDH30 地图。直接输入为 `mean_annual_EDH30_1986_2025.tif` 与 `world_boundaries.gpkg`；绘图入口为 `plot_Figure1a.R`，输出 `Figure1a.pdf`。三个编号 Python 脚本分别用于气候导出、面积提取和均值计算。

#### 指标与数据范围

EDH30 为每小时超过 30°C 的温度部分之和：`sum(max(T_hour_C - 30, 0))`。年度值按 UTC 日历年累计；对各年先应用正玉米收获面积掩膜，再对像元有效年份取多年平均。它不是平均气温、生育期暴露或按热日数归一化的指标。

GeoTIFF 含 599,246 个有效像元，其中 28,070 个为零；数值范围为 0–32,135.62 °C·h。零值与缺失值不同。文件不含逐年值或每个像元的有效年份数，不能默认每个像元均有完整的 40 年。地图颜色没有再按面积加权。

## Figure1b

洲际年度 EDH30。核心输入 `continent_annual_EDH30_1986_2025.csv` 为六洲 × 40 年的 240 行记录。`plot_Figure1b.R` 重拟合六个 GLS AR(1)/REML 模型，生成 PDF 和 `GLS_models.csv`；五个 sheet 的 `Figure1b_Source_Data.xlsx` 作为归档输出随附。`maize_continent_crosswalk.csv.gz` 固定原空间归属，空白归属不自动补入任何洲。

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

## Figure1c

格点趋势地图。`EDH30_slope_per_decade.tif` 和 `EDH30_HAC_BH_adjusted_pvalue.tif` 为直接作图输入；`Figure1c_Source_Data.csv.gz` 保留完整逐像元结果。`prepare_grid_EDH30_trends.py` 用于上游分析，`plot_Figure1c.py` 用于现有数据重绘。

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

## Figure1d

区域夏季热图与趋势图。核心输入为 `region_summer_1985_2025.csv` 和 `regions.csv`；`plot_Figure1d.py` 调用 `plot_heatmap.R`，生成两份 PDF；五个 sheet 的 Source Data 工作簿作为归档输出随附，不由此入口改写。`maize_region_crosswalk.csv.gz` 固定原区域归属。

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

## 数据来源与使用条件

来源页面核查日期：2026-11-08。本节为来源与访问记录，不代表已完成公开再分发许可确认。

| 数据来源 | 本仓库用途与处理 | 当前边界 |
| --- | --- | --- |
| [ERA5-Land / Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY) | 小时气温计算 EDH30，保留派生结果和处理代码 | 不随附原始小时/月度数据；遵循提供方许可与署名要求 |
| [Natural Earth](https://www.naturalearthdata.com/about/terms-of-use/) | 国家背景及部分空间归属；随附仅含几何的 GeoPackage | 官方声明其地图数据为 public domain |
| [GADM](https://gadm.org/license.html) | 原区域归属使用其国家范围 | 不分发原始边界；对派生 ID–区域表不声明无限制再分发，公开前仍需澄清适用条件 |
| [面积数据 DOI](https://doi.org/10.11888/Terre.tpdc.300446)、[来源论文](https://doi.org/10.1038/s41597-023-02526-z) | 使用 `Harvested_area_1961-2020.h5` 的 Maize 数据组作掩膜与权重 | 原 HDF5 和逐年面积表不随附；本次无法访问 DOI 门户，论文 CC BY 4.0 不能单独证明数据文件的再分发许可 |
| [Zenodo 7408852](https://zenodo.org/records/7408852) | 同名数据记录的辅助入口 | 元数据可见，文件仍受限，不能视为匿名下载镜像 |

气候数据署名：**Generated using Copernicus Climate Change Service Information (2026).** EDH30 为作者的派生分析；European Commission 和 ECMWF 不对其使用承担责任。2026 指本次归档年份，不代表新的气候数据版本。地图背景署名：**Made with Natural Earth.**

作者尚未指定代码开源许可或全部输出的统一开放许可；私有审稿访问不等于授予无限制再分发权。公开前仍需确认面积及派生空间数据适用条款、确定作者代码/输出许可，并安排审稿访问。本次未联系数据提供方或代为申请许可。

本仓库以代码及其输入输出供查阅为目的，不包含额外的投稿工作簿组装和历史验证报告。随附 Source Data 工作簿为已归档结果。保留的数据处理、模型分析和绘图代码用于说明这些结果如何得到，不要求使用者重跑整个流程。

<a id="figure2"></a>

## Figure2 文件与运行命令

Figure2 共 96 个文件，仅保留分析代码、输入、绘图需要的拟合模型和结果。Source Data 工作簿供直接查阅；工作簿组装脚本、排版资产、重复格式化回归表及核验报告未收录。

| 目录 | 文件数 | 内容 |
| --- | ---: | --- |
| `Figure2/Figure2a/` | 14 | 比较输入、总体 REML 代码与结果、直方图、地图及空间输入；另含可选面积栅格准备代码。 |
| `Figure2/Figure2b/` | 17 | 施氮分组分析的拟合、绘图、简单斜率代码，比较与模型输入、变换参数、拟合模型、系数/协方差/Wald/曲线结果、PDF 和 Source Data。 |
| `Figure2/Figure2c/` | 15 | 环境类型与施氮量分析的拟合、绘图代码，以及对应输入、模型、结果、PDF 和 Source Data。 |
| `Figure2/Figure2d/` | 15 | 三因素交互分析的拟合、绘图代码，以及对应输入、模型、结果、PDF 和 Source Data。 |
| `Figure2/Extended_Data_Figure/` | 8 | 三份输入表、绘图代码、三份 PDF 和 Source Data。 |
| `Figure2/Supplementary_Figure_heat_classification/` | 6 | 热暴露曲线输入、分类及诊断代码、分类结果、两份 PDF 和 Source Data。 |
| `Figure2/Supplementary_Table_SD_imputation/` | 4 | 当前 SD 输入、插补代码、重插补 SD 表及分组中位 CV 工作簿。 |
| `Figure2/SD_sensitivity/` | 8 | 拟合与汇总代码，以及 Baseline、CV05、CV10、CV20 四情景的六份核心结果 CSV。 |
| `Figure2/Weather_preprocessing/` | 9 | 下载及指标计算代码、下载任务和环境/R1 输入、四份派生天气表。 |

Figure2 额外使用的 R 包为 `Matrix`、`metafor`、`patchwork`、`scales`、`openxlsx`、`mclust`、`cluster`、`lubridate`，按需安装。天气下载另需 Python `cdsapi` 和使用者自己的 CDS 凭据；凭据与原始小时天气文件未收录。

以下命令从仓库根目录执行。代码根据自身位置定位已收录输入，不依赖整理者的本地绝对路径。结果可直接查阅；重新运行可能覆盖派生输出，如需保留归档版本，请在独立副本运行。

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

按需重新拟合模型时，默认读取已收录的比较输入表：

```text
Rscript Figure2/Figure2a/fit_overall_REML.R
Rscript Figure2/Figure2b/fit_Figure2b.R
Rscript Figure2/Figure2c/fit_Figure2c.R
Rscript Figure2/Figure2d/fit_Figure2d.R
Rscript Figure2/Figure2b/calculate_simple_slopes.R
Rscript Figure2/SD_sensitivity/run_sensitivity.R
```

敏感性分析入口内部调用 `summarize_sensitivity.py`，除六份核心结果外还生成运行状态 CSV；归档未保留旧状态表。重新拟合可能较耗时。三份 `model_REML.rds` 是面板绘图代码的读取输入，因此保留。Source Data 工作簿作为现成输出供查阅，未收录的组装脚本不属于运行步骤。

查阅结果不需要原始天气文件。若要重建天气指标，将以下占位符换为自己的目录；这些示例不含整理者的实际机器路径。下载需要 CDS 授权，数据量较大；`--check-only` 仅检查本地文件，不下载。

```text
python Figure2/Weather_preprocessing/01_download_weather.py --output-dir <weather_directory> --check-only
python Figure2/Weather_preprocessing/01_download_weather.py --output-dir <weather_directory>
Rscript Figure2/Weather_preprocessing/02_calculate_weather_metrics.R --raw-dir <weather_directory>/raw_utc --output-dir <derived_output_directory>
python Figure2/Figure2a/prepare_maize_area.py --area-h5 <maize_area_hdf5> --output <mean_area_tif>
```

下载任务表对应通过 CDS 获取 ERA5-Land 数据。气候、作物面积和边界数据沿用各自来源的署名与使用条件；本仓库不赋予这些外部数据额外授权。本次推送准备未下载天气数据，也未重新拟合 Figure2 模型。

<a id="figure3"></a>

## Figure3 文件与运行命令

Figure3 共 33 个文件，仅保留分析代码、输入、拟合模型和数值/图件输出。AI、TIFF 拼版、重复格式化回归表、表格组装脚本及旧运行状态表未收录。

| 目录 | 文件数 | 内容 |
| --- | ---: | --- |
| `Figure3/Figure3a/` | 9 | 环境与热暴露曲线输入、区域多边形及世界边界、`reproduce_Figure3a.py`、两份区域汇总 CSV、PDF 和 Source Data。 |
| `Figure3/Figure3b/` | 17 | 拟合/绘图代码、比较与模型输入、变换参数、拟合 RDS、国家编码、区域支持、系数/协方差/方差分量/模型参数/Wald 检验表、预测曲线、PDF 和 Source Data。 |
| `Figure3/SD_sensitivity/` | 7 | 分析入口及六份核心 CSV：系数、区域热响应系数、联合检验、方差分量、模型参数和曲线摘要。 |

Figure3a 的环境输入为 856 个环境，曲线输入为 610 个非 Low 环境、61,000 行 bin；空间归属后，12 个区域实际纳入 596 个非 Low 环境。两份区域汇总各 1,200 行。彩色曲线表示平滑强度均值及 Q25–Q75 分布带，虚线表示平均累计暴露。Figure3b 按另一套区域纳入规则使用 4,541 对比较、846 个环境、277 篇文献；108 个系数的 REML 模型对应 3,636 行预测（12 区域 × 3 个施氮情景 × 101 个热暴露点）。两面板的纳入集合不同。

Figure3a 使用 Python `numpy`、`pandas`、`geopandas`、`matplotlib`、`openpyxl`；Figure3b 使用 R `data.table`、`Matrix`、`metafor`、`ggplot2`。敏感性分析另使用 Figure2 共用 Python 汇总代码的依赖。图中文字使用 Arial。以下命令从仓库根目录运行：

```text
python Figure3/Figure3a/reproduce_Figure3a.py --output-dir <output_directory>
Rscript Figure3/Figure3b/plot_Figure3b.R
```

Figure3b 绘图读取随附的 `model_REML.rds`、`model_input.csv`、区域支持/检验表和变换参数。比较输入与模型输入分别服务拟合和绘图，均需保留；国家编码及协方差结果用于理解模型。Source Data 工作簿可直接查阅；Figure3a 分析脚本也会生成自己的工作簿。

按需重新拟合或执行 SD 敏感性分析：

```text
Rscript Figure3/Figure3b/fit_Figure3b.R
Rscript Figure3/SD_sensitivity/run_sensitivity.R
```

敏感性分析入口调用 `Figure2/SD_sensitivity/run_sensitivity.R`，复用 `Figure2/Supplementary_Table_SD_imputation/` 的 SD 输入，并由 `Figure2/SD_sensitivity/summarize_sensitivity.py` 汇总。上述文件及共用代码中的 Figure3 分支均为必要依赖。四种情景为 Baseline、CV05、CV10、CV20。归档保留六份核心 CSV 供查阅，不保存临时模型及旧状态报告；入口运行时生成新的状态表，并拒绝用不同内容覆盖现有结果。重新拟合可能较耗时；绘图或拟合可能覆盖面板派生输出，如需保留归档版本，请在独立副本运行。

Figure3a 的上游天气曲线及热分类来源见 Figure2；区域和底图文件保留在本目录，以便独立使用，沿用相同来源署名。各面板已提供必要数值输入，默认拟合/绘图不依赖整理者的原始本地主表。本次推送准备未重新拟合 Figure3 模型，也未下载天气数据。

<a id="figure4"></a>

## Figure4 文件与运行命令

Figure4 共 36 个文件，保留模拟相关输入、无损压缩逐日模拟输出、后处理/绘图代码、数值结果和 PDF。重复 PNG/SVG/TIFF、旧排版、核验报告、输入准备/表格组装脚本及格式化 Word 表未收录。

| 位置 | 文件数 | 内容 |
| --- | ---: | --- |
| `Figure4/` | 1 | 当前完整主图 `Figure4_combined.pdf`。 |
| `Figure4/APSIM_inputs/` | 4 | 81 站点土壤/管理配置、两份回提日气象及 ERA5-Land 小时下载代码。 |
| `Figure4/APSIM_outputs/` | 2 | `baseline_daily.csv.gz`、`tmax_plus2_daily.csv.gz`，解压后完整恢复原 CSV 字节。 |
| `Figure4/N2O_analysis/` | 25 | 7 个分析/绘图脚本、14 份数值 CSV、4 份 PDF。 |
| `Figure4/Yield_analysis/` | 4 | 产量代码、逐站点季节及区域年度表、PDF。 |

两份模拟输出压缩后约为 33.20、31.55 MiB，代码通过 pandas 自动识别 gzip；解压不改变成熟收获端点使用的源行号。回提气象分别为 134,198 和 127,152 个站点日，日期覆盖不完整，不能视为完整原始 APSIM 驱动气象。

模拟覆盖 81 站点、2011–2022 播种年度。增温为日最高温增加 2°C、最低温不变。每季总施氮 240 kg N ha^-1，分别采用全部基肥或基肥/End of juvenile 追肥各 50%；灌溉管理在 End of juvenile 和 Flag leaf 灌溉。包括 Rainfed 在内的所有处理均有促出苗灌溉，具体水量尚未独立追溯。`xIW70` 为作者确认的初始土壤水分相当于田间持水量 70%。上述管理定义来自研究作者确认；原始 APSIM 工程、版本和品种参数仍缺失。本目录提供模拟输出与后处理归档，不是可完整重跑 APSIM 的工程。

配对比较有 2,913 对基线–增温季节；`excluded_pairs.csv` 记录 APSIM_012 在 2013 年排除的三种管理季节。主图累计柱按情景分别筛选完整站点，基线 81 个、增温 80 个；配对分析则使用配对端点表，两者纳入规则不同。现成 CSV 已提供数值来源，无需另建 Source Data 工作簿。

Python 使用 `numpy`、`pandas`、`matplotlib`、`geopandas`、`shapely`、`pyproj`、`pypdf`、`pymupdf`、`Pillow`；R 面板使用 `ggplot2`、`scales`、`ragg`、`svglite` 及基础 R 的 `grid`，其他依赖以脚本为准。以下命令从仓库根目录执行；绘图可能覆盖派生文件，如需保留归档版本，请在独立副本运行。

```text
python Figure4/N2O_analysis/plot_Figure4_combined.py
python Figure4/N2O_analysis/plot_N2O_cumulative.py
python Figure4/N2O_analysis/plot_supplementary_N2O_annual.py
Rscript Figure4/N2O_analysis/plot_warming_management_comparison.R
python Figure4/Yield_analysis/plot_supplementary_yield.py
```

组合图通过 PATH 定位 `Rscript`，也可用 `--rscript <Rscript_executable>` 指定。共用绘图模块读取已归档的 `Figure2/Figure2a/study_regions.geojson` 和 `world_boundaries.gpkg`。`Extended_Data_Figure4.pdf` 为累计过程，`Extended_Data_N2O_increments.pdf` 为增量过程，二者均保留。

按需重算上游后处理或下载外部小时天气：

```text
python Figure4/N2O_analysis/analyze_baseline.py --hourly-dir <hourly_weather_directory>
python Figure4/N2O_analysis/analyze_warming.py
python Figure4/APSIM_inputs/ERA5_hourly/scripts/download_era5_land_hourly.py --dry-run
python Figure4/APSIM_inputs/ERA5_hourly/scripts/download_era5_land_hourly.py --output-dir <external_weather_directory>
```

基线分箱需要外部 `<xSite>_2010_2025_hourly_utc.csv` 文件，本仓库未收录。下载使用 CDS ERA5-Land timeseries，需使用者自行配置凭据、接受数据条款，并安装 `requests`、`xarray`、`h5netcdf`、`h5py`、`cdsapi` 等依赖；数据量较大。本次归档未下载天气或运行 APSIM。查阅现成结果无需这些外部输入。

<a id="extended-data-figure3"></a>

## Extended Data Figure3

`Extended_Data_Figure3/` 为精简版归档，共 229 个文件、约 194.5 MiB。保留完整扩展图 PDF 与代码、五波段平均 EDH30 栅格、固定参数工作簿、RF 训练输入/折分/调参/验证结果及模型、SHAP 数值与绘图输入、96 份年度 NPZ、30 份时期平均产量栅格、区域结果和补充图 PDF。数据 TIFF、NPZ 属于数值输入输出，不能作为重复图片排除。

98 份上游气候/土壤/玉米面积源栅格及 576 份年度产量 TIFF 仅本地保留；重复 PNG/SVG/图件 TIFF、纯核验代码与报告、子目录 README、表格组装脚本和根目录重复验证图未收录。验证图保留 `RF_model/figures/` 版本。六份 JSON 提供运行配置、变量顺序、方法参数或必需校验信息，继续保留；其中本地绝对来源路径已改为外部文件标识。

| 位置 | 内容 |
| --- | --- |
| `RF_model/` | 精确的 6,040 行训练表及 17 个预测变量、原始训练输入、折分/调参/验证表、`rf_model.joblib`、训练/预测/验证绘图代码与 PDF。模型保存版本为 scikit-learn 1.8.0。 |
| `RF_model/SHAP/` | 完整 SHAP、重要性和预测分解、计算代码、RCS/直方图结果、绘图代码及 PDF；解释训练集上的拟合模型，不等于 OOF 或因果效应。 |
| `prediction/` | 96 份 NPZ、共享参数和掩膜、30 份时期均值、区域输出、预测/准备/汇总代码。 |
| `Extended_Data_Figure3_a-e/` | 年均 EDH30 栅格、地图代码、必要配置及 PDF。 |
| `Extended_Data_Figure3_f-g/` | 区域氮肥增产热图代码与 PDF。 |
| `Supplementary_Figures/` | 基于现成时期均值的历史产量及氮肥增产地图。 |

历史为 2011–2022 年，未来为 SSP2-4.5/SSP5-8.5 的 2030–2050、2070–2090 年。历史 ERA5-Land 小时累计与未来 ACCESS-CM2 日温正弦估算沿用不同方法；未来输入仍为暂定状态，完整性检查不代表外推可靠性验证。地图采用 35,267 个共同支持格点，区域增产表则保留各包自身支持域；这些结果不覆盖全部全球玉米面积。

从仓库根目录运行：

```text
python Extended_Data_Figure3/assemble_Extended_Data_Figure3.py
python Extended_Data_Figure3/Extended_Data_Figure3_f-g/plot_yield_heatmaps.py
python Extended_Data_Figure3/Supplementary_Figures/plot_supplementary.py
python Extended_Data_Figure3/RF_model/plot_validation.py
python Extended_Data_Figure3/RF_model/SHAP/scripts/plot_extended_shap.py
python Extended_Data_Figure3/prediction/summarize_outputs.py
python Extended_Data_Figure3/prediction/run_prediction.py --package historical_2011_2022
```

绘图可能覆盖派生输出，请在独立副本运行。最后两条命令默认仅检查已收录文件及预测输入，不训练或执行 RF 推断。批量推断需显式加 `--execute --output-dir <new_external_package_directory>`，未来情景还需 `--accept-provisional-climate`；包名及参数见 `--help`。现成 NPZ 与模型可重新生成年度预测，不依赖未上传的源气候栅格。

基于外部年度预测重建时期汇总时，使用 `summarize_outputs.py --annual-root <external_root>`，文件应位于 `<external_root>/<package>/outputs/`。重建 NPZ 或共享土壤输入时，准备脚本要求 `--source-dir <external_sources_directory>`，目录布局应对应 `prediction/manifest.json` 的 `external_source_files`；这些上游文件不在仓库内。GEE 导出代码需将 `YOUR_MAIZE_ASSET_ID` 换为自己有权限访问的资产，不包含原资产 ID 或凭据。

Python 依赖包括 numpy、pandas、scipy、scikit-learn、joblib、matplotlib、geopandas、shapely、pyproj、rasterio、statsmodels、patsy；完整 SHAP 重算另需 shap（归档方法版本 0.52.0）。图件使用 Arial，空间边界复用 `Figure2/Figure2a/`。本次整理未训练 RF、重算 SHAP、批量预测或导出气候数据。
