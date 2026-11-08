# Heat exposure reduces maize yield gains from nitrogen worldwide

[English README](README.md)

本仓库保存论文的数据与分析代码。目前仅包含 **Figure1**，其他图尚未纳入。仓库保持私有，编辑和审稿人的访问需另行安排。

根目录仅保留两份 README 和 `Figure1/` 主体目录。代码、数据与分图保持原始归档版本；本地 `manuscripts/` 原件未修改。

## 目录

- [环境与安装](#环境与安装)
- [运行步骤](#运行步骤)
- [Figure1a](#figure1a)
- [Figure1b](#figure1b)
- [Figure1c](#figure1c)
- [Figure1d](#figure1d)
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
