> ????????????????????????? [????](../../VERIFICATION.md)?????????? [???????](../../RIGHTS.md)???????????????????

# Figure1d：区域夏季 EDH30 年际变化与趋势

本目录保留 10 个必要文件，分别复现上方热图与下方趋势图。不迁入探索性百分比变化图、区域 GLS 图、配色表、聚类缓存、日志或 QC 报告。

## 文件详情

| 文件 | 内容与用途 |
| --- | --- |
| `region_summer_1985_2025.csv` | 核心作图输入，15 区域 × 41 年，共 615 行。合并原 1985 年夏季表和 1986–2025 年区域表，保留实际季节、面积覆盖及加权均值；不保留未采用的百分比变化列。 |
| `regions.csv` | 15 区域的 ID、原英文名称、半球及原始顺序。保留聚类前顺序，删除未使用的配色列。 |
| `maize_region_crosswalk.csv.gz` | 原区域归属表的 `maize_id`、`region_id` 两列，保留全部 840,477 行，不重新分配边界附近网格。 |
| `prepare_regional_summer_EDH30.py` | 从原版 01、04 号脚本裁剪，按季节提取月度 EDH30、按当年玉米面积加权，输出唯一的 615 行核心表。 |
| `plot_heatmap.R` | 从原版 07 号脚本裁剪；区域内样本标准差标准化、欧氏距离和 Ward.D2 聚类，保留原热图布局。由下列主脚本调用，仅在临时目录输出热图和实际作图值。 |
| `plot_Figure1d.py` | 调用 R 绘制热图，沿用原版 09 号脚本的 OLS/HAC/BH 分析和趋势图版式，同时生成 Source Data。共用 1c 的 HAC 函数及上一级工作簿写入函数，不依赖原 E 盘代码。 |
| `Figure1d_heatmap.pdf` | 1985–2025 年热图，78 × 42 mm，保留原分图版式。 |
| `Figure1d_trend.pdf` | 1986–2025 年区域趋势及 95% 区间，135 × 78 mm，顺序与热图一致。 |
| `Figure1d_Source_Data.xlsx` | 五个全英文数据 sheet，仅含输入、热图值、模型结果、拟合参数和系数协方差。 |
| `README.md` | 中文文件说明、方法、字段及复现步骤，不放入 Source Data。 |

## 两个子图的区别

热图使用 **1985–2025 年，共 41 年**；每个区域用自身的均值和样本标准差（分母 `n-1`）计算 Z-score。先对 15 个区域的 41 维 Z-score 序列计算欧氏距离，再用 R `hclust(method='ward.D2')` 得到区域排列。年份保持时间顺序；不聚类年份，不显示聚类树。

下方趋势使用 **1986–2025 年，共 40 年**，不是热图 Z-score，而是原始面积加权夏季 EDH30。各区域分别拟合含截距的 OLS，采用 Newey–West/Bartlett HAC 标准误、3 阶滞后和 `40/38` 有限样本修正。双侧 t 检验与点态 95% 斜率区间均使用自由度 38。15 个斜率 P 值共同做一次 BH 校正，星号采用 `* q<0.05`、`** q<0.01`、`*** q<0.001`。区间没有进行多重比较调整。

1985 年数据仅用于热图，不再承担原探索性分析中的百分比基准。保留原稿两个子图不同的年份范围，不擅自统一，也不改用区域 GLS。

## 季节与区域定义

季节按 UTC 月份累计，以结束年份标记：北半球为当年 6–8 月；南半球为上年 12 月至当年 2 月。1985 年南半球夏季需要 **1984 年 12 月**。Brazil、Argentina、South Africa 的整个指定区域均按南半球季节处理，不按单个网格纬度重新划季节，沿用原分析设定。

每季先对每个玉米网格累计三个月 EDH30，只保留三个月均有效的记录，再计算 `sum(season_EDH30 * maize_area) / sum(valid_maize_area)`。面积使用季节结束年份，2021–2025 年沿用 2020 年。真实零暴露保留；缺失月份不补零。面积网格中心映射到实际气候栅格仿射坐标，不假定精确 0.1°。

原固定空间归属来自 GADM 4.1 国家范围及 Natural Earth admin-1 中国省级范围。中国四区域的原始省名如下，未按其他分区体系调整：

- `CHN_NE`：Heilongjiang、Liaoning、Jilin、Inner Mongol。
- `CHN_NCP`：Beijing、Tianjin、Hebei、Henan、Shandong。
- `CHN_SW`：Sichuan、Yunnan、Guizhou、Chongqing、Guangxi。
- `CHN_NW`：Shaanxi、Shanxi、Gansu、Ningxia、Qinghai、Xinjiang。

其余区域为 USA、IND、MEX、UKR、NGA、ETH、CAN、FRA、BRA、ARG、ZAF。原归属表使用网格中心与区域边界相交判断；349,414 行归属目标区域，491,063 行没有目标区域归属，没有多重归属记录。空白不等于温度缺失或零暴露，包含分析范围之外的网格及可能的边界缺口。为固定实际输入并避免复制约 295 MB 的边界文件，直接保留原两列归属，不把它替换成 1b 的洲别表。

指标是夏季累计超过 30°C 的温度量，单位 °C·h；**不是除以热日数的暴露强度**。最终拼图中 “on hot days” 的措辞不能据此解释为按热日归一化。趋势同时涉及气象变化和逐年玉米面积空间分布变化，不是固定权重的纯气候效应。

## Source Data 与字段

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

## 从核心数据重绘

在本目录运行：

```powershell
python plot_Figure1d.py `
  --rscript Rscript `
  --overwrite
```

Rscript 不在 PATH 时指定完整路径，如 `D:/R-4.5.3/bin/Rscript.exe`。Python 依赖 `numpy`、`pandas`、`scipy`、`statsmodels`、`matplotlib`、`rasterio`、`openpyxl`；R 依赖 `data.table`、`pheatmap`。脚本只在最终目录写入两份 PDF 和一份 XLSX，临时表自动清理，不修改输入 CSV。首次输出可省略 `--overwrite`；原分图中的字体分别沿用 R sans 和 Arial。

## 从 GEE 开始复现

共享数据来源与认证步骤见 [Figure1a README](../Figure1a/README.md)。从 `Figure1/` 运行：

```powershell
python Figure1a/01_export_era5land_heat_annual_tifs.py `
  --project YOUR_PROJECT_ID `
  --start-year 1984
```

这会提交 1984–2025 年共 42 个 GEE 月度导出任务，需手动从 Google Drive 下载。这里全年导出 1984 年，但本面板仅使用其 12 月；默认不加该参数仍只导出 1986–2025 年，1a–1c 行为不变。该命令是从头重建示例：已有 1986–2025 年结果时不要直接重复运行全量提交。脚本不自动去重、识别所缺年份或下载。

下载后替换路径，使用新的工作目录：

```powershell
$monthly = 'D:/Figure1_inputs/monthly'
$areaH5 = 'D:/Figure1_inputs/Harvested_area_1961-2020.h5'
$work = 'D:/Figure1d_rebuild'

python Figure1a/02_prepare_maize_area.py `
  --area-h5 "$areaH5" `
  --output-dir "$work/area" `
  --start-year 1985
python Figure1d/prepare_regional_summer_EDH30.py `
  --tif-dir "$monthly" `
  --area-dir "$work/area" `
  --output "$work/region_summer_1985_2025.csv"
python Figure1d/plot_Figure1d.py `
  --input "$work/region_summer_1985_2025.csv" `
  --output-dir "$work/plots" `
  --rscript Rscript
```

面积提取另需 `h5py`，下载另需 `earthengine-api`。面积表只输出 1985–2020 年共 36 份；后五年由聚合脚本读取 2020 年表。既有完整面积表可复用，跳过提取步骤。聚合脚本拒绝覆盖已有输出，无按年缓存或日志写入。

原版来源均在 `Figure1_v3/country_summer_EDH30/`：脚本 01、04、07、09；核心输入来自 `table/region_summer_annual_1986_2025.csv`、`region_summer_baseline_1985.csv`、`region_style.csv`、`grid_region_crosswalk.csv.gz`。大型气候与面积数据不随目录提供；外部数据访问、空间输入再分发条件仍需在公开前核实。两个 PDF 保留原分图，不自动重建最终 Figure1 的人工拼版。

## 验证结果

2026-10-15 已完成以下验证，临时重算结果均未覆盖原分析目录：

- 从原始 HDF5 重新提取 1985 年面积，795,423 个正面积网格的四个保留字段与原表逐值一致；1986–2020 年提取逻辑未改变，沿用此前验证。
- 从现有 GEE 月度栅格完整重算 1985–2025 年所有夏季，包含南半球 1984 年 12 月，共 615 行、全部字段与原输入表在数值容差内一致（`rtol=1e-12, atol=1e-10`）。
- 使用新聚合数据重新生成两张图及 Source Data，五个 sheet 全部回读并与随附工作簿核对一致；原始 615 个热图值、聚类顺序及 15 个趋势模型结果也分别与原版输出对照通过。
- 15 个区域均通过独立 statsmodels HAC 对照；参数和协方差可重建柱图斜率、标准误及区间。BH 校正后有 9 个区域满足 `p_BH < 0.05`，与原版一致。
- 两份随附 PDF 与原分图在 216 dpi 下逐像素一致，页面尺寸一致，并已进行目视检查。未宣称与人工拼图版式一致。
- 工作簿仅有五类数据 sheet，所有字段和文本均为英文，无 Excel 公式或错误单元格。
- 跨年季节、真零值保留和缺月排除的小型测试通过。GEE 默认 40 任务及新增 42 任务范围采用离线模拟验证，未联网提交或重新下载气候数据；离线检查不替代远程服务可用性验证。

本次只增加 1d 必要文件，并对共享脚本增加默认行为不变的起始年份选项和工作簿 P 值显示支持。新代码补充中文步骤注释，数据列名仍为英文；未创建 Shell 包装或改造无关历史脚本。
