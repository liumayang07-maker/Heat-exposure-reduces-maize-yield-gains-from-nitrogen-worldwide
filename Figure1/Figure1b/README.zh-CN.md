> ????????????????????????? [????](../../VERIFICATION.md)?????????? [???????](../../RIGHTS.md)???????????????????

# Figure1b：六大洲年度 EDH30 变化趋势

本目录保留 8 个必要文件，用于从 GEE 月度栅格和年度面积表生成洲际年序列、重新拟合模型并绘制 Figure1b。Source Data 工作簿保存必要的拟合值和模型参数；不复制大型原始气候栅格、面积数据、移动平均图或日志。

## 文件详情

| 文件 | 内容与用途 |
| --- | --- |
| `continent_annual_EDH30_1986_2025.csv` | 核心输入：6 个洲 × 40 年，共 240 行、9 列。保留原年度均值及面积覆盖信息，删除未使用的 `moving_average_5y` 列；其余字段值与来源表完全一致。 |
| `plot_Figure1b.R` | 独立 R 脚本，从同目录 CSV 拟合六个 GLS AR(1) 模型，生成 PDF 和模型表。无跨目录运行依赖，不进行平滑或重新计算格点气候数据。 |
| `GLS_models.csv` | 六行模型结果，包含每十年斜率、标准误、95% 置信区间、原始 P 值、六洲 BH 校正 P 值、自由度及 AR(1) 参数。 |
| `prepare_continent_EDH30.py` | 从原 `plot_continent_annual_EDH30_moving_average.py` 的 `compute()` 裁剪，保留原始年度聚合，仅输出 240 行核心表。支持当前 12 波段和原版 48 波段 GEE 数据，按波段名取 EDH30。 |
| `maize_continent_crosswalk.csv.gz` | 从原空间归属表无损提取的 `maize_id`、`continent` 两列，共 840,477 行，约 1.94 MB；包含原项目 1980–2025 年出现过玉米的网格并集，按本次实际面积网格 ID 连接。 |
| `Figure1b.pdf` | 年度观测线、GLS 趋势线及点态 95% 拟合均值置信带。保留原分图尺寸（89 × 81 mm）与样式，未加入最终拼图中人工添加的图例和统计标注。 |
| `README.md` | 文件说明、字段定义、数据来源、模型方法与运行方式。 |
| `Figure1b_Source_Data.xlsx` | 投稿用多 sheet 工作簿，包含全部年度输入、实际作图值、模型结果、截距与残差参数、系数协方差。由实际绘图脚本重拟合生成，不从历史结果表补录。 |

## Source Data 与仓库完整性

工作簿直接位于本目录，也可作为独立投稿附件。仅保留以下 5 个英文 sheet，全部字段和文本为英文；不含 Readme、配色、校验值、栅格信息、汇总统计或字段说明 sheet。字段定义和方法说明留在本 README。

| Sheet | 内容 |
| --- | --- |
| `Annual_data` | 240 行、9 列，完整保留核心年度输入及面积覆盖信息。 |
| `Plot_data` | 240 行，年度观测值、GLS 拟合值、拟合均值标准误及 95% 区间。 |
| `GLS_models` | 六洲斜率、标准误、置信区间、原始及 BH 校正 P 值、自由度与 AR(1) 参数。 |
| `Model_parameters` | 六行，中心化基准年、截距、年度斜率、残差标准差、AR(1) 参数与自由度。 |
| `Coefficient_covariance` | 每洲 2 × 2 系数协方差，共 24 行，可与参数表重建拟合均值标准误。 |

从 `Figure1/` 运行：

```powershell
python build_source_data.py --rscript Rscript --overwrite
```

`--rscript` 可指定 Rscript 完整路径，首次生成可省略 `--overwrite`。依赖 Python 的 `pandas`、`openpyxl` 及本面板已有 R 依赖。脚本调用 `plot_Figure1b.R --source-data-dir=临时目录` 重拟合，只写入本目录 `Figure1b_Source_Data.xlsx`；临时 CSV 自动清理，不修改 PDF 或现有模型 CSV。生成时核对年度输入、六洲模型表和 Excel 回读数据。1a 不生成 Excel。

`gls_fit`、`gls_mean_se`、`gls_ci_low`、`gls_ci_high` 分别为拟合均值、均值标准误和点态 95% 区间上下限，单位均为 °C·h。`year_center` 为 2005.5，`intercept` 是该中心年份的拟合值；`slope_per_year` 单位为 °C·h/年，`residual_sigma` 为残差标准差，`ar1_phi` 为残差一阶相关参数。`Coefficient_covariance` 的 `term_row`、`term_col` 标明 `(Intercept)` 与 `year_c`，`covariance` 为 `vcov(fit)` 的对应元素，单位是两系数单位的乘积。

此前缺少的拟合值与完整线性模型参数现已补入工作簿。核心作图输入、分析代码、分析后数据及分图已具备；GEE 下载和面积提取代码仍与 1a 共用。大型上游数据获取、公开许可及最终人工拼图排版的限制与 1a 相同，尚不能宣称整个 Figure1 已完成公开发布准备。

2026-10-15 已核对工作簿年度输入和六洲模型表与现有 CSV 一致；由工作簿的截距、斜率和系数协方差独立重建全部 240 个拟合值、标准误及置信区间，均在数值容差内一致。临时重绘 PDF 的页面尺寸及绘图指令与原分图一致，原 PDF 和模型 CSV 未被覆盖。Excel 保留数值，不对极小 P 值作字符串截断，P 值以科学计数法显示。

## 核心数据字段

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

## 模型与图注

每个洲独立拟合 `area_weighted_mean ~ year_c`，其中 `year_c = year - mean(year)`；使用 `nlme::gls`、`corAR1(~ year)` 和 REML。每洲 40 个年度值，不使用移动平均。斜率及其标准误乘以 10，单位为 °C·h/十年。

斜率置信区间和拟合均值置信带使用 t 分布（自由度 38）。阴影是点态 95% 拟合均值区间，不是预测区间或同时置信带。六个斜率检验的双侧 P 值进行 BH 校正；原值保存在 `p_value`，校正值保存在 `p_BH_six_continents`。

最终 `Figure1.jpg` 中标为 `P` 的数字对应 BH 校正值，例如 Asia 的原始 P 约为 0.074，而图中 0.106 为校正值；正式图注需明确。当前 PDF 是用于拼图的无图例分图，洲别颜色和统计结果由脚本及模型表记录，不等同于已经完成最终拼图排版。

趋势同时受到气候变化及逐年玉米面积空间分布变化影响，不能解释为固定空间范围的纯气候趋势，也不是平均气温变化或产量效应。

## 数据来源与复现范围

- 输入来源：项目根目录下 `Figure1_v3/table/continent_annual_EDH30_moving_average.csv`。虽然原文件名含移动平均，当前模型只使用原始年均值。
- 上游聚合脚本：`Figure1_v3/scripts/plot_continent_annual_EDH30_moving_average.py`；模型与绘图代码源于 `Figure1_v3/scripts/plot_continent_annual_EDH30_GLS.R`。
- 上游气候项目：`E:/global_maize_edh_1980_2025`。计算读取 `GEE/era5land_heat_monthly_YEAR.tif`、`data/intermediate/maize_grid/global_maize_area_YEAR.csv.gz` 和 `data/intermediate/continent_crosswalk/global_maize_grid_continent_crosswalk.csv.gz`，并与该项目全球年度汇总表进行校验。
- 气候下载代码和玉米面积数据获取说明见相邻 [Figure1a README](../Figure1a/README.md)。为避免重复，本目录不再复制共享的导出脚本或大型数据。

完整运行顺序及可直接修改路径使用的命令见 [Figure1a 的上游处理流程](../Figure1a/README.md#从下载数据到两个面板)：GEE 导出及手动下载 → 原始 HDF5 面积提取 → 本目录的洲际年聚合 → GLS 拟合与成图。面积提取代码共用 `../Figure1a/02_prepare_maize_area.py`。Figure1a 的多年平均栅格不含逐年信息，不能替代本面板的年度输入。

洲别表来自原气候项目 `data/intermediate/continent_crosswalk/global_maize_grid_continent_crosswalk.csv.gz`，其两列及全部记录与来源完全一致。原表由 Natural Earth 国家边界按 `CONTINENT` 合并后进行点在多边形内匹配，并对部分未匹配点作最近洲别填补；835,278 条直接归属、4,372 条最近邻填补、827 条仍未归属。未归属记录的洲别为空，聚合时不计入任何洲，不能自动分配或当作气候零值。

当前原版归属脚本默认最近邻半径 75 km，保存表中实际填补距离最大约 9.98 km。最大观测距离不是搜索半径，不能据此断言参数不一致；本次额外检查发现，827 个未归属点在当前 75 km 设置下也不会新增归属。为固定实际采用的空间输入并减少边界处理依赖，本目录直接保留原归属表的必要两列。缺少网格 ID 会报错。该表并非新的分析结果或探索性文件。

## 运行与验证

在本目录运行：

```powershell
Rscript plot_Figure1b.R --overwrite
```

已有输出时需 `--overwrite`；首次生成时可省略。脚本只写入 `Figure1b.pdf` 和 `GLS_models.csv`。依赖 `data.table`、`ggplot2`、`nlme` 及 Cairo PDF 支持；无自动安装或联网操作。

验证环境：R 4.5.3、data.table 1.18.2.1、ggplot2 4.0.2、nlme 3.1.168。2026-10-14 已从原始 HDF5 重新提取 35 年面积，读取全部 40 年既有 GEE 月度数据，重算出与随附核心表逐值一致的 240 行年度数据；重拟合 GLS 结果也完全一致。重绘 PDF 的页面尺寸与绘图指令流一致。此前亦独立核查了面积加权均值、覆盖比例及 BH P 值。

数据读取采用一次提取 12 个 EDH30 波段，减少原 48 波段 TIFF 的重复解压；月度累计顺序和加权运算保留原版，完整数值对照已经通过。临时重算文件已清理；未提交新 GEE 任务。
