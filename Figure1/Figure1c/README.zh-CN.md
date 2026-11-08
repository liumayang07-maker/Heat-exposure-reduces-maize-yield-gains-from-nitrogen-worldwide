> ????????????????????????? [????](../../VERIFICATION.md)?????????? [???????](../../RIGHTS.md)???????????????????

# Figure1c：玉米区域年度 EDH30 的格点趋势

本目录只保留 7 个必要文件。作图对象是 1986–2025 年始终具有有效玉米掩膜和完整年度 EDH30 的气候像元；颜色为 OLS 每十年斜率，黑点为全体受检像元统一 BH 校正后 `p < 0.05`。不复制探索性配色图、PNG/SVG、日志、QC 报告或可由 P 值重新生成的显著性掩膜。

## 文件详情

| 文件 | 内容与用途 |
| --- | --- |
| `prepare_grid_EDH30_trends.py` | 从原版 17、18 号脚本裁剪，重建 OLS 斜率、HAC 标准误与双侧 P 值，再对全部受检像元做一次 BH 校正。只输出下面两份 TIFF 和一份 Source Data。 |
| `plot_Figure1c.py` | 从原版 18 号脚本裁剪，保留离散色阶、黑点大小、地图范围和版式，只输出 PDF。共用 `../Figure1a/world_boundaries.gpkg`，无联网下载或原 E 盘路径依赖。 |
| `EDH30_slope_per_decade.tif` | 直接作图输入：单波段 float32，每像元的 OLS 斜率，单位 °C·h/十年。未纳入分析的像元为 NaN。 |
| `EDH30_HAC_BH_adjusted_pvalue.tif` | 直接作图输入：同网格 float64 BH 校正 P 值；黑点由 `< 0.05` 得到，不另存重复掩膜。 |
| `Figure1c_Source_Data.csv.gz` | 英文列名的逐像元模型结果表，按全局栅格行、列顺序记录中心经纬度、斜率、HAC 标准误、95% 斜率置信区间、原始及 BH 校正 P 值。不含说明、色阶或汇总统计。压缩 CSV 避免将约 54 万行结果膨胀为大型 Excel。 |
| `Figure1c.pdf` | 原分图版式，155 mm × 6.3 inch。热图和显著性标记为栅格化图层。与最终 Figure1.jpg 的人工标题和图例位置不同，不等同于最终拼图排版。 |
| `README.md` | 文件、方法、字段、来源和复现命令；不作为 Source Data 工作表。 |

## 方法与字段

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

## 重绘现有分图

在本目录运行：

```powershell
python plot_Figure1c.py --overwrite
```

依赖 `numpy`、`rasterio`、`matplotlib`、`geopandas`，沿用 Arial 字体。首次输出可省略 `--overwrite`。绘图不读取 Source Data，也不重新拟合，仅读取两份 TIFF 和共享边界。

## 从 GEE 到分析结果

GEE 导出、手动下载及面积 HDF5 获取见 [Figure1a README](../Figure1a/README.md)。共用 `01_export_era5land_heat_annual_tifs.py`、`02_prepare_maize_area.py` 和 `03_prepare_mean_EDH30.py` 中的月累计、玉米掩膜函数；不重复复制共享代码。

从 `Figure1/` 运行，替换输入路径，使用新的输出目录：

```powershell
$monthly = 'D:/Figure1_inputs/monthly'
$areaH5 = 'D:/Figure1_inputs/Harvested_area_1961-2020.h5'
$work = 'D:/Figure1c_rebuild'

python Figure1a/02_prepare_maize_area.py --area-h5 "$areaH5" --output-dir "$work/area"
python Figure1c/prepare_grid_EDH30_trends.py --tif-dir "$monthly" --area-dir "$work/area" --output-dir "$work/results"
python Figure1c/plot_Figure1c.py --data-dir "$work/results" --output "$work/results/Figure1c.pdf"
```

若已有完整面积表，可以直接复用并跳过提取步骤。处理还需 `pandas`、`scipy`、`statsmodels`；面积提取另需 `h5py`。月度文件按 EDH30 波段名读取，兼容 12 波段导出和历史 48 波段文件。计算在内存中保存候选像元的 40 年序列，不写年度缓存。已有输出时拒绝覆盖；失败后的不完整输出不能用于作图，应使用新目录重跑。

若已有原气候项目的 40 份年度玉米掩膜 TIFF，可使用 `--annual-dir` 替代 `--tif-dir` 和 `--area-dir`。年度文件必须命名为 `gee_era5land_annual_heat_maize_area_YEAR.tif`，含 `annual_EDH30_YEAR` 波段及 `maize_mask_applied=True` 标签；该模式不把未掩膜的年度数据当作合格输入。

## 原始来源与边界

- 原项目：`E:/global_maize_edh_1980_2025`，仅作为历史溯源，不是公开仓库默认运行路径。
- 分析原版：`scripts/python/17_plot_edh30_long_term_linear_trend.py`、`18_edh30_trend_significance.py`。
- 月累计及掩膜原版：`scripts/python/11_plot_multiyear_mean_global_heat_maps.py`，相关必要函数已在 1a 中裁剪保留。
- 结果参照：`data/derived/grid_trends/EDH30_1986_2025_min40/`；图名为 `EDH30_linear_trend_1986_2025_HAC_BH_p005_discrete_split0.pdf`。
- 地图背景与国界复用 1a 的 Natural Earth 几何，不参与模型计算；未改变原图的像元、色阶或黑点抽样方式。

大型气候数据与面积 HDF5 不打包。外部面积数据的下载权限、再分发条件和正式公开仓库地址仍需在发布前落实。当前代码复现分图，不自动重建整个 Figure1 的人工拼版。

## 验证结果

2026-10-15 已在临时目录读取全部 40 年原年度玉米掩膜 TIFF，完整重算后才迁入本目录；原 E 盘结果未覆盖。

- 540,269 个完整像元；BH 校正后显著像元为 308,069 个，与原结果一致。
- 重算的斜率和 BH P 值两份 TIFF 与原版逐像元一致，包括 NaN、网格位置与 CRS。
- Source Data 的坐标、斜率、HAC 标准误、原始及 BH P 值全部核对；斜率按上文原版容差比较。全部斜率置信区间及全局 BH 校正亦独立复算通过。
- 新 PDF 与原分图页面尺寸相同，144 dpi 渲染逐像素一致，另已查看预览。最终拼图的人工排版未作为自动重绘目标。
- 1986 年从原月度 GEE 文件重新累计并应用玉米掩膜，与该年原年度 TIFF 的所有像元一致。本次没有重跑其余 39 年的月度累计；共享月度和面积处理链的此前验证见 1a README，不能将本次验证表述为重新运行了全部 GEE 下载。
- 小型测试验证了真零值保留、恒定序列 P=1、不完整 40 年记录排除及拒绝覆盖既有结果；HAC 数值另与 `statsmodels.OLS(...).fit(cov_type='HAC', use_t=True)` 独立对照。

本次环境：Python 3.12.5、numpy 2.4.4、pandas 2.3.3、rasterio 1.5.0、scipy 1.17.1、statsmodels 0.14.6、matplotlib 3.10.8、geopandas 1.1.3。没有提交 GEE 任务，也没有新增软件安装。
