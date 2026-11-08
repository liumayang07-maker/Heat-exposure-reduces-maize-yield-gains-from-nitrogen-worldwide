> ????????????????????????? [????](../../VERIFICATION.md)?????????? [???????](../../RIGHTS.md)???????????????????

# Figure1a：全球玉米区域多年平均年度 EDH30

本目录保留 Figure1a 的核心数据、绘图脚本，以及从 GEE 月度导出和原始玉米面积 HDF5 开始的处理代码。面积提取步骤由 Figure1a、Figure1b 共用；大型原始数据与逐年中间数据不随投稿目录提供。

## 文件详情

| 文件 | 内容与用途 |
| --- | --- |
| `plot_Figure1a.R` | 独立绘图脚本，只读取同目录 GeoTIFF 和 GeoPackage，输出 `Figure1a.pdf`。采用原分图的离散色阶、地图范围及版式；不重新计算气候指标或玉米掩膜。 |
| `mean_annual_EDH30_1986_2025.tif` | 单波段 `mean_annual_EDH30`，表示 1986–2025 年玉米区域的多年平均年度 EDH30，单位 °C·h。WGS84 / EPSG:4326，1800 × 3600 像元，实际分辨率约 0.1000005°。从原多指标栅格无损提取，保留原网格与缺失值。 |
| `world_boundaries.gpkg` | Natural Earth 1:110m 国家边界，仅保留几何，用作背景和轮廓。由原 `ne_110m_admin_0_countries.shp` 转换，无简化或重投影，不参与气候值计算。 |
| `Figure1a.pdf` | 脚本生成的单页分图；热暴露图层为栅格，边界和文字为矢量。沿用源分图样式，与最终拼图中的人工色标布局及文字排版存在差异。 |
| `01_export_era5land_heat_annual_tifs.py` | 精简后的 GEE 导出脚本。虽然沿用原文件名，当前仅导出 1986–2025 年逐月 EDH30：每年一个 GeoTIFF、12 个 `EDH30_YYYY_MM` 波段，共 40 个任务；不含其他热指标或年度导出模式。 |
| `02_prepare_maize_area.py` | 从原版 `01_prepare_global_maize_area.py` 和 `edh_config.py` 裁剪；读取原始 HDF5 的 `Maize`，输出 1986–2020 年的 35 份面积表，仅含两个面板所需的 ID、经纬度和面积四列。 |
| `03_prepare_mean_EDH30.py` | 从原版 `11_plot_multiyear_mean_global_heat_maps.py` 裁剪；逐月求和、施加逐年玉米掩膜、计算 1986–2025 年有效年份均值，仅输出 EDH30 单波段 GeoTIFF。 |
| `README.md` | 本目录文件说明、数据获取来源及运行方法。 |

## Source Data 与仓库完整性

1a 的直接作图数据就是同目录 `mean_annual_EDH30_1986_2025.tif`，没有统计拟合模型。不再单独提供仅含说明、色阶、文件校验值、栅格信息或汇总统计的 Excel，也不将像元重复展开为大表。原 `Figure1a_Source_Data.xlsx` 已移除；TIFF、边界和原分图保持不变。

Source Data 仅保留直接作图数据和必要模型结果，工作簿中的 sheet 名、字段及文本均为英文；说明和字段解释留在中文 README 中。1b 的工作簿生成方式见 [Figure1b README](../Figure1b/README.md#source-data-与仓库完整性)。

目前可从随附输入重绘分图，并具备从 GEE 导出到本地聚合的代码。大型逐年气候 TIFF、面积 HDF5 和年度面积表不随附。公开前仍需核实面积数据访问权限、文件许可和稳定下载链接。若投稿要求地图数值数据，应提供现有均值 TIFF 或其持久公开链接，不以描述性 Excel 代替。

## 指标与数据范围

EDH30 为每小时超过 30°C 的温度部分之和：`sum(max(T_hour_C - 30, 0))`。年度值按 UTC 日历年累计；对各年先应用正玉米收获面积掩膜，再对像元有效年份取多年平均。它不是平均气温、生育期暴露或按热日数归一化的指标。

GeoTIFF 含 599,246 个有效像元，其中 28,070 个为零；数值范围为 0–32,135.62 °C·h。零值与缺失值不同。文件不含逐年值或每个像元的有效年份数，不能默认每个像元均有完整的 40 年。地图颜色没有再按面积加权。

## 玉米收获面积数据下载来源

- 数据集：*Global crop-specific nitrogen fertilization dataset in 1961–2020*。
- 正式数据入口：[国家青藏高原科学数据中心，DOI: 10.11888/Terre.tpdc.300446](https://doi.org/10.11888/Terre.tpdc.300446)。该 DOI 由来源论文的数据引用给出；本次未能直接打开数据门户，未验证当前下载权限。
- 来源论文：[Scientific Data，DOI: 10.1038/s41597-023-02526-z](https://doi.org/10.1038/s41597-023-02526-z)，其中 Data Records 明确说明收获面积文件。
- 所需文件：`Harvested_area_1961-2020.h5`；读取其中 `Maize` 数据组。数据覆盖 1961–2020 年，全球 5 角分网格、WGS84；论文将收获面积单位定义为 ha。
- [Zenodo 同名数据记录](https://zenodo.org/records/7408852) 可辅助检索；截至 2026-10-14，该记录标注文件访问受限，不能作为已验证可匿名下载的镜像。

本分析使用 1986–2020 年玉米面积，2021–2025 年沿用 2020 年面积。原项目年度表中的后五年坐标和面积已核对与 2020 年一致。掩膜条件为 `maize_area > 0`，将符合条件的玉米网格中心映射到气候栅格；不是用面积大小对本图像元加权。

原项目将 HDF5 转为 `global_maize_area_YEAR.csv.gz`；这些表是本地派生文件，不是上述门户的原始下载文件。为控制投稿文件大小，本目录不包含 HDF5、年度面积表或其压缩包。重绘现有 PDF 不需要下载它们。

## 其他数据来源

- 气候：GEE 数据集 `ECMWF/ERA5_LAND/HOURLY` 的 `temperature_2m`；脚本将 Kelvin 转为 °C，再计算 EDH30。
- 边界：Natural Earth，`ne_110m_admin_0_countries`。
- 本地溯源：原气候项目 `E:/global_maize_edh_1980_2025`；面积文件原存于 `E:/N_applicattion_rate_maize/Harvested_area_1961-2020.h5`。这些路径仅说明来源，不是重绘依赖。

## 运行方法

在本目录运行，R 需安装 `terra`、`sf`、`ggplot2` 并支持 Cairo：

```powershell
Rscript plot_Figure1a.R --overwrite
```

`--overwrite` 表示明确替换已有 PDF；首次生成时可以省略。已验证环境为 R 4.5.3、terra 1.9.1、sf 1.1.0、ggplot2 4.0.2。

需要重新获取上游气候数据时，先安装 `earthengine-api`，完成认证并使用具有 Earth Engine 权限的 Google Cloud 项目：

```powershell
earthengine authenticate
python 01_export_era5land_heat_annual_tifs.py --project YOUR_PROJECT_ID
```

该命令提交 40 个导出任务至 Google Drive 的 `global_maize_edh_era5land_exports` 文件夹，不直接下载到本地。加 `--dry-run` 只构建任务、不提交，但仍需要认证。导出使用 EPSG:4326、11132 m 尺度，后续应保留实际导出的仿射变换，不能假定精确 0.1° 网格。

为共用 Figure1d 的跨年夏季输入，导出脚本新增可选 `--start-year 1984`（42 个任务），面积提取脚本新增 `--start-year 1985`（36 年面积表）。默认仍为气候 1986–2025、面积 1986–2020，1a–1c 的原命令与分析范围不变；详见 [Figure1d README](../Figure1d/README.md)。

## 从下载数据到两个面板

以下 PowerShell 示例从 `Figure1/` 运行。先完成上述 GEE 导出，并从 Google Drive 将 40 份文件下载至自行指定的 `monthly` 文件夹，保持文件名 `era5land_heat_monthly_1986.tif` 至 `era5land_heat_monthly_2025.tif`。Python 处理脚本按波段名称读取，兼容当前 12 波段导出和原版 48 波段导出；缺失月份、重复波段或不一致的年份网格会报错。

另外按上文来源获取 `Harvested_area_1961-2020.h5`。Python 本地处理依赖 `h5py`、`numpy`、`pandas`、`rasterio`；GEE 导出另需 `earthengine-api`。在自己的 Python 环境安装这些依赖后运行。下列三个路径必须替换为实际路径；`work` 使用新目录，避免覆盖随附结果。

```powershell
$monthly = 'D:/Figure1_inputs/monthly'
$areaH5 = 'D:/Figure1_inputs/Harvested_area_1961-2020.h5'
$work = 'D:/Figure1_rebuild'

python Figure1a/02_prepare_maize_area.py --area-h5 "$areaH5" --output-dir "$work/area"
python Figure1a/03_prepare_mean_EDH30.py --tif-dir "$monthly" --area-dir "$work/area" --output "$work/Figure1a/mean_annual_EDH30_1986_2025.tif"
python Figure1b/prepare_continent_EDH30.py --tif-dir "$monthly" --area-dir "$work/area" --output "$work/Figure1b/continent_annual_EDH30_1986_2025.csv"

Copy-Item -LiteralPath Figure1a/plot_Figure1a.R,Figure1a/world_boundaries.gpkg -Destination "$work/Figure1a/"
Copy-Item -LiteralPath Figure1b/plot_Figure1b.R -Destination "$work/Figure1b/"
Rscript "$work/Figure1a/plot_Figure1a.R"
Rscript "$work/Figure1b/plot_Figure1b.R"
```

结果为工作目录下的 `Figure1a/Figure1a.pdf`、`Figure1b/Figure1b.pdf` 和 `Figure1b/GLS_models.csv`，原投稿目录不变。处理脚本拒绝覆盖既有输出；中途失败后应选用新的工作目录重跑，而不是把不完整文件当作完成结果。不会自动生成日志、年度气候缓存、PNG 或额外 QC 表。

### 原版逻辑的保留范围

- 面积数据索引为 `[year - 1961, longitude_index, latitude_index]`，预期维度 `(60, 4320, 2160)`；坐标中心和 `maize_id` 算法保持原版规则，仅保留有限且大于零的面积。2021–2025 年读取 2020 年表，不重复输出五份相同数据。
- 1a 将面积网格中心按气候栅格实际仿射变换映射到像元，以是否存在正面积网格生成布尔掩膜。月度累计保留原版 float32 运算，跨年累加使用 float64，最终均值写为 float32。缺失月份不补零，有效零值保留。
- 1b 保留原版 float64 年累计与面积加权顺序，不能直接用 1a 的 float32 年累计代替。洲别使用随 Figure1b 提供的原始空间归属表，固定实际使用的空间输入。
- 删除了与两个面板无关的 HH30、HH35、EDH35、五年移动平均、其他绘图、年度缓存及冗余审计输出。没有更改温度阈值、时间窗口或面积复用规则。

### 完整验证结果

2026-10-14 已在临时目录从原始 HDF5 和现有 GEE 月度文件完整运行上述本地链条：35 年面积表的全部保留字段与原表一致；40 年重算的 1a 均值栅格逐像元一致（包括缺失值）；1b 的 240 行年度表和六个 GLS 模型结果完全一致；两个重绘 PDF 的页面尺寸和绘图指令流与随附分图一致。临时数据已清理，投稿目录原有核心结果没有被覆盖。

另已验证 12 波段与原 48 波段格式产生相同的年度 EDH30、真实零值保留、缺月不补零，以及错误波段年份被拒绝。GEE 导出采用离线模拟检查了 40 个任务及每年 12 个波段，未重新运行远程 GEE 计算或提交任务。新下载的数据若被上游更新，可能无法与现有核心数据逐值一致。
