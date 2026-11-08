"""Reproduce the regional heat-exposure profiles in Figure3a."""
import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl.styles import Font, PatternFill


ROOT = Path(__file__).resolve().parent
OUTSIDE = "Outside_12_map_regions"
SMOOTHING_BINS = 5
FIGURE_SIZE_MM = (103, 60)
PANEL_POSITIONS = {
    "North_America_Irrigated": (.04, .60, "NAMI"),
    "North_America_Rainfed": (.225, .64, "NAMR"),
    "Europe": (.44, .79, "EUR"),
    "Northwest_China": (.65, .80, "NWC"),
    "Northeast_China": (.85, .78, "NEC"),
    "North_China_Plain": (.85, .57, "NCP"),
    "West_Asia_North_Africa": (.49, .56, "WAS"),
    "South_Asia": (.69, .31, "SAS"),
    "Southwest_China": (.68, .54, "SWC"),
    "West_Africa": (.34, .38, "WAF"),
    "East_Africa": (.51, .27, "EAF"),
    "South_America": (.19, .13, "SAM"),
}


def prepare_data():
    # 1. 读取最简输入，按原版多边形 intersects 规则重新确定区域归属。
    environment = pd.read_csv(ROOT / "environment_input.csv")
    profiles = pd.read_csv(ROOT / "profiles_input.csv")
    regions = gpd.read_file(ROOT / "study_regions.geojson").to_crs(4326)
    assert len(environment) == 856 and environment.env_ID.is_unique
    assert not environment.isna().any().any()
    assert set(environment.heat_exposure_code) == {0, 1, 2, 3}
    assert len(regions) == 12 and regions.region_id.is_unique
    assert regions.is_valid.all() and set(regions.region_id) == set(PANEL_POSITIONS)
    for i in range(len(regions)):
        for j in range(i):
            assert regions.geometry.iloc[i].intersection(regions.geometry.iloc[j]).area <= 1e-10
    points = gpd.GeoDataFrame(
        environment,
        geometry=gpd.points_from_xy(environment.lon, environment.lat),
        crs=4326,
    )
    mapping = gpd.sjoin(
        points,
        regions[["region_id", "region_name", "geometry"]],
        how="left",
        predicate="intersects",
    )
    assert len(mapping) == 856 and mapping.env_ID.is_unique
    mapping = pd.DataFrame(mapping.drop(columns=["geometry", "index_right"]))
    mapping[["region_id", "region_name"]] = mapping[["region_id", "region_name"]].fillna(OUTSIDE)
    inside = mapping[(mapping.heat_exposure_code > 0) & mapping.region_id.ne(OUTSIDE)]
    assert len(inside) == 596
    assert mapping.region_id.eq(OUTSIDE).sum() == 20
    hot_ids = set(environment.loc[environment.heat_exposure_code > 0, "env_ID"])
    assert len(hot_ids) == 610 and set(profiles.env_ID) == hot_ids
    assert not profiles.duplicated(["env_ID", "bin"]).any()
    assert profiles.groupby("env_ID").bin.apply(lambda x: set(x) == set(range(1, 101))).all()
    assert np.isfinite(profiles[["bin_hours", "EDH30_per_hour"]].to_numpy()).all()
    assert (profiles.bin_hours > 0).all() and (profiles.EDH30_per_hour >= 0).all()

    # 2. 沿用原版：先逐环境累计暴露及 5-bin 平滑，再按区域等权汇总。
    bins = profiles.merge(inside, on="env_ID", how="inner", validate="many_to_one")
    bins = bins.sort_values(["env_ID", "bin"])
    assert len(bins) == 59600
    bins["EDH30_sum"] = bins.EDH30_per_hour * bins.bin_hours
    bins["EDH30_cumulative"] = bins.groupby("env_ID").EDH30_sum.cumsum()
    bins["EDH30_per_hour_5bin"] = bins.groupby("env_ID").EDH30_per_hour.transform(
        lambda x: x.rolling(window=SMOOTHING_BINS, center=True, min_periods=1).mean()
    )
    intensity = bins.groupby(["region_id", "region_name", "bin"]).EDH30_per_hour_5bin.agg(
        environment_n="count",
        mean="mean",
        q25=lambda x: x.quantile(.25),
        q75=lambda x: x.quantile(.75),
    ).reset_index()
    cumulative = bins.groupby(["region_id", "region_name", "bin"]).EDH30_cumulative.agg(
        mean="mean",
    ).reset_index()
    assert len(intensity) == len(cumulative) == 1200
    assert cumulative.groupby("region_id")["mean"].apply(lambda x: x.diff().dropna().ge(-1e-10).all()).all()
    counts = pd.crosstab(mapping.region_id, mapping.heat_exposure_code).reindex(
        index=regions.region_id.tolist() + [OUTSIDE],
        columns=[0, 1, 2, 3],
        fill_value=0,
    )
    counts.columns = ["low_n", "early_n", "middle_n", "late_n"]
    counts["nonLow_n"] = counts[["early_n", "middle_n", "late_n"]].sum(axis=1)
    counts["all_environment_n"] = counts.low_n + counts.nonLow_n
    counts = counts.reset_index()
    assert counts.nonLow_n.sum() == 610 and counts.all_environment_n.sum() == 856
    expected_n = counts.set_index("region_id").nonLow_n
    assert intensity.environment_n.eq(intensity.region_id.map(expected_n)).all()
    return regions, mapping, bins, counts, intensity, cumulative


def plot_profiles(regions, intensity, cumulative, cumulative_ymax, output):
    # 3. 保留原版地图式布局、双纵轴和区域配色，不输出单区图或说明性页脚。
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": .6,
    })
    ymax = max(intensity.q75.max(), intensity["mean"].max()) * 1.08
    fig = plt.figure(figsize=tuple(size / 25.4 for size in FIGURE_SIZE_MM))
    world = gpd.read_file(ROOT / "world_boundaries.gpkg").to_crs(4326)
    background = fig.add_axes([0, 0, 1, 1], zorder=0)
    world.dissolve().plot(ax=background, color="#eeeeee", edgecolor="none", linewidth=0)
    background.set(xlim=(-130, 155), ylim=(-60, 85))
    background.set_aspect("auto")
    background.set_axis_off()
    for region in regions.itertuples():
        left, bottom, label = PANEL_POSITIONS[region.region_id]
        # 为双轴单位预留左侧空间，保持原版各区域的相对排列。
        ax = fig.add_axes([.04 + .96 * left, bottom, .96 * .105, .18],
                          zorder=1, facecolor="none")
        sub = intensity[intensity.region_id.eq(region.region_id)].sort_values("bin")
        cum = cumulative[cumulative.region_id.eq(region.region_id)].sort_values("bin")
        right = ax.twinx()
        right.set_zorder(ax.get_zorder() - .1)
        ax.patch.set_visible(False)
        right.plot(
            np.r_[0, cum.bin],
            np.r_[0, cum["mean"]],
            color="#4d4d4d",
            linestyle="--",
            linewidth=.8,
        )
        right.set_ylim(0, cumulative_ymax)
        right.set_yticks([0, 2000, 4000])
        right.spines[["left", "top", "bottom"]].set_visible(False)
        right.spines["right"].set_linewidth(.4)
        right.tick_params(axis="y", labelsize=4.5, length=1.5, width=.4, pad=1)
        ax.fill_between(sub.bin, sub.q25, sub.q75, color=region.outline_color, alpha=.20, linewidth=0)
        ax.plot(sub.bin, sub["mean"], color=region.outline_color, linewidth=.65)
        ax.set(xlim=(0, 100), ylim=(0, ymax), yticks=[0, 1, 2, 3])
        ax.text(.06, .96, label, transform=ax.transAxes, ha="left", va="top", fontsize=5.5, clip_on=True)
        environment_n = sub.environment_n.unique()
        assert len(environment_n) == 1
        ax.text(.06, .66, f"n={int(environment_n[0])}", transform=ax.transAxes,
                ha="left", va="top", fontsize=4.5, clip_on=True)
        ax.xaxis.set_visible(False)
        ax.spines[["right", "top", "bottom"]].set_visible(False)
        ax.spines["left"].set_linewidth(.4)
        ax.tick_params(axis="y", labelsize=4.5, length=1.5, width=.4, pad=1)
        if region.region_id == "North_America_Irrigated":
            ax.set_ylabel("EDH30 intensity\n(°C h per hour)", fontsize=4.5, labelpad=1)
            right.set_ylabel("Cumulative EDH30 (°C h)", fontsize=4.5, labelpad=0)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        bounds = ax.get_tightbbox(renderer)
        assert bounds.x0 >= 0 and bounds.y0 >= 0
        assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
    fig.savefig(output / "Figure3a.pdf", dpi=400)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    args = parser.parse_args()
    regions, mapping, bins, counts, intensity, cumulative = prepare_data()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # 原版右轴按区域累计暴露 Q75 的最大值留出 8% 边距，不另行输出未绘制的分位数。
    cumulative_ymax = max(
        bins.groupby(["region_id", "bin"]).EDH30_cumulative.quantile(.75).max(),
        cumulative["mean"].max(),
    ) * 1.08
    plot_profiles(regions, intensity, cumulative, cumulative_ymax, args.output_dir)

    # 4. 仅保存实际曲线数据及全英文 Source Data，不增加模型或审计附件。
    intensity.to_csv(args.output_dir / "region_pooled_heat_profile_summary.csv", index=False)
    cumulative.to_csv(args.output_dir / "region_pooled_cumulative_EDH30_summary.csv", index=False)
    sheets = {
        "Environment_data": mapping,
        "Profile_input": bins[["env_ID", "bin", "bin_hours", "EDH30_per_hour"]],
        "Intensity_curves": intensity,
        "Cumulative_curves": cumulative,
        "Region_counts": counts,
    }
    with pd.ExcelWriter(args.output_dir / "Figure3a_Source_Data.xlsx", engine="openpyxl") as writer:
        for name, data in sheets.items():
            assert not data.isna().any().any()
            for value in data.select_dtypes(include="object").to_numpy().flat:
                assert isinstance(value, str) and value.isascii()
            data.to_excel(writer, sheet_name=name, index=False)
            ws = writer.sheets[name]
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for row in ws:
                for cell in row:
                    cell.font = Font(name="Arial", size=10, bold=cell.row == 1)
                    if cell.row == 1:
                        cell.fill = PatternFill("solid", fgColor="E8EDF2")
            for column in ws.columns:
                width = max(len(str(cell.value)) for cell in column[:20]) + 2
                ws.column_dimensions[column[0].column_letter].width = min(55, max(12, width))
    print("Completed: 856 environments; 596 included non-Low environments; 12 regions; 1,200 points per curve table.")


if __name__ == "__main__":
    main()
