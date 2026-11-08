"""Reproduce three supplementary maps from the archived RF predictions."""

from pathlib import Path
import argparse
import hashlib
import json
import string

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.path import Path as PlotPath
from matplotlib.patches import PathPatch
import numpy as np
from pyproj import Transformer
import rasterio
from rasterio.transform import from_bounds
from rasterio.warp import reproject, Resampling


# 1. 直接引用归档输入，不复制 TIFF、不重新训练或预测。
HERE = Path(__file__).resolve().parent
PREDICTION = HERE.parent / "prediction"
BOUNDARIES = HERE.parents[1] / "Figure2" / "Figure2a"
OUTPUT = HERE / "figures"
PACKAGES = [
    "historical_2011_2022",
    "ssp245_2030_2050",
    "ssp585_2030_2050",
    "ssp245_2070_2090",
    "ssp585_2070_2090",
]
SCENARIOS = [
    "Historical\n2011\u20132022",
    "SSP2-4.5\n2030\u20132050",
    "SSP5-8.5\n2030\u20132050",
    "SSP2-4.5\n2070\u20132090",
    "SSP5-8.5\n2070\u20132090",
]
COLORS = [
    "#FFF7BC", "#B5D76A", "#38A65A", "#00441B",
]
YIELD_BREAKS = [0, 2, 4, 6, 8, 10, 12, 14, 16]
GAIN_BREAKS = [0, .5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7]
UNIT = r"Mg ha$^{-1}$"
SOURCE_HASHES = {
    entry["path"]: entry["sha256"]
    for entry in json.loads((PREDICTION / "output_manifest.json").read_text(encoding="utf-8"))["files"]
}
CHECKED = set()
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": .6,
    "xtick.major.width": .6,
    "ytick.major.width": .6,
    "mathtext.default": "regular",
})


def checked_path(path):
    if path not in CHECKED:
        relative = path.relative_to(PREDICTION).as_posix()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == SOURCE_HASHES[relative], relative
        CHECKED.add(path)
    return path


def mean_raster(package, rate, water):
    period = "_".join(package.split("_")[-2:])
    path = checked_path(PREDICTION / package / "summary" / f"mean_yield_{period}_N{rate}_{water}.tif")
    with rasterio.open(path) as src:
        assert src.shape == (720, 1440) and src.crs.to_epsg() == 4326
        layer = src.read(1, masked=True).astype("float64").filled(np.nan)
        return layer, src.transform


def map_geometry():
    # 2. 复用扩展图 a–e 的投影、裁切边界、网线及区域描边。
    crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"
    world = gpd.read_file(BOUNDARIES / "world_boundaries.gpkg").to_crs(crs)
    regions = gpd.read_file(BOUNDARIES / "study_regions.geojson").to_crs(crs)
    assert len(regions) == 12 and regions.outline_color.notna().all()
    project = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform
    lon = np.r_[
        np.linspace(-150, 150, 361), np.full(121, 150),
        np.linspace(150, -150, 361), np.full(121, -150),
    ]
    lat = np.r_[
        np.full(361, -60), np.linspace(-60, 75, 121),
        np.full(361, 75), np.linspace(75, -60, 121),
    ]
    x, y = project(lon, lat)
    boundary = PlotPath(np.column_stack([x, y]), closed=True)
    extent = (min(x), max(x), min(y), max(y))
    target = from_bounds(extent[0], extent[2], extent[1], extent[3], 1800, 760)
    return crs, world, regions, project, boundary, extent, target


def draw_map(ax, layer, transform, cmap, norm, geometry):
    crs, world, regions, project, boundary, extent, target = geometry
    clip = PathPatch(boundary, transform=ax.transData, facecolor="white", edgecolor="none")
    ax.add_patch(clip)
    world.plot(ax=ax, facecolor="#eeeeee", edgecolor="#aaaaaa", linewidth=.25, zorder=1)
    dest = np.full((760, 1800), np.nan, dtype="float32")
    reproject(
        source=layer,
        destination=dest,
        src_transform=transform,
        src_crs="EPSG:4326",
        src_nodata=np.nan,
        dst_transform=target,
        dst_crs=crs,
        dst_nodata=np.nan,
        resampling=Resampling.nearest,
    )
    image = ax.imshow(
        dest,
        extent=extent,
        origin="upper",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        zorder=2,
    )
    for meridian in (-120, -60, 0, 60, 120):
        ax.plot(*project(np.full(100, meridian), np.linspace(-60, 75, 100)),
                color="#bbbbbb", lw=.3, ls=":", zorder=3)
    for parallel in (-30, 0, 30, 60):
        ax.plot(*project(np.linspace(-180, 180, 181), np.full(181, parallel)),
                color="#bbbbbb", lw=.3, ls=":", zorder=3)
    regions.boundary.plot(ax=ax, color=regions.outline_color.tolist(), linewidth=.5, zorder=4)
    for artist in [*ax.collections, *ax.images, *ax.lines]:
        artist.set_clip_path(clip)
    ax.add_patch(PathPatch(boundary, transform=ax.transData, facecolor="none",
                           edgecolor="#777777", linewidth=.5, zorder=5))
    ax.set(xlim=extent[:2], ylim=extent[2:], aspect="equal")
    ax.set_axis_off()
    assert np.isfinite(dest).sum() > 1000
    return image


def save_figure(fig, stem, expected_panels, extensions=("png", "pdf", "svg")):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = list(fig.texts)
    for ax in fig.axes:
        labels.extend(ax.texts)
        # 地图关闭坐标轴后，内部经纬刻度对象仍存在，但不属于实际输出。
        if ax.axison:
            labels.extend(ax.get_xticklabels())
            labels.extend(ax.get_yticklabels())
            labels.extend([ax.xaxis.label, ax.yaxis.label])
    for text in labels:
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0, (stem, text.get_text())
        assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, (stem, text.get_text())
    for extension in extensions:
        fig.savefig(OUTPUT / f"{stem}.{extension}", dpi=400, facecolor="white")
    plt.close(fig)
    print(f"Saved {stem}: {expected_panels} panels; {'/'.join(extensions).upper()}", flush=True)


def plot_maps(rate, geometry, continuous_preview=False):
    historical = rate is None
    assert not continuous_preview or historical
    columns = 3 if historical else 2
    rows = 2 if historical else 5
    width = 11.5 if historical else 210 / 25.4
    height = 5.3 if historical else 297 / 25.4
    fig = plt.figure(figsize=(width, height))
    left, right, gap = (.56, .12, .15) if historical else (.45, .45, .20)
    panel_width = (width - left - right - gap * (columns - 1)) / columns
    extent = geometry[5]
    panel_height = panel_width * (extent[3] - extent[2]) / (extent[1] - extent[0])
    breaks = YIELD_BREAKS if historical else GAIN_BREAKS
    # 增产图将主要色相变化集中于 0–4，数值刻度仍保持线性。
    cmap = LinearSegmentedColormap.from_list("yellow_lime_green", COLORS, N=256)
    if not historical:
        cmap = LinearSegmentedColormap.from_list(
            "gain_yellow_cyan_green",
            [(0.0, "#FFF3A6"), (2 / 7, "#36B8C5"),
             (4 / 7, "#238B45"), (1.0, "#00441B")],
            N=256,
        )
    norm = Normalize(vmin=breaks[0], vmax=breaks[-1])
    ranges = []
    for water_index, water in enumerate(("rainfed", "irrigated")):
        if historical:
            top = height - .70 - water_index * (panel_height + .33)
            fig.text(.18 / width, (top - panel_height / 2) / height,
                     water.capitalize(), rotation=90, ha="center", va="center", fontsize=11)
        else:
            x = left + water_index * (panel_width + gap)
            fig.text((x + panel_width / 2) / width, (height - .75 + .06) / height,
                     water.capitalize(), va="baseline", ha="center", fontsize=11)
        for item_index in range(3 if historical else 5):
            row = water_index if historical else item_index
            column = item_index if historical else water_index
            top = (height - .70 - row * (panel_height + .33)) if historical else (
                height - .75 - row * (panel_height + .23)
            )
            package = PACKAGES[0] if historical else PACKAGES[item_index]
            n = (0, 120, 240)[item_index] if historical else rate
            layer, transform = mean_raster(package, n, water)
            if not historical:
                control, control_transform = mean_raster(package, 0, water)
                assert transform == control_transform
                assert np.array_equal(np.isfinite(layer), np.isfinite(control))
                layer = layer - control
            low, high = float(np.nanmin(layer)), float(np.nanmax(layer))
            assert breaks[0] <= low <= high <= breaks[-1], (package, low, high)
            ranges.append((low, high))
            x = left + column * (panel_width + gap)
            ax = fig.add_axes([x / width, (top - panel_height) / height,
                               panel_width / width, panel_height / height])
            image = draw_map(ax, layer, transform, cmap, norm, geometry)
            ax.text(-.015 if historical else 0, 1.035 if historical else .98,
                    string.ascii_lowercase[row * columns + column],
                    transform=ax.transAxes, fontsize=12, fontfamily="Arial",
                    fontweight="bold" if historical else "normal", fontstyle="normal",
                    va="baseline" if historical else "top")
            if historical and row == 0:
                fig.text((x + panel_width / 2) / width, (height - .10) / height,
                         f"N{n}", va="top", ha="center", fontsize=11)
            elif not historical and column == 0:
                # 五行情景、两列水分条件；共用行标签以避免重复文字。
                fig.text(.5, (top + .06) / height,
                         SCENARIOS[item_index].replace("\n", " | "),
                         va="baseline" if item_index == 0 else "bottom", ha="center", fontsize=9.5)
    # 区域增产热图色阶为 314.5 × 6 pt；这里保留相同物理尺寸和边界刻度样式。
    bar_width = 314.5 / 72
    bar_height = 6 / 72
    bar_bottom = .60
    if not historical:
        # 色阶顶边距最后一行地图底边 0.12 英寸，两张增产图采用相同间距。
        last_panel_bottom = height - .75 - (rows - 1) * (panel_height + .23) - panel_height
        bar_bottom = last_panel_bottom - .12 - bar_height
    bar = fig.colorbar(
        image,
        cax=fig.add_axes([(.5 * width - .5 * bar_width) / width, bar_bottom / height,
                          bar_width / width, bar_height / height]),
        orientation="horizontal",
        spacing="uniform",
        drawedges=False,
        ticks=breaks,
    )
    bar.set_ticklabels([f"{value:g}" for value in breaks])
    bar.minorticks_off()
    bar.ax.tick_params(labelsize=9, length=2, width=.45, pad=3)
    bar.dividers.set_color("#999999")
    bar.dividers.set_linewidth(.35)
    bar.outline.set_edgecolor("#777777")
    bar.outline.set_linewidth(.45)
    title = "Mean predicted yield" if historical else f"Predicted yield gain (N{rate} \u2212 N0)"
    bar.set_label(f"{title} ({UNIT})", fontsize=10, labelpad=4)
    stem = "Supplementary_historical_yield" if historical else f"Supplementary_N{rate}_gain"
    if continuous_preview:
        save_figure(fig, stem + "_continuous_preview", rows * columns, extensions=("png",))
    else:
        save_figure(fig, stem, rows * columns)
    print(f"  Grid range: {min(v[0] for v in ranges):.6f} to {max(v[1] for v in ranges):.6f} Mg/ha")




def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gains-only",
        action="store_true",
        help="Export only the N120 and N240 gain PNG/PDF/SVG figures.",
    )
    parser.add_argument(
        "--historical-only",
        action="store_true",
        help="Export only the formal historical-yield PNG/PDF/SVG figure.",
    )
    parser.add_argument(
        "--historical-continuous-preview",
        action="store_true",
        help="Export only a continuous-color historical-yield PNG preview.",
    )
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    geometry = map_geometry()
    if args.historical_continuous_preview:
        plot_maps(None, geometry, continuous_preview=True)
    elif args.historical_only:
        plot_maps(None, geometry)
    elif args.gains_only:
        for rate in (120, 240):
            plot_maps(rate, geometry)
    else:
        for rate in (None, 120, 240):
            plot_maps(rate, geometry)
    print(f"Verified {len(CHECKED)} archived input checksums; no RF fitting or prediction.")


if __name__ == "__main__":
    main()
