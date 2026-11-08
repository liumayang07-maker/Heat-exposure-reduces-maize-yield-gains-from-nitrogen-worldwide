"""Reproduce the original harvested-area map and study locations for Figure2a."""
import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, to_rgba_array
import numpy as np
import pandas as pd
import rasterio

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "Figure2a_map.pdf")
    args = parser.parse_args()
    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    # 1. 读取面积栅格、固定区域边界与逐条比较元数据。
    with rasterio.open(ROOT / "mean_maize_area_2011_2020.tif") as src:
        raster = src.read(1)
    outlines = gpd.read_file(ROOT / "study_regions.geojson")
    meta = pd.read_csv(ROOT / "comparison_input.csv")
    env = meta[["env_ID", "lon", "lat"]].drop_duplicates()
    assert env.env_ID.is_unique and meta.data_ID.is_unique
    sites = meta.groupby(["lon", "lat"], as_index=False).agg(
        n_environments=("env_ID", "nunique"),
        n_references=("reference_ID", "nunique"),
        n_comparisons=("data_ID", "size"),
    )
    maximum = int(sites.n_environments.max())
    labels = ["[1, 5)", "[5, 10)", "[10, 15)", "[15, 20)", f"[20, {maximum}]"]
    palette = np.array(["#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"])
    diameters_mm = np.array([1.8, 2.7, 3.6, 4.5, 5.4]) * (2 / 3)
    classes = np.searchsorted([5, 10, 15, 20], sites.n_environments, side="right")
    sites["point_fill"] = palette[classes]
    sites["point_diameter_mm"] = diameters_mm[classes]
    fill_alpha = 0.8
    edge_color = "black"
    edge_width = 0.25
    # 2. 沿用原区域描边版地图的图层顺序、符号与版式。
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.linewidth': 0.6})
    fig, ax = plt.subplots(figsize=(210 / 25.4, 210 / 25.4 * 6.3 / 6.1))
    fig.subplots_adjust(left=.09, right=.98, top=.97, bottom=.36)
    world = gpd.read_file(ROOT / "world_boundaries.gpkg").to_crs(4326)
    world.plot(ax=ax, color='#eeeeee', edgecolor='none', zorder=0)
    cmap = plt.get_cmap('YlGn').copy()
    cmap.set_bad((0, 0, 0, 0))
    vmax = float(np.nanmax(raster))
    im = ax.imshow(np.ma.masked_where(~np.isfinite(raster) | (raster <= 0), raster),
                   extent=(-180, 180, -90, 90), origin='upper', interpolation='nearest',
                   cmap=cmap, norm=LogNorm(vmin=1, vmax=vmax), zorder=1)
    world.boundary.plot(ax=ax, color='#aaaaaa', linewidth=.35, zorder=2)
    outlines.boundary.plot(ax=ax, color=outlines.outline_color.tolist(), linewidth=0.5, zorder=2.5)
    # Match original five-class diameters and fills; draw smaller points last.
    plot_sites = sites.sort_values('n_environments', ascending=False)
    ax.scatter(plot_sites.lon, plot_sites.lat,
               s=(plot_sites.point_diameter_mm * 72 / 25.4)**2,
               facecolors=to_rgba_array(plot_sites.point_fill, alpha=fill_alpha),
               edgecolors=edge_color, linewidths=edge_width, zorder=3)
    ax.set(xlim=(-130, 155), ylim=(-60, 85), xlabel='', ylabel='')
    ax.set_xticks([-120, -60, 0, 60, 120], ['120 W', '60 W', '0', '60 E', '120 E'])
    ax.set_yticks([-60, -30, 0, 30, 60], ['60 S', '30 S', '0', '30 N', '60 N'])
    ax.tick_params(length=3, color='#888888')
    ax.grid(False)
    cax = fig.add_axes([.15, .32, .79, .019])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='min')
    cb.set_ticks([1, 10, 100, 1000, 10000], labels=['1', '10', '100', '1,000', '10,000'])
    cb.ax.minorticks_off()
    cb.set_label('Mean annual maize harvested area\n(ha per grid cell), 2011-2020', fontsize=8)
    handles = [ax.scatter([], [], s=(d*72/25.4)**2, facecolors=to_rgba_array(c, alpha=fill_alpha), edgecolors=edge_color, linewidths=edge_width, label=label)
               for d, c, label in zip(diameters_mm, palette, labels)]
    fig.legend(handles=handles, title='Environments per location', loc='lower right',
               bbox_to_anchor=(.94, .12), ncol=3, frameon=False, fontsize=7, title_fontsize=8,
               columnspacing=.9, labelspacing=1.2)
    fig.savefig(output)

    plt.close(fig)
    print(f"Saved {output}; {len(sites)} locations, {len(env)} environments.")


if __name__ == "__main__":
    main()
