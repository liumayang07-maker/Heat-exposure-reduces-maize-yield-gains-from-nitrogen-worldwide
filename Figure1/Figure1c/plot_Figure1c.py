"""Figure1c: original OLS/HAC/BH map, trimmed from source script 18.

Uses the shared Figure1a boundary; writes only the panel PDF.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import geopandas as gpd
import numpy as np
import rasterio

ALPHA = .05
ROOT = Path(__file__).resolve().parent


def plot(data_dir, boundary, output):
    with rasterio.open(data_dir/'EDH30_slope_per_decade.tif') as src:
        slope = src.read(1,masked=True)
        transform = src.transform
        extent = [src.bounds.left,src.bounds.right,src.bounds.bottom,src.bounds.top]
        signature = (src.shape,src.transform,src.crs)
    with rasterio.open(data_dir/'EDH30_HAC_BH_adjusted_pvalue.tif') as src:
        if signature!=(src.shape,src.transform,src.crs):
            raise ValueError('P-value grid mismatch')
        p = src.read(1,masked=True).filled(np.nan)
    boundaries = [-1600, -800, -400, -200, -100, -50, 0, 50, 100, 200, 400, 800, 1600]
    colors = ['#474b9e', '#4575b4', '#659ac7', '#88b9d7', '#abd9e9', '#d5edf2',
              '#fff0a6', '#fee090', '#fdae61', '#f46d43', '#d73027', '#bd0026']
    plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42})
    fig = plt.figure(figsize=(155/25.4, 6.3), facecolor='white')
    ax = fig.add_axes([.07, .38, .92, .50])
    cmap = ListedColormap(colors, name='EDH30_discrete_trend')
    cmap.set_bad((0, 0, 0, 0))
    cmap.set_under('#313695')
    cmap.set_over('#a50026')
    world = gpd.read_file(boundary)
    world.plot(ax=ax, color='#EEEEEE', edgecolor='none', zorder=0)
    im = ax.imshow(slope,extent=extent,origin='upper',cmap=cmap,
                   norm=BoundaryNorm(boundaries, cmap.N),interpolation='nearest',zorder=1)
    rows, cols = np.where(np.isfinite(p)&(p<ALPHA)&~np.ma.getmaskarray(slope))
    lon, lat = rasterio.transform.xy(transform,rows,cols)
    ax.scatter(lon,lat,s=.015,c='black',marker='o',linewidths=0,alpha=.65,
               rasterized=True,zorder=2)
    world.boundary.plot(ax=ax, color='#555555', linewidth=.25, zorder=3)
    ax.set(xlim=(-130,155),ylim=(-60,85))
    ax.set_xticks([-120, -60, 0, 60, 120], labels=['120\u00b0W', '60\u00b0W', '0\u00b0', '60\u00b0E', '120\u00b0E'])
    ax.set_yticks([-60, -30, 0, 30, 60], labels=['60\u00b0S', '30\u00b0S', '0\u00b0', '30\u00b0N', '60\u00b0N'])
    ax.tick_params(axis='both', which='major', direction='out', length=1.5/25.4*72,
                   width=.25, colors='#404040', labelsize=7.5, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#777777')
        spine.set_linewidth(.35)
    cax = fig.add_axes([.11, .205, .84, 4/25.4/6.3])
    bar = fig.colorbar(im,cax=cax,orientation='horizontal',
                      extend='both', boundaries=boundaries, ticks=boundaries,
                      spacing='uniform', drawedges=True)
    bar.ax.tick_params(labelsize=7.2, length=0, pad=3)
    bar.dividers.set_color('white')
    bar.dividers.set_linewidth(.4)
    bar.set_label(r'EDH30 trend ($^\circ$C h per decade)')
    bar.ax.xaxis.set_label_position('top')
    bar.ax.xaxis.label.set_size(9)
    bar.ax.xaxis.labelpad = 8
    fig.text(.5, .115, 'Black dots: BH-adjusted p < 0.05 (two-sided HAC, 3 lags; all tested grids)',
             ha='center', va='center', fontsize=7.2, color='#303030')
    fig.savefig(output, dpi=300, facecolor='white')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, default=ROOT / 'Figure1c.pdf')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError('PDF exists; use --overwrite explicitly.')
    boundary = ROOT.parent / 'Figure1a' / 'world_boundaries.gpkg'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot(args.data_dir, boundary, args.output)


if __name__ == '__main__':
    main()
