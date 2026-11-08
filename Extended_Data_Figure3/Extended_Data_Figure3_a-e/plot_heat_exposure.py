"""Prepare five mean-EDH30 rasters and shared map helpers from saved RF inputs."""
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
from matplotlib.path import Path as PlotPath
from matplotlib.patches import PathPatch
import numpy as np
from pyproj import Transformer
import rasterio
from rasterio.transform import from_bounds, from_origin
from rasterio.warp import reproject, Resampling


HOME = Path(__file__).resolve().parent
PREDICTION = HOME.parent / 'prediction'
INPUT_MANIFEST = PREDICTION / 'manifest.json'
WORLD = HOME.parents[1] / 'Figure2/Figure2a/world_boundaries.gpkg'
REGIONS = HOME.parents[1] / 'Figure2/Figure2a/study_regions.geojson'
SHAPE = (720, 1440)
TRANSFORM = from_origin(-180, 90, .25, .25)
FEATURE = 'annual_EDH30_C_h'
NODATA = -9999.
ABSOLUTE_BREAKS = [100, 200, 400, 800, 1600, 3200, 6400, 12800]
ABSOLUTE_COLORS = ['#68BDC4', '#BFE3C7', '#FFF7BC', '#FEE391', '#FEC44F',
                   '#FE9929', '#F03B20', '#D7191C', '#B10026', '#67000D']
ABSOLUTE_LABELS = ['0', '(0,100]', '(100,200]', '(200,400]', '(400,800]',
                   '(800,1600]', '(1600,3200]', '(3200,6400]', '(6400,12800]', '>12800']
GROUPS = {
    'historical_2011_2022': (PREDICTION / 'historical_2011_2022', range(2011, 2023)),
    'ssp245_2030_2050': (PREDICTION / 'ssp245_2030_2050', range(2030, 2051)),
    'ssp245_2070_2090': (PREDICTION / 'ssp245_2070_2090', range(2070, 2091)),
    'ssp585_2030_2050': (PREDICTION / 'ssp585_2030_2050', range(2030, 2051)),
    'ssp585_2070_2090': (PREDICTION / 'ssp585_2070_2090', range(2070, 2091)),
}
def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def absolute_classes(values):
    # 与 Figure1a 一致：零值独立，正值区间左开右闭，最高档不截断数据。
    result = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    assert (values[finite] >= 0).all()
    result[finite & (values == 0)] = 0
    positive = finite & (values > 0)
    result[positive] = np.searchsorted(ABSOLUTE_BREAKS, values[positive], side='left') + 1
    return result


def read_means():
    # 1. 只提取实际 RF 输入的年度 EDH30；逐年核验像元和区域身份。
    means = {}
    provenance = {}
    manifest = json.loads(INPUT_MANIFEST.read_text(encoding='utf-8'))
    records = {record['path']: record['sha256'] for record in manifest['files']}
    packages = {record['name']: record for record in manifest['packages']}
    for name, (package, years) in GROUPS.items():
        assert packages[name]['years'] == list(years)
        provenance[name] = {'years': list(years), 'inputs_checked': len(years)}
        total = None
        for year in years:
            path = package / f'data/inputs_{year}_median_P50.npz'
            expected = records[path.relative_to(PREDICTION).as_posix()]
            signature = digest(path)
            assert signature == expected, path
            with np.load(path, allow_pickle=False) as data:
                features = data['features'].tolist()
                assert features.count(FEATURE) == 1
                ids = data['pixel_id'].astype('int64')
                regions = data['region_code']
                values = data['X'][:, features.index(FEATURE)].astype('float64')
                assert len(ids) == len(np.unique(ids)) == len(values)
                assert ((ids >= 0) & (ids < np.prod(SHAPE))).all()
                assert np.isfinite(values).all() and (values >= 0).all()
                if total is None:
                    reference_ids = ids.copy()
                    reference_regions = regions.copy()
                    total = np.zeros(len(ids), dtype='float64')
                np.testing.assert_array_equal(ids, reference_ids)
                np.testing.assert_array_equal(regions, reference_regions)
                total += values
        means[name] = (reference_ids, reference_regions, total / len(years))
        print(f'{name}: {len(years)} years, {len(reference_ids)} cells', flush=True)
    # 2. 五组数据均限制到共同有效域，避免不同图之间混入空间覆盖变化。
    common = next(iter(means.values()))[0]
    for ids, _, _ in means.values():
        common = np.intersect1d(common, ids)
    assert len(common) > 0
    arrays = {}
    support = {}
    common_regions = None
    for name, (ids, regions, values) in means.items():
        order = np.argsort(ids)
        index = order[np.searchsorted(ids[order], common)]
        np.testing.assert_array_equal(ids[index], common)
        if common_regions is None:
            common_regions = regions[index]
        np.testing.assert_array_equal(regions[index], common_regions)
        layer = np.full(SHAPE, np.nan, dtype='float64')
        layer.ravel()[common] = values[index]
        arrays[name] = layer
        support[name] = {'original_cells': len(ids), 'common_cells': len(common), 'excluded_cells': len(ids) - len(common)}
    print(f'All-map common support: {len(common)} cells', flush=True)
    return arrays, support, provenance


def map_geometry():
    crs = '+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs'
    world = gpd.read_file(WORLD).to_crs(crs)
    regions = gpd.read_file(REGIONS)
    assert len(regions) == 12 and regions.region_id.is_unique
    assert regions.outline_color.notna().all()
    regions = regions.to_crs(crs)
    project = Transformer.from_crs('EPSG:4326', crs, always_xy=True).transform
    # 仅裁切显示边界，保留原始全球栅格和统计范围。
    lon = np.r_[np.linspace(-150, 150, 361), np.full(121, 150), np.linspace(150, -150, 361), np.full(121, -150)]
    lat = np.r_[np.full(361, -60), np.linspace(-60, 75, 121), np.full(361, 75), np.linspace(75, -60, 121)]
    x, y = project(lon, lat)
    boundary = PlotPath(np.column_stack([x, y]), closed=True)
    extent = (min(x), max(x), min(y), max(y))
    target = from_bounds(extent[0], extent[2], extent[1], extent[3], 1800, 760)
    return crs, world, regions, project, boundary, extent, target


def draw_map(ax, layer, label, caption, cmap, norm, geometry):
    crs, world, regions, project, boundary, extent, target = geometry
    clip = PathPatch(boundary, transform=ax.transData, facecolor='white', edgecolor='none')
    ax.add_patch(clip)
    world.plot(ax=ax, facecolor='#eeeeee', edgecolor='#aaaaaa', linewidth=.25, zorder=1)
    dest = np.full((760, 1800), np.nan, dtype='float64')
    reproject(
        source=np.where(np.isfinite(layer), layer, NODATA),
        destination=dest,
        src_transform=TRANSFORM,
        src_crs='EPSG:4326',
        src_nodata=NODATA,
        dst_transform=target,
        dst_crs=crs,
        dst_nodata=np.nan,
        resampling=Resampling.nearest,
    )
    image = ax.imshow(dest, extent=extent, origin='upper', cmap=cmap, norm=norm, interpolation='nearest', zorder=2)
    for meridian in (-120, -60, 0, 60, 120):
        ax.plot(*project(np.full(100, meridian), np.linspace(-60, 75, 100)), color='#bbbbbb', lw=.3, ls=':', zorder=3)
    for parallel in (-30, 0, 30, 60):
        ax.plot(*project(np.linspace(-180, 180, 181), np.full(181, parallel)), color='#bbbbbb', lw=.3, ls=':', zorder=3)
    # 沿用 Figure2a 各区域的原始描边色与 0.5 pt 线宽，仅叠加边界、不填充区域。
    regions.boundary.plot(
        ax=ax,
        color=regions.outline_color.tolist(),
        linewidth=.5,
        zorder=4,
    )
    for artist in list(ax.collections) + list(ax.images) + list(ax.lines):
        artist.set_clip_path(clip)
    ax.add_patch(PathPatch(boundary, transform=ax.transData, facecolor='none', edgecolor='#777777', linewidth=.5, zorder=5))
    ax.set(xlim=extent[:2], ylim=extent[2:], aspect='equal')
    ax.set_axis_off()
    ax.text(-.025, 1.035, label, transform=ax.transAxes, fontsize=12, fontweight='bold')
    ax.text(.5, 1.035, caption, transform=ax.transAxes, fontsize=9, ha='center')
    return image


def class_legend(fig, image, position, labels, title):
    # 刻度标在绝对值色块边界，零值独立保留为第一档。
    assert labels == ABSOLUTE_LABELS
    bar = fig.colorbar(
        image,
        cax=fig.add_axes(position),
        orientation='horizontal',
        drawedges=True,
        extend='neither',
    )
    ticks = [.5, *np.arange(1.5, 9.5, 1)]
    tick_labels = ['0', *map(str, ABSOLUTE_BREAKS)]
    bar.set_ticks(ticks)
    bar.set_ticklabels(tick_labels)
    bar.minorticks_off()
    bar.ax.tick_params(labelsize=7.2, length=2, width=.45, pad=3)
    bar.dividers.set_color('#999999')
    bar.dividers.set_linewidth(.35)
    bar.outline.set_edgecolor('#777777')
    bar.outline.set_linewidth(.45)
    bar.set_label(title + r' ($^{\circ}$C h)', labelpad=4, fontsize=8)
    return bar



def write_raster(arrays):
    # 3. 五个波段仅保留五组时期均值，并逐像元回读验证。
    path = HOME / 'EDH30_map_data.tif'
    with rasterio.open(
        path,
        'w',
        driver='GTiff',
        height=SHAPE[0],
        width=SHAPE[1],
        count=len(arrays),
        dtype='float64',
        crs='EPSG:4326',
        transform=TRANSFORM,
        nodata=NODATA,
        compress='deflate',
        predictor=3,
    ) as dst:
        for band, (name, layer) in enumerate(arrays.items(), 1):
            dst.write(np.where(np.isfinite(layer), layer, NODATA), band)
            dst.set_band_description(band, name)
            dst.update_tags(
                band,
                units='degree_C_hour',
                statistic='mean annual cumulative EDH30',
            )
    with rasterio.open(path) as src:
        assert src.shape == SHAPE and src.crs.to_epsg() == 4326
        assert src.transform == TRANSFORM and src.count == len(GROUPS) == 5
        assert src.descriptions == tuple(arrays)
        for band, layer in enumerate(arrays.values(), 1):
            np.testing.assert_array_equal(src.read(band, masked=True).filled(np.nan), layer)
    return path


def main():
    arrays, support, provenance = read_means()
    raster = write_raster(arrays)
    # 4. 保存输入追溯和五波段验收记录；成图由上级拼版脚本统一生成。
    report = {
        'feature': FEATURE,
        'units': 'degree_C_hour',
        'input_manifest': '../prediction/manifest.json',
        'input_manifest_sha256': digest(INPUT_MANIFEST),
        'inputs': provenance,
        'support': support,
        'domain': 'Intersection of all five RF input supports within the 12 study regions; not all global maize land',
        'bands': {str(i): name for i, name in enumerate(arrays, 1)},
        'display': {'longitude': [-150, 150], 'latitude': [-60, 75], 'projection': 'Robinson'},
        'boundaries': {
            'world': {'path': '../../Figure2/Figure2a/world_boundaries.gpkg', 'sha256': digest(WORLD)},
            'regions': {'path': '../../Figure2/Figure2a/study_regions.geojson', 'sha256': digest(REGIONS)},
        },
        'values': {name: {'min': float(np.nanmin(layer)), 'max': float(np.nanmax(layer))}
                   for name, layer in arrays.items()},
        'raster_roundtrip': 'PASS',
        'script_sha256': digest(Path(__file__)),
        'outputs_sha256': {raster.name: digest(raster)},
        'limitations': [
            'Historical hourly ERA5-Land and future daily-sine ACCESS-CM2 retain their distinct calculation methods.',
            'Single climate model; existing provisional future inputs; no new physical validation, RF inference or significance testing.',
        ],
    }
    (HOME / 'calculation.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print('PASS: 96 inputs checked; five period-mean raster bands. Run ../assemble_Extended_Data_Figure3.py to plot.', flush=True)


if __name__ == '__main__':
    main()
