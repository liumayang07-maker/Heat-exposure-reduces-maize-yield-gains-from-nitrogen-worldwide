"""Assemble five EDH30 period-mean maps and two regional yield-gain heatmaps."""
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import rasterio
from shapely.geometry import Polygon


HOME = Path(__file__).resolve().parent
HEAT = HOME / 'Extended_Data_Figure3_a-e'
PAGE_WIDTH = 673.8
MAP_LEFT = 19.2
MAP_RIGHT = 650.6
MAP_GAP = 14.0
ROW_GAP = 24.0
HISTORICAL_SCALE = .94
MAP_COLORBAR_WIDTH = 255.0
YIELD_COLORBAR_WIDTH = 314.5
COLORBAR_HEIGHT = 6.0
EXPORT_DPI = 300
PANELS = [
    ('historical_2011_2022', 'Historical | 2011-2022'),
    ('ssp245_2030_2050', 'SSP2-4.5 | 2030-2050'),
    ('ssp245_2070_2090', 'SSP2-4.5 | 2070-2090'),
    ('ssp585_2030_2050', 'SSP5-8.5 | 2030-2050'),
    ('ssp585_2070_2090', 'SSP5-8.5 | 2070-2090'),
]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


exposure = load_module('extended_data_figure3_exposure', HEAT / 'plot_heat_exposure.py')
heatmaps = load_module('extended_data_figure3_heatmaps', HOME / 'Extended_Data_Figure3_f-g/plot_yield_heatmaps.py')

# 标签放在邻近空白处；缩写与区域身份直接复用热图定义。
REGION_LABEL_POSITIONS = {
    'NCP': (143, 34),
    'NEC': (137, 60),
    'NWC': (78, 59),
    'SWC': (131, 9),
    'NAMI': (-134, 22),
    'NAMR': (-65, 51),
    'EUR': (-12, 61),
    'SAM': (-88, -28),
    'SAS': (85, 0),
    'WAS': (47, 15),
    'EAF': (58, -16),
    'WAF': (-18, 8),
}


def annotate_regions(ax, geometry):
    regions = geometry[2].set_index('region_id')
    project = geometry[3]
    assert set(regions.index) == set(heatmaps.ORDER)
    assert set(REGION_LABEL_POSITIONS) == set(heatmaps.LABELS)
    for region_id, label in zip(heatmaps.ORDER, heatmaps.LABELS):
        region = regions.loc[region_id].geometry
        # 多部件区域取最大连通部分中的内部点，防止引线落入海洋。
        if region.geom_type == 'MultiPolygon':
            region = max(region.geoms, key=lambda part: part.area)
        anchor = region.representative_point()
        text = ax.annotate(
            label,
            xy=(anchor.x, anchor.y),
            xytext=project(*REGION_LABEL_POSITIONS[label]),
            ha='center',
            va='center',
            fontsize=8,
            fontfamily='Arial',
            fontweight='normal',
            color='#222222',
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .8, 'pad': .5},
            arrowprops={'arrowstyle': '-', 'color': '#555555', 'lw': .45},
            zorder=10,
        )
        text.set_gid(f'region-{label}')


def map_layout(geometry):
    # 1. 按手动版布局：历史图在左上略缩小，右上留给热暴露色阶。
    extent = geometry[5]
    ratio = (extent[3] - extent[2]) / (extent[1] - extent[0])
    full_width = MAP_RIGHT - MAP_LEFT
    small_width = (full_width - MAP_GAP) / 2
    small_height = small_width * ratio
    historical_width = small_width * HISTORICAL_SCALE
    historical_height = historical_width * ratio
    top = 28.0
    second_top = top + historical_height + ROW_GAP
    third_top = second_top + small_height + ROW_GAP
    right_left = MAP_LEFT + small_width + MAP_GAP
    rectangles = [
        (MAP_LEFT + (small_width - historical_width) / 2, top,
         historical_width, historical_height),
        (MAP_LEFT, second_top, small_width, small_height),
        (right_left, second_top, small_width, small_height),
        (MAP_LEFT, third_top, small_width, small_height),
        (right_left, third_top, small_width, small_height),
    ]
    bar_top = top + historical_height * .55
    return rectangles, bar_top


def draw_maps(fig, arrays, geometry, rectangles, bar_top):
    # 2. 五幅图使用时期均值，并共享绝对值色阶。
    page_height = fig.get_figheight() * 72
    cmap = ListedColormap(exposure.ABSOLUTE_COLORS)
    norm = BoundaryNorm(np.arange(-.5, len(exposure.ABSOLUTE_COLORS) + .5, 1), cmap.N)
    axes = []
    for index, ((name, caption), (x, y, width, height)) in enumerate(zip(PANELS, rectangles)):
        ax = fig.add_axes([
            x / PAGE_WIDTH,
            1 - (y + height) / page_height,
            width / PAGE_WIDTH,
            height / page_height,
        ])
        ax.set_gid(chr(97 + index))
        image = exposure.draw_map(
            ax, exposure.absolute_classes(arrays[name]), chr(97 + index),
            caption, cmap, norm, geometry,
        )
        for text in ax.texts:
            text.set_visible(False)
        annotate_regions(ax, geometry)
        label_x = MAP_LEFT - 11.6 if index == 0 else x - 11.6
        fig.text(label_x / PAGE_WIDTH, 1 - (y - 7) / page_height,
                 chr(97 + index), fontsize=14, fontweight='normal')
        fig.text((x + width / 2) / PAGE_WIDTH, 1 - (y - 7) / page_height,
                 caption, ha='center', fontsize=11, fontweight='normal')
        axes.append(ax)
    bar = exposure.class_legend(
        fig, image,
        [(rectangles[2][0] + rectangles[2][2] / 2 - MAP_COLORBAR_WIDTH / 2) / PAGE_WIDTH,
         1 - (bar_top + COLORBAR_HEIGHT) / page_height,
         MAP_COLORBAR_WIDTH / PAGE_WIDTH, COLORBAR_HEIGHT / page_height],
        exposure.ABSOLUTE_LABELS, 'Mean annual EDH30',
    )
    bar.ax.tick_params(labelsize=8, pad=3)
    bar.set_label('Mean annual EDH30 (\N{DEGREE SIGN}C h)', fontsize=9, labelpad=2)
    return axes, bar


def validate_layout(fig, axes, colorbars):
    # 验收最终绘图区，而非仅检查分配给坐标轴的矩形。
    fig.canvas.draw()
    boxes = np.array([ax.get_position().extents for ax in axes])
    boxes *= [PAGE_WIDTH, fig.get_figheight() * 72, PAGE_WIDTH, fig.get_figheight() * 72]
    np.testing.assert_allclose(boxes[1, [0, 2]], boxes[3, [0, 2]], atol=.01)
    np.testing.assert_allclose(boxes[2, [0, 2]], boxes[4, [0, 2]], atol=.01)
    np.testing.assert_allclose(boxes[1, [1, 3]], boxes[2, [1, 3]], atol=.01)
    np.testing.assert_allclose(boxes[3, [1, 3]], boxes[4, [1, 3]], atol=.01)
    np.testing.assert_allclose(boxes[0, [0, 2]].mean(), boxes[1, [0, 2]].mean(), atol=.01)
    np.testing.assert_allclose(boxes[0, 2] - boxes[0, 0],
                               (boxes[1, 2] - boxes[1, 0]) * HISTORICAL_SCALE, atol=.01)
    np.testing.assert_allclose(boxes[2, 0] - boxes[1, 2], MAP_GAP, atol=.01)
    for bar, expected_width in zip(colorbars, (MAP_COLORBAR_WIDTH, YIELD_COLORBAR_WIDTH)):
        box = bar.ax.get_position(original=True)
        assert np.isclose(box.width * PAGE_WIDTH, expected_width)
        assert np.isclose(box.height * fig.get_figheight() * 72, COLORBAR_HEIGHT)
    renderer = fig.canvas.get_renderer()
    # 逐图检查全部区域缩写，且只比较文字本身，不将引线计入文字边框。
    for ax in axes[:5]:
        labels = [text for text in ax.texts if (text.get_gid() or '').startswith('region-')]
        assert {text.get_text() for text in labels} == set(heatmaps.LABELS)
        bounds = [Text.get_window_extent(text, renderer) for text in labels]
        for index, box in enumerate(bounds):
            assert ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
            assert not any(box.overlaps(other) for other in bounds[index + 1:])
    for text in fig.texts:
        box = text.get_window_extent(renderer)
        assert 0 <= box.x0 <= box.x1 <= fig.bbox.width, text.get_text()
        assert 0 <= box.y0 <= box.y1 <= fig.bbox.height, text.get_text()
    return boxes.tolist()


def export_figure(fig, stem):
    fig.savefig(f'{stem}.pdf', dpi=EXPORT_DPI, facecolor='white')
    fig.savefig(f'{stem}.svg', dpi=EXPORT_DPI, facecolor='white')
    paths = [stem.with_suffix('.pdf'), stem.with_suffix('.svg')]
    for extension in ('png', 'tif'):
        path = stem.with_suffix('.' + extension)
        options = {'pil_kwargs': {'compression': 'tiff_lzw'}} if extension == 'tif' else {}
        fig.savefig(path, dpi=EXPORT_DPI, facecolor='white', **options)
        paths.append(path)
    plt.close(fig)
    return paths


def main():
    # 3. 复用经校验的五波段数据；不重跑气象处理或 RF。
    raster = HEAT / 'EDH30_map_data.tif'
    heat_manifest = json.loads((HEAT / 'calculation.json').read_text(encoding='utf-8'))
    assert exposure.digest(raster) == heat_manifest['outputs_sha256'][raster.name]
    assert exposure.digest(HEAT / 'plot_heat_exposure.py') == heat_manifest['script_sha256']
    for path in (exposure.WORLD, exposure.REGIONS):
        record = next(row for row in heat_manifest['boundaries'].values() if Path(row['path']).name == path.name)
        assert exposure.digest(path) == record['sha256']
    with rasterio.open(raster) as src:
        assert src.shape == (720, 1440) and src.count == 5
        assert src.descriptions == tuple(name for name, _ in PANELS)
        arrays = {name: src.read(src.descriptions.index(name) + 1, masked=True).filled(np.nan)
                  for name, _ in PANELS}
    mask = np.isfinite(arrays[PANELS[0][0]])
    for layer in arrays.values():
        np.testing.assert_array_equal(np.isfinite(layer), mask)
        assert (layer[mask] >= 0).all()
    summary = heatmaps.read_summary()
    heatmaps.configure_style()
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.weight': 'normal',
        'axes.labelweight': 'normal',
        'axes.titleweight': 'normal',
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
    })
    geometry = list(exposure.map_geometry())
    # 先裁切底图几何，避免 PDF 中保留显示域以外的隐藏边界路径；栅格不变。
    display_boundary = Polygon(geometry[4].vertices)
    for index in (1, 2):
        geometry[index] = geometry[index].set_geometry(geometry[index].geometry.make_valid())
        geometry[index] = geometry[index].clip(display_boundary)
    rectangles, bar_top = map_layout(geometry)
    map_bottom = rectangles[-1][1] + rectangles[-1][3]
    map_height = map_bottom + 12

    # 4. 地图独立输出与整图复用同一绘图函数，保证面板内容和尺寸一致。
    fig = plt.figure(figsize=(PAGE_WIDTH / 72, map_height / 72))
    axes, bar = draw_maps(fig, arrays, geometry, rectangles, bar_top)
    validate_layout(fig, axes, [bar])
    outputs = export_figure(fig, HEAT / 'Extended_Data_Figure3a_e')

    # 热图表头紧接末行地图；保留单元格尺寸、区域顺序及全部数值。
    shift = map_bottom + 55 - heatmaps.HEATMAP_RECTS[0][1]
    page_height = 703.8 + shift
    fig = plt.figure(figsize=(PAGE_WIDTH / 72, page_height / 72))
    axes, map_bar = draw_maps(fig, arrays, geometry, rectangles, bar_top)
    yield_bar, plotted = heatmaps.draw_heatmaps(
        fig, summary, PAGE_WIDTH, page_height,
        offset_y=-shift,
        panel_labels=('f', 'g'),
    )
    boxes = validate_layout(fig, axes, [map_bar, yield_bar])
    outputs.extend(export_figure(fig, HOME / 'Extended_Data_Figure3'))
    report = {
        'panels': dict(zip('abcdefg', [*[name for name, _ in PANELS], 'rainfed_N_gain', 'irrigated_N_gain'])),
        'map_cells_per_panel': int(mask.sum()),
        'heatmap_cells': len(plotted),
        'layout': {
            'figure_points': [PAGE_WIDTH, page_height],
            'map_only_points': [PAGE_WIDTH, map_height],
            'map_plot_bounds_points': boxes,
            'map_rows': [['a'], ['b', 'c'], ['d', 'e']],
            'historical_width_relative_to_future': HISTORICAL_SCALE,
            'heat_exposure_colorbar_position': 'Upper right, alongside historical map',
            'font': 'Arial',
            'font_weight': 'normal',
            'colorbar_points': {
                'heat_exposure': [MAP_COLORBAR_WIDTH, COLORBAR_HEIGHT],
                'yield_gain': [YIELD_COLORBAR_WIDTH, COLORBAR_HEIGHT],
            },
            'raster_dpi': EXPORT_DPI,
        },
        'map_display_bounds': {'longitude': [-150, 150], 'latitude': [-60, 75], 'projection': 'Robinson'},
        'inputs_sha256': {str(p.relative_to(HOME.parent.parent)): exposure.digest(p) for p in [
            raster, heatmaps.SOURCE, HEAT / 'plot_heat_exposure.py',
            heatmaps.HOME / 'plot_yield_heatmaps.py', exposure.WORLD, exposure.REGIONS,
        ]},
        'script_sha256': exposure.digest(Path(__file__)),
        'outputs_sha256': {p.relative_to(HOME).as_posix(): exposure.digest(p) for p in outputs},
        'note': 'Five period means. Historical ERA5-Land and future ACCESS-CM2 retain their respective EDH30 methods. Yield heatmaps retain their package-specific spatial supports.',
    }
    (HOME / 'Extended_Data_Figure3_reproduction.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS: five mean-EDH30 maps, 240 yield-gain cells; Extended_Data_Figure3 PDF/PNG/SVG/TIFF.')


if __name__ == '__main__':
    main()
