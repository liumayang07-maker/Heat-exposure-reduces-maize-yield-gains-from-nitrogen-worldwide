"""Reproduce Extended Data Figure 3 f-g from the archived regional nitrogen-gain source."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd

HOME = Path(__file__).resolve().parent
SOURCE = HOME.parent / 'prediction/summary/regional_N_gain_source.csv'
INPUT_MANIFEST = HOME.parent / 'prediction/output_manifest.json'
HEATMAP_RECTS = [(40.425, 443.27, 288.2, 195.2), (366.3158, 443.27, 288.2, 195.2)]
ORDER = ['North_China_Plain', 'Northeast_China', 'Northwest_China',
         'Southwest_China', 'North_America_Irrigated', 'North_America_Rainfed',
         'Europe', 'South_America', 'South_Asia', 'West_Asia_North_Africa',
         'East_Africa', 'West_Africa']
LABELS = ['NCP', 'NEC', 'NWC', 'SWC', 'NAMI', 'NAMR', 'EUR', 'SAM', 'SAS', 'WAS', 'EAF', 'WAF']
GROUPS = ['historical_2011_2022', 'ssp245_2030_2050', 'ssp585_2030_2050',
          'ssp245_2070_2090', 'ssp585_2070_2090']
TITLES = ['Historical\n2011-2022', 'SSP2-4.5\n2030-2050', 'SSP5-8.5\n2030-2050',
          'SSP2-4.5\n2070-2090', 'SSP5-8.5\n2070-2090']


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_summary():
    # 1. 直接引用共享作图输入，核对校验值及全部 240 种组合。
    manifest = json.loads(INPUT_MANIFEST.read_text(encoding='utf-8'))
    record = next(row for row in manifest['files'] if row['path'] == 'summary/regional_N_gain_source.csv')
    assert digest(SOURCE) == record['sha256']
    summary = pd.read_csv(SOURCE, float_precision='round_trip')
    keys = ['region_id', 'dataset', 'irrigated', 'N_rate']
    assert len(summary) == 240 and not summary.duplicated(keys).any()
    expected = pd.MultiIndex.from_product([ORDER, GROUPS, [0, 1], [120, 240]], names=keys)
    assert summary.set_index(keys).index.sort_values().equals(expected.sort_values())
    assert np.isfinite(summary.mean_Mg_ha).all()
    assert summary.mean_Mg_ha.between(0, 7).all()
    assert (summary.years == np.where(summary.dataset.eq('historical_2011_2022'), 12, 21)).all()
    return summary


def configure_style():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 9.5,
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
    })


def draw_heatmaps(fig, summary, page_width, page_height, offset_y=0, panel_labels=('f', 'g')):
    # 2. 提取原最终拼图中的热图步骤；独立出图只平移纵坐标，不改变单元格尺寸。
    PAGE_WIDTH, PAGE_HEIGHT = page_width, page_height
    colorbar_width = 314.5 / PAGE_WIDTH
    colorbar_left = (341.85 - 314.5 / 2) / PAGE_WIDTH
    colorbar_height = 6 / PAGE_HEIGHT
    bounds = np.linspace(0, 7, 11)
    cmap = plt.get_cmap('viridis', 10)
    norm = BoundaryNorm(bounds, cmap.N)
    plotted = []
    for col, water in enumerate(('Rainfed', 'Irrigated')):
        ax = fig.add_axes([HEATMAP_RECTS[col][0] / PAGE_WIDTH,
                                    1 - (HEATMAP_RECTS[col][1] - offset_y + HEATMAP_RECTS[col][3]) / PAGE_HEIGHT,
                                    HEATMAP_RECTS[col][2] / PAGE_WIDTH,
                                    HEATMAP_RECTS[col][3] / PAGE_HEIGHT])
        subset = summary.loc[summary.irrigated.eq(col)].set_index(['region_id', 'dataset', 'N_rate'])
        matrix = np.array([[subset.loc[(r, g, n), 'mean_Mg_ha'] for g in GROUPS for n in (120, 240)] for r in ORDER])
        assert np.isfinite(matrix).all() and matrix.min() >= 0 and matrix.max() <= 7
        plotted.extend(matrix.ravel().tolist())
        image = ax.imshow(matrix, cmap=cmap, norm=norm, aspect='auto', interpolation='nearest')
        ax.set_yticks(range(12), LABELS)
        ax.set_xticks(range(10), ['N120', 'N240'] * 5)
        ax.tick_params(axis='y', length=0, pad=4, labelsize=9)
        ax.tick_params(axis='x', length=0, pad=4, labelsize=8)
        ax.xaxis.tick_top()
        for j, title in enumerate(TITLES):
            ax.text((j * 2 + 1) / 10, 1.085, title, transform=ax.transAxes,
                    ha='center', va='bottom', fontsize=8, linespacing=1.2)
        for y in np.arange(.5, 12, 1):
            ax.axhline(y, color='white', lw=.5)
        for x in np.arange(.5, 10, 1):
            ax.axvline(x, color='white', lw=2.8 if x == 1.5 else (1.7 if x % 2 == 1.5 else .5))
        for i in range(12):
            for j in range(10):
                luminance = np.dot(cmap(norm(matrix[i, j]))[:3], [.2126, .7152, .0722])
                ax.text(j, i, f'{matrix[i, j]:.2f}', ha='center', va='center', fontsize=7.3932,
                        color='white' if luminance < .46 else '#151515')
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.text((7.5889 if col == 0 else 339.3089) / PAGE_WIDTH,
                 1 - (426.7385 - offset_y) / PAGE_HEIGHT, panel_labels[col], fontsize=14)
        ax.set_xlabel(f'{water}: \N{GREEK CAPITAL LETTER DELTA}Grain yield from N0 to N120 and N240',
                      fontsize=9, labelpad=4)
    bar = fig.colorbar(image, cax=fig.add_axes([colorbar_left, 1 - (668.2 - offset_y) / PAGE_HEIGHT, colorbar_width, colorbar_height]),
                       orientation='horizontal', boundaries=bounds, ticks=bounds, drawedges=True)
    bar.ax.tick_params(labelsize=8, length=2)
    bar.set_label('Predicted yield difference relative to N0 (Mg/ha)', fontsize=9, labelpad=2)
    assert len(plotted) == 240
    return bar, plotted



def save_figure(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.texts:
        box = text.get_window_extent(renderer)
        assert 0 <= box.x0 <= box.x1 <= fig.bbox.width
        assert 0 <= box.y0 <= box.y1 <= fig.bbox.height
    paths = []
    for suffix in ('pdf', 'png', 'svg'):
        path = HOME / f'{stem}.{suffix}'
        fig.savefig(path, dpi=400, facecolor='white')
        paths.append(path)
    plt.close(fig)
    return paths


def main():
    summary = read_summary()
    configure_style()
    fig = plt.figure(figsize=(673.8 / 72, 303.8 / 72))
    bar, plotted = draw_heatmaps(fig, summary, 673.8, 303.8, offset_y=400)
    fig.canvas.draw()
    box = bar.ax.get_position(original=True)
    assert np.isclose(box.width * 673.8, 314.5) and np.isclose(box.height * 303.8, 6)
    # 3. 仅输出图件及最小复现记录，不另存一份 CSV。
    paths = save_figure(fig, 'Extended_Data_Figure3f_g')
    report = {
        'input': '../prediction/summary/regional_N_gain_source.csv',
        'input_sha256': digest(SOURCE),
        'input_manifest_sha256': digest(INPUT_MANIFEST),
        'cells': len(plotted),
        'region_order': ORDER,
        'datasets': GROUPS,
        'N_rates': [120, 240],
        'range_Mg_ha': [min(plotted), max(plotted)],
        'layout': {'figure_points': [673.8, 303.8], 'colorbar_points': [314.5, 6], 'font': 'Arial'},
        'script_sha256': digest(Path(__file__)),
        'outputs_sha256': {path.name: digest(path) for path in paths},
    }
    (HOME / 'calculation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('PASS: Extended Data Figure 3 f-g PNG/PDF/SVG; 240 verified heatmap values; no RF inference.')


if __name__ == '__main__':
    main()
