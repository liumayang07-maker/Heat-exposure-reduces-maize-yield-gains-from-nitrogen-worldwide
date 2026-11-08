# 绘图入口：共享归档结果，直接运行 Python。
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.font_manager import FontProperties, findfont
from shapely.geometry import box
from pyproj import Transformer
from pypdf import PdfReader
import xml.etree.ElementTree as ET

# 1. 复用归档季节端点，按情景筛选完整12年的站点并汇总；曲线及地图不变。
root = Path(__file__).resolve().parent.parent
project = root.parent.parent
base = root / 'N2O_analysis'
inputs = [base/'regional_bins.csv', base/'stage_positions.csv',
          root/'N2O_analysis/site_metadata.csv',
          root.parent/'Figure2/Figure2a/study_regions.geojson',
          root.parent/'Figure2/Figure2a/world_boundaries.gpkg',
          root/'N2O_analysis/regional_12year_sum.csv',
          base/'paired_season_endpoints.csv', base/'season_baseline.csv']
hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
summary = pd.read_csv(inputs[0])
stages = pd.read_csv(inputs[1])
sites = pd.read_csv(inputs[2])
boundaries = gpd.read_file(inputs[3]).to_crs(4326)
world = gpd.read_file(inputs[4])
previous_totals = pd.read_csv(inputs[5]).set_index(['region', 'mode'])
# 两情景分别使用完整站点；基线恢复配对表曾排除的有效2013年记录。
paired_endpoints = pd.read_csv(inputs[6])
baseline_endpoints = pd.read_csv(inputs[7]).rename(columns={'mzN2Oatm': 'seasonal_N2O'})
baseline_endpoints['climate'] = 'Baseline'
season_keys = ['xSite', 'region', 'mode', 'climate']
columns = season_keys + ['sowing_year', 'seasonal_N2O']
seasonal = pd.concat([baseline_endpoints[columns],
                      paired_endpoints.loc[paired_endpoints.climate.eq('Tmax+2'), columns]], ignore_index=True)
assert not seasonal.duplicated(season_keys + ['sowing_year']).any()
years = set(range(2011, 2023))
complete = seasonal.groupby(season_keys).sowing_year.agg(lambda values: set(values) == years)
site_complete = complete.groupby(level=['xSite', 'climate']).agg(['all', 'size'])
eligible = site_complete.index[site_complete['all'] & site_complete['size'].eq(3)]
excluded = site_complete.index.difference(eligible)
excluded_sites = sorted(set(excluded.get_level_values('xSite')))
selected = pd.MultiIndex.from_frame(seasonal[['xSite', 'climate']]).isin(eligible)
site_totals = seasonal[selected].groupby(season_keys).seasonal_N2O.sum().reset_index()
site_totals = site_totals.rename(columns={'seasonal_N2O': 'cumulative_N2O'})
assert len(site_totals) == len(eligible) * 3
regional = site_totals.groupby(['region', 'mode', 'climate']).cumulative_N2O.agg(
    mean='mean', sd='std', n_sites='count').reset_index()
regional['se'] = regional['sd'] / np.sqrt(regional['n_sites'])
assert len(regional) == 72 and regional.n_sites.ge(2).all()
site_totals.to_csv(base/'site_12year_cumulative.csv', index=False)
regional.to_csv(base/'regional_12year_cumulative_mean_sd.csv', index=False)
totals = regional.pivot(index=['region', 'mode'], columns='climate', values='mean')
se = regional.pivot(index=['region', 'mode'], columns='climate', values='se')
totals = totals.join(se.add_suffix('_se'))
totals = totals.join(regional.pivot(index=['region', 'mode'], columns='climate', values='n_sites').add_suffix('_n'))
totals['percent_change'] = 100 * (totals['Tmax+2'] - totals['Baseline']) / totals['Baseline']
unaffected = ~previous_totals.index.get_level_values('region').isin(seasonal.loc[seasonal.xSite.isin(excluded_sites), 'region'])
np.testing.assert_allclose(totals.loc[previous_totals.index[unaffected], ['Baseline', 'Tmax+2']],
                           previous_totals.loc[unaffected, ['Baseline', 'Tmax+2']])
assert totals.index.is_unique and len(totals) == 36
assert np.isfinite(totals[['Baseline', 'Tmax+2']]).all().all()
assert (totals[['Baseline', 'Tmax+2']] >= 0).all().all()
target_mean_height = .50
tick_length_pt = 1.5 * 72 / 25.4
axis_linewidth_pt = 1.0
modes = ['F1_mzI0', 'F1_mzI2', 'F2_mzI2']
colors = ['#66A5ED', '#F6C45C', '#7758A3']
labels = ['Basal-Rainfed', 'Basal-Irrigated', 'Split-Irrigated']
heat_color = '#CE4C58'
projection = '+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs'
points = gpd.GeoDataFrame(sites, geometry=gpd.points_from_xy(sites.lon, sites.lat), crs=4326).to_crs(projection)
assert len(points) == 81 and points.xSite.nunique() == 81
edges = np.arange(101)

# 2. 中央地图与周围12个区域小图；布局接近地理位置，避免连接线穿越小图。
panel_width = .175
panel_height = .165
layout = {
    'NAMR': (.065, .735, panel_width, panel_height, 'top'),
    'EUR':  (.300, .735, panel_width, panel_height, 'top'),
    'NWC':  (.535, .735, panel_width, panel_height, 'top'),
    'NEC':  (.770, .735, panel_width, panel_height, 'top'),
    'NAMI': (.065, .525, panel_width, panel_height, 'left'),
    'SAM':  (.065, .325, panel_width, panel_height, 'left'),
    'NCP':  (.770, .525, panel_width, panel_height, 'right'),
    'SWC':  (.770, .325, panel_width, panel_height, 'right'),
    'WAF':  (.065, .120, panel_width, panel_height, 'bottom'),
    'EAF':  (.300, .120, panel_width, panel_height, 'bottom'),
    'WAS':  (.535, .120, panel_width, panel_height, 'bottom'),
    'SAS':  (.770, .120, panel_width, panel_height, 'bottom'),
}
assert set(layout) == set(summary.region.unique())
arial_path = findfont(FontProperties(family='Arial'), fallback_to_default=False)
plt.rcParams.update({'font.family':'Arial', 'font.size':10,
                     'font.weight':'normal', 'axes.labelweight':'normal',
                     'axes.titleweight':'normal', 'text.parse_math':False,
                     'svg.fonttype':'none', 'pdf.fonttype':42})
def draw_figure(kind, return_figure=False):
    assert kind in {'totals', 'increments', 'cumulative'}
    n2o_metric = 'mzN2Oatm_cumulative' if kind == 'cumulative' else 'mzN2Oatm_increment'
    heat_metric = 'EDH30_cumulative' if kind == 'cumulative' else 'EDH30_intensity'
    out = root if kind == 'totals' else base
    fig = plt.figure(figsize=(297/25.4, 210/25.4))
    # 按经度裁切实际地理要素，避免高纬度仍显示150度以外的区域。
    map_ax = fig.add_axes([.260 - 10 / 297, .315, .49, .39])
    longitude_window = box(-150, -90, 150, 90)
    cropped_world = world.to_crs(4326).copy()
    cropped_world.geometry = cropped_world.geometry.make_valid().intersection(longitude_window)
    cropped_regions = boundaries.copy()
    cropped_regions.geometry = cropped_regions.geometry.make_valid().intersection(longitude_window)
    cropped_world[~cropped_world.is_empty].to_crs(projection).plot(
        ax=map_ax, color='#ededed', edgecolor='#cccccc', linewidth=.3)
    cropped_regions.to_crs(projection).plot(ax=map_ax, color=boundaries.outline_color, alpha=.13, edgecolor='none')
    assert sites.lon.between(-150, 150).all()
    for region in layout:
        group = points[points.figure2_region_code.eq(region)]
        group.plot(ax=map_ax, color=group.outline_color, markersize=19,
                   edgecolor='white', linewidth=.45, zorder=5)
    project_map = Transformer.from_crs(4326, projection, always_xy=True)
    map_half_width = project_map.transform(150, 0)[0]
    map_ax.set_xlim(-map_half_width, map_half_width)
    map_ax.set_ylim(-6400000, 7800000)
    map_ax.set_axis_off()
    panel_axes = {}
    inset_axes = []
    counts = []

    def draw_steps(axis, block, metric, color, heat=False):
        mean = block[metric+'_mean'].to_numpy()
        sd = block[metric+'_std'].to_numpy()
        axis.stairs(mean, edges, baseline=None, color=color, lw=.8, zorder=3)
        axis.fill_between(edges, np.r_[mean-sd, (mean-sd)[-1]],
                          np.r_[mean+sd, (mean+sd)[-1]], step='post',
                          color=color, alpha=.055 if heat else .14, linewidth=0)

    def set_limits(axis, block, metric):
        mean = block[metric+'_mean']
        sd = block[metric+'_std'].fillna(0)
        lo = min(0, float((mean-sd).min()))
        hi = max(.001, float((mean+sd).max()))
        pad = (hi-lo)*.06
        ticks = MaxNLocator(nbins=5, steps=[1,2,2.5,5,10]).tick_values(lo-pad if lo<0 else 0, hi+pad)
        axis.set_yticks(ticks)
        # 沿用原区域图范围，不为内嵌图额外扩展双Y轴。
        axis.set_ylim(ticks[0], ticks[-1])

    for index, (region, (x,y,w,h,side)) in enumerate(layout.items()):
        if kind != 'totals':
            # 保持面板中心不变，为相邻双 Y 轴的刻度标签留出间距。
            side_inset = 6 / (fig.get_figwidth() * 72)
            x += side_inset
            w -= 2 * side_inset
        ax = fig.add_axes([x,y,w,h])
        block = summary[summary.region.eq(region)]
        group = block[block['mode'].eq(modes[0])].sort_values('bin')
        nsites = int(group.n_sites.iloc[0])
        n = int(group.mzN2Oatm_increment_count.iloc[0])
        same_heat = None
        if kind != 'totals':
            right = ax.twinx()
            # 主曲线及SD位于内嵌柱图上层；父轴背景透明，避免覆盖柱图。
            ax.set_zorder(9)
            right.set_zorder(10)
            ax.patch.set_visible(False)
            right.patch.set_visible(False)
            block = summary[summary.region.eq(region)]
            groups = [block[block['mode'].eq(m)].sort_values('bin') for m in modes]
            for group in groups:
                assert np.array_equal(group.bin.to_numpy(), np.arange(1,101))
            heat_columns = [heat_metric+'_'+s for s in ['mean','std','count']]
            same_heat = all(np.allclose(groups[0][heat_columns], g[heat_columns], equal_nan=True) for g in groups[1:])
            for mode, color, group in zip(modes, colors, groups):
                draw_steps(ax, group, n2o_metric, color)
                if not same_heat:
                    draw_steps(right, group, heat_metric, heat_color, True)
                for stage, marker in [('end_of_juvenil','v'), ('flag_leaf','^')]:
                    values = stages.loc[stages.region.eq(region)&stages['mode'].eq(mode)&stages.event.eq(stage), 'progress']
                    if len(values):
                        ax.plot(values.mean(), .025, marker=marker, color=color, ms=3,
                                transform=ax.get_xaxis_transform(), clip_on=True)
            if same_heat:
                draw_steps(right, groups[0], heat_metric, heat_color, True)
            set_limits(ax, block, n2o_metric)
            set_limits(right, block, heat_metric)
            ax.set_xlim(0,100)
            ax.set_xticks(np.arange(0, 101, 20))
            ax.set_xlabel('')
            ax.set_ylabel('')
            right.set_ylabel('')
            for axis in [ax, right]:
                axis.spines['top'].set_visible(False)
                axis.spines['bottom'].set_visible(axis is ax)
                axis.tick_params(axis='both', which='both', length=tick_length_pt,
                                 width=axis_linewidth_pt, labelsize=10, pad=3)
                axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
                for spine in axis.spines.values():
                    spine.set_linewidth(axis_linewidth_pt)
            right.tick_params(axis='y', colors='#333333')
            nsites = int(groups[0].n_sites.iloc[0])
            n = int(groups[0].mzN2Oatm_increment_count.iloc[0])
            assert nsites == int(sites.figure2_region_code.eq(region).sum())
            ax.text(.045, .97, region, transform=ax.transAxes,
                    ha='left', va='top', fontsize=11)
            ax.text(.045, .84, f'n={n}', transform=ax.transAxes,
                    ha='left', va='top', fontsize=10)
        else:
            # 各区域按相同规则独立截断，下方留出数据极差的25%，平均可见柱高约50%。
            inset = ax
            inset.patch.set_visible(False)
            values_region = totals.loc[region, ['Baseline', 'Tmax+2']].to_numpy()
            se_region = totals.loc[region, ['Baseline_se', 'Tmax+2_se']].to_numpy()
            error_min = (values_region - se_region).min()
            error_max = (values_region + se_region).max()
            data_range = error_max - error_min
            candidates = []
            for exponent in range(-1, 4):
                for factor in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8):
                    step = factor * 10 ** exponent
                    lower = np.floor((error_min - .10 * data_range) / step) * step
                    if error_min >= 0:
                        lower = max(0, lower)
                    upper = lower + 3 * step
                    if lower < error_min and error_max < upper - .25 * (upper - lower):
                        fraction = (values_region.mean() - lower) / (upper - lower)
                        candidates.append((abs(fraction - target_mean_height), lower, upper, step))
            assert candidates, region
            _, bar_ymin, ymax, step = min(candidates)
            span = ymax - bar_ymin
            for j, (mode, color) in enumerate(zip(modes, colors)):
                values = totals.loc[(region, mode)]
                inset.bar(j-.18, values['Baseline'], width=.32, facecolor='white',
                          edgecolor=color, linewidth=.65, yerr=values['Baseline_se'],
                          error_kw={'elinewidth': .7, 'capsize': 2, 'capthick': .7, 'ecolor': '#333333'})
                inset.bar(j+.18, values['Tmax+2'], width=.32, color=color, linewidth=0,
                          yerr=values['Tmax+2_se'],
                          error_kw={'elinewidth': .7, 'capsize': 2, 'capthick': .7, 'ecolor': '#333333'})
                percent = 100*(values['Tmax+2']-values['Baseline'])/values['Baseline'] if values['Baseline']>0 else np.nan
                if np.isfinite(percent):
                    assert np.isclose(percent, values['percent_change'])
                percent_text = f'{percent:+.1f}%' if np.isfinite(percent) else 'NA'
                inset.text(j, max(values['Baseline'] + values['Baseline_se'], values['Tmax+2'] + values['Tmax+2_se']) + span * .035, percent_text,
                           ha='center', va='bottom', fontsize=10)
            inset.set_xlim(-.6, 2.6)
            # 仅截断显示范围，保留原始累计排放值及变化率。
            ticks = bar_ymin + np.arange(4) * step
            inset.set_yticks(ticks)
            inset.set_ylim(bar_ymin, ymax)
            inset.set_xticks([0,1,2], ['BR','BI','SI'])
            inset.tick_params(axis='both', which='both', length=tick_length_pt,
                              width=axis_linewidth_pt, labelsize=10, pad=1)
            inset.spines[['top','right']].set_visible(False)
            for spine in inset.spines.values():
                spine.set_linewidth(axis_linewidth_pt)
            ax.text(.045, .97, region, transform=ax.transAxes, ha='left', va='top', fontsize=11)
        panel_axes[region] = ax
        counts.append(dict(region=region, sites=nsites, site_years_per_mode=n, heat_profiles_merged=same_heat))

    # 3. 连接线指向区域内站点中心；所有线在小图底层，数据曲线不受遮盖。
    fig.canvas.draw()
    panel_sizes = [ax.get_window_extent().size for ax in panel_axes.values()]
    assert np.allclose(panel_sizes, panel_sizes[0], atol=1e-8)
    if kind == 'totals':
        assert sum(len(ax.patches) for ax in panel_axes.values()) == 72
        assert len(fig.axes) == 13
    else:
        assert len(fig.axes) == 25
    # 各区域使用独立通道；底部外侧线先转弯，避免穿过内侧线的竖段。
    leader_lanes = {
        'NAMR': .709,
        'EUR': .715,
        'NWC': .709,
        'NEC': .703,
        'WAF': .300,
        'EAF': .291,
        'WAS': .291,
        'SAS': .300,
        'NAMI': .265,
        'SAM': .270,
        'NCP': .740,
        'SWC': .735,
    }
    # 主图沿用用户手动PDF的连线；坐标单位pt，原点为A4画布左上角。
    manual_routes_pt = {
        'NAMI': [(274.178, 233.380), (223.100, 233.380), (223.100, 233.646), (212.156, 233.646)],
        'SWC': [(533.877, 249.989), (619.115, 249.989), (619.115, 352.700), (626.315, 352.700)],
        'WAF': [(404.075, 281.929), (404.075, 412.650), (193.214, 412.650), (193.214, 416.692)],
        'EAF': [(443.715, 305.413), (443.715, 419.275), (326.233, 419.275)],
        'WAS': [(449.609, 240.363), (449.609, 419.571), (524.076, 419.571)],
        'SAS': [(503.074, 255.720), (503.074, 416.871), (657.095, 416.871), (657.095, 420.860)],
        'NAMR': [(294.147, 223.802), (294.147, 177.497), (126.334, 177.497), (126.334, 172.609)],
        'EUR': [(325.000, 172.609), (325.000, 176.905), (410.333, 176.905), (410.333, 215.103)],
        'NWC': [(520.333, 172.609), (520.333, 237.646)],
        'NEC': [(721.921, 172.609), (721.921, 178.387), (541.889, 178.387), (541.889, 223.238)],
        'NCP': [(631.878, 233.380), (539.704, 233.380)],
        'SAM': [(210.240, 352.005), (315.222, 352.005)],
    }
    leader_routes = {}
    for region, (x,y,w,h,side) in layout.items():
        group = points[points.figure2_region_code.eq(region)]
        center = (float(group.geometry.x.mean()), float(group.geometry.y.mean()))
        start = fig.transFigure.inverted().transform(map_ax.transData.transform(center))
        sx,sy = start
        if side == 'top':
            end = (x+w*.5, y-.008)
            lane = leader_lanes[region]
            route = [(sx,sy),(sx,lane),(end[0],lane),end]
        elif side == 'bottom':
            # 外侧面板从靠地图的一端接入，避开两侧面板的X轴刻度文字。
            anchor = {'WAF': .94, 'SAS': .06}.get(region, .5)
            end = (x+w*anchor, y+h+.008)
            lane = leader_lanes[region]
            route = [(sx,sy),(sx,lane),(end[0],lane),end]
        elif side == 'left':
            end = (x+w+.012, y+h*.5)
            lane = leader_lanes[region]
            route = [(sx,sy),(lane,sy),(lane,end[1]),end]
        else:
            end = (x-.012, y+h*.5)
            lane = leader_lanes[region]
            route = [(sx,sy),(lane,sy),(lane,end[1]),end]
        route = [(px / (297 / 25.4 * 72), 1 - py / (210 / 25.4 * 72))
                 for px, py in manual_routes_pt[region]]
        leader_routes[region] = route
        fig.add_artist(Line2D(*zip(*route), transform=fig.transFigure,
                             color='#999999', lw=.55, alpha=.8, zorder=.5))
    # 检查整条引导线之间既无共线重叠，也无交叉。
    from shapely.geometry import LineString
    from itertools import combinations
    for first, second in combinations(leader_routes, 2):
        assert not LineString(leader_routes[first]).intersects(LineString(leader_routes[second])), (first, second)
    # 仅保留各自图件需要的图例，两图地图数据与连接线相同。
    management_handles = [Line2D([], [], color=c, lw=1.5, label=f'{label} ({code})')
                          for c, label, code in zip(colors, labels, ['BR', 'BI', 'SI'])]
    handles = [Line2D([], [], color='#555555', marker='o', ls='', ms=4, label='Selected sites (n=81)')]
    if kind != 'totals':
        handles += [Line2D([], [], color='#555555', marker='v', ls='', ms=4, label='End of juvenile'),
                    Line2D([], [], color='#555555', marker='^', ls='', ms=4, label='Flag leaf')]
        heat_label = 'Cumulative EDH30' if kind == 'cumulative' else 'EDH30 intensity'
        handles += management_handles + [Line2D([], [], color=heat_color, lw=1.2, label=heat_label)]
        left_title = 'Cumulative N2O-N (kg N/ha)' if kind == 'cumulative' else 'N2O-N increment (kg N/ha/bin)'
        right_title = 'Cumulative EDH30 (°C h)' if kind == 'cumulative' else 'Mean EDH30 intensity (°C·h per hour)'
        fig.text(.018, .52, left_title,
                 rotation=90, va='center', ha='center', fontsize=11)
        fig.text(.989, .52, right_title,
                 rotation=270, va='center', ha='center', fontsize=11)
        fig.text(.5, .067, 'Relative thermal progress: 0–100% (left to right)', ha='center', fontsize=10)
    else:
        handles += [Patch(facecolor='white', edgecolor='#555555', label='Baseline'),
                    Patch(facecolor='#777777', label='Tmax +2°C')]
        handles += management_handles
        fig.text(.018, .52, '12-year cumulative N2O-N emissions (kg N/ha)',
                 rotation=90, va='center', ha='center', fontsize=11)
    fig.legend(handles=handles, loc='upper left',
               bbox_to_anchor=(.605 - 10 / 297, .515),
               ncol=1, frameon=False, fontsize=10, handlelength=1.5,
               handletextpad=.5, labelspacing=.8, borderaxespad=0)

    # 组合图入口直接复用数据绘图对象，不读取已导出的位图。
    if return_figure:
        return fig

    stem = {'totals': 'Figure4', 'increments': 'Extended_Data_N2O_increments',
            'cumulative': 'Extended_Data_Figure4'}[kind]
    for ext in ['png', 'pdf', 'svg', 'tif']:
        temporary = out / f'{stem}.tmp.{ext}'
        options = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tif' else {}
        fig.savefig(temporary, dpi=300, facecolor='white', **options)
        temporary.replace(out / f'{stem}.{ext}')
    plt.close(fig)
    assert sum(c['sites'] for c in counts) == 81
    assert sum(c['site_years_per_mode'] for c in counts) * 3 == 2916
    assert len(PdfReader(out / f'{stem}.pdf').pages) == 1
    ET.parse(out / f'{stem}.svg')
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in hashes.items())
    report = dict(status='passed', regions=12, sites=81, seasons=2916,
                  bars=72 if kind == 'totals' else 0, page_mm=[297, 210],
                  panel_mm=[297 * ax.get_position().width, 210 * panel_height],
                  font='Arial', leader_routes=leader_routes, input_hashes=hashes)
    if kind == 'totals':
        report.update(statistic='Mean +/- SE of scenario-specific complete-site 12-year totals',
                      analysis_sites=site_totals.groupby('climate').xSite.nunique().to_dict(),
                      excluded_incomplete_site_scenarios=[list(pair) for pair in excluded],
                      comparison_note='NAMI compares different site sets: baseline n=6, warming n=5; percentages are not paired warming effects.',
                      site_totals_rows=len(site_totals), regional_summary_rows=len(regional),
                      error_bars=72, error_bar_type='SE = sample SD / sqrt(n_sites)', sd_ddof=1)
    (base / f'{stem}_verification.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(f'{stem}: PASS; 81 sites, 12 regions, unchanged inputs')


if __name__ == '__main__':
    draw_figure('totals')
    draw_figure('increments')

