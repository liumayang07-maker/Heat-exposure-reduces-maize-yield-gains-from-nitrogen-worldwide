"""Plot annual regional mean seasonal N2O-N emissions for two climates."""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# 1. 复用年度站点均值；单季总排放取成熟收获累计值，不累加逐日累计列。
HERE = Path(__file__).resolve().parent
YEARS = np.arange(2011, 2023)
REGIONS = ['NCP','NEC','NWC','SWC','NAMI','NAMR','EUR','SAM','SAS','WAS','EAF','WAF']
MODES = ['F1_mzI0','F1_mzI2','F2_mzI2']
LABELS = ['Basal-Rainfed','Basal-Irrigated','Split-Irrigated']
COLORS = ['#66A5ED','#F6C45C','#7758A3']
SCENARIOS = ['Baseline','Tmax+2']
annual = pd.read_csv(HERE/'regional_annual_site_mean_sd.csv')
seasons = pd.read_csv(HERE/'paired_season_endpoints.csv')
keys = ['region','mode','sowing_year','climate']
assert len(annual) == 864 and not annual.duplicated(keys).any()
assert len(seasons) == 5826
assert not seasons.duplicated(['xSite','mode','sowing_year','climate']).any()
assert np.isfinite(annual.annual_site_mean).all() and annual.annual_site_mean.ge(0).all()
assert annual.groupby(['region','mode','sowing_year']).n_sites.nunique().eq(1).all()
assert set(annual.climate) == set(SCENARIOS)
assert set(annual.region) == set(REGIONS) and set(annual['mode']) == set(MODES)

# 2. 从两情景各2,913个季末端点独立复核年度均值、样本SD及站点数。
check = seasons.groupby(keys).seasonal_N2O.agg(['mean','std','count'])
stored = annual.set_index(keys)[['annual_site_mean','annual_site_sd','n_sites']]
assert np.allclose(check.sort_index(),stored.sort_index(),equal_nan=True)
assert seasons.groupby(['xSite','mode','sowing_year']).climate.nunique().eq(2).all()

# 3. 与产量图保持相同布局、编号和配色；两情景共用相同的纵轴范围。
plt.rcParams.update({
    'font.family':'Arial',
    'font.size':8,
    'mathtext.fontset':'custom',
    'mathtext.rm':'Arial',
    'mathtext.it':'Arial:italic',
    'pdf.fonttype':42,
    'svg.fonttype':'none',
})
ymax = max(2., float(np.ceil(annual.annual_site_mean.max()/2)*2))
# 4. 颜色编码管理措施，线型编码气候情景；每个区域保留六条完整曲线。
climates = [('Baseline','-'),('Tmax+2','--')]
fig, axes = plt.subplots(
    3,
    4,
    figsize=(297/25.4,210/25.4),
    sharex=True,
    sharey=True,
)
for panel, (axis, region) in enumerate(zip(axes.flat, REGIONS)):
    block = annual[annual.region.eq(region)]
    for mode, color in zip(MODES,COLORS):
        for climate, linestyle in climates:
            values = block[block['mode'].eq(mode)&block.climate.eq(climate)]
            values = values.set_index('sowing_year').reindex(YEARS)
            assert values.annual_site_mean.notna().all()
            axis.plot(
                YEARS,
                values.annual_site_mean,
                color=color,
                linestyle=linestyle,
                lw=1.25,
                marker='o',
                ms=2.5,
                markerfacecolor=color if climate == 'Baseline' else 'white',
                markeredgewidth=.7,
            )
    assert len(axis.lines) == 6
    assert [line.get_linestyle() for line in axis.lines] == ['-','--']*3
    low, high = int(block.n_sites.min()), int(block.n_sites.max())
    count = str(low) if low == high else f'{low}\u2013{high}'
    # 5. 对十二个年度站点均值等权平均，透明三线表避开实际曲线范围。
    mean_table = block.groupby(['climate','mode']).annual_site_mean.mean()
    assert block.groupby(['climate','mode']).sowing_year.nunique().eq(12).all()
    table_values = [
        ['', 'BR', 'BI', 'SI'],
        ['Baseline'] + [f'{mean_table.loc[("Baseline", mode)]:.2f}' for mode in MODES],
        ['Tmax+2\u00b0C'] + [f'{mean_table.loc[("Tmax+2", mode)]:.2f}' for mode in MODES],
    ]
    table_bottom = .07 if block.annual_site_mean.max() > ymax * .75 else .76
    table_height = .23
    # NAMI 的变化幅度较大，使用底部较窄的空白，避免遮挡最低排放点。
    if region == 'NAMI':
        table_bottom = .015
        table_height = .15
    assert (
        (table_bottom + table_height) * ymax < block.annual_site_mean.min()
        or table_bottom * ymax > block.annual_site_mean.max()
    ), f'Table overlaps curves: {region}'
    table = axis.table(
        cellText=table_values,
        cellLoc='center',
        colWidths=[.34, .22, .22, .22],
        bbox=[.05, table_bottom, .90, table_height],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    for (row, column), cell in table.get_celld().items():
        cell.visible_edges = 'TB' if row == 0 else 'B' if row == 2 else ''
        cell.set_edgecolor('#444444')
        cell.set_linewidth(.5)
        cell.set_facecolor('none')
        cell.PAD = .03
        cell.get_text().set_fontfamily('Arial')
        if column == 0:
            cell.get_text().set_ha('left')
    axis.text(0,1.035,chr(ord('a')+panel),transform=axis.transAxes,
              fontsize=10,fontweight='bold',ha='left',va='bottom')
    axis.set_title(f'{region} | sites={count}',loc='left',x=.09,fontsize=8.5,pad=5)
    axis.set_xlim(2010.6,2022.4)
    axis.set_xticks([2011,2014,2017,2020,2022])
    axis.set_ylim(0,ymax)
    axis.set_yticks(np.arange(0,ymax+.1,2))
    axis.tick_params(labelsize=7,length=3,pad=2)
    axis.grid(axis='y',color='#eeeeee',lw=.5)
    # 表格保持透明，仅隐藏穿过表格高度范围的辅助网格线。
    for tick, gridline in zip(axis.get_yticks(), axis.get_ygridlines()):
        if table_bottom <= tick / ymax <= table_bottom + table_height:
            gridline.set_visible(False)
    axis.spines[['top','right']].set_visible(False)
fig.supxlabel('Sowing year',y=.055,fontsize=9)
fig.supylabel(r'Mean seasonal N$_2$O-N emissions (kg N ha$^{-1}$)',x=.02,fontsize=9)
handles = [Line2D([],[],color=c,lw=1.25,label=f'{l} ({abbr})')
           for c,l,abbr in zip(COLORS,LABELS,['BR','BI','SI'])]
handles += [
    Line2D([],[],color='#444444',linestyle='-',lw=1.25,label='Baseline'),
    Line2D([],[],color='#444444',linestyle='--',lw=1.25,label='Tmax +2\u00b0C'),
]
fig.legend(
    handles=handles,
    loc='lower center',
    bbox_to_anchor=(.5,.009),
    ncol=5,
    frameon=False,
    fontsize=8,
)
fig.subplots_adjust(left=.075,right=.985,top=.96,bottom=.11,wspace=.16,hspace=.22)
stem = 'Supplementary_N2O_annual'
for extension in ['png','pdf','svg']:
    fig.savefig(HERE/f'{stem}.{extension}',dpi=300,facecolor='white')
fig.savefig(HERE/f'{stem}.tif',dpi=300,facecolor='white',
            pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
print(f'{stem}: 12 regions, 6 curves each; solid baseline, dashed warming; Y range 0-{ymax:g}')
print('Verified: 5826 seasonal endpoints; 864 annual means; no input tables rewritten')

