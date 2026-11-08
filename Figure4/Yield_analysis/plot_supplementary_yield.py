"""Extract mature-harvest yield and plot regional annual site means."""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# 1. 共享已核验的成熟收获端点；单位经用户确认为 t ha^-1，不做换算。
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
YEARS = np.arange(2011, 2023)
REGIONS = ['NCP','NEC','NWC','SWC','NAMI','NAMR','EUR','SAM','SAS','WAS','EAF','WAF']
MODES = ['F1_mzI0','F1_mzI2','F2_mzI2']
LABELS = ['Basal-Rainfed','Basal-Irrigated','Split-Irrigated']
COLORS = ['#66A5ED','#F6C45C','#7758A3']
SCENARIOS = [
    ('Baseline', 'baseline_daily.csv.gz'),
    ('Tmax+2', 'tmax_plus2_daily.csv.gz'),
]
KEYS = ['xSite','region','mode','sowing_year']
endpoints = pd.read_csv(ROOT/'N2O_analysis/paired_season_endpoints.csv')
assert len(endpoints) == 5826
assert not endpoints.duplicated(KEYS+['climate']).any()
records = []
for climate, filename in SCENARIOS:
    source = pd.read_csv(
        ROOT/'APSIM_outputs'/filename,
        usecols=['xSite','mzManagement','Date()','Maize.stage','mzGrain','xIW','mzNrate'],
    )
    selected = endpoints[endpoints.climate.eq(climate)].reset_index(drop=True)
    rows = source.iloc[selected.source_harvest_row.to_numpy()-2].reset_index(drop=True)
    assert len(rows) == 2913
    assert rows.xSite.equals(selected.xSite)
    assert rows.mzManagement.str.replace('mzF','F',regex=False).equals(selected['mode'])
    assert rows['Maize.stage'].eq('harvest_ripe').all()
    assert rows.xIW.eq('xIW70').all() and rows.mzNrate.eq('mzN240').all()
    assert pd.to_datetime(rows['Date()'],format='%d/%m/%Y').equals(pd.to_datetime(selected.harvest_date))
    assert np.isfinite(rows.mzGrain).all() and rows.mzGrain.ge(0).all()
    selected = selected[KEYS+['climate','harvest_date','source_harvest_row']].copy()
    selected['yield_t_ha'] = rows.mzGrain.to_numpy()
    records.append(selected)
site_yields = pd.concat(records, ignore_index=True)
paired = site_yields.groupby(KEYS).climate.nunique()
assert len(paired) == 2913 and paired.eq(2).all()

# 2. 每个区域、播种年度和管理下对站点等权平均；有效零产量保留。
group_keys = ['region','mode','sowing_year','climate']
annual = site_yields.groupby(group_keys).agg(
    mean_yield_t_ha=('yield_t_ha','mean'),
    sd_yield_t_ha=('yield_t_ha','std'),
    n_sites=('xSite','nunique'),
).reset_index()
assert len(annual) == 864
assert annual.groupby(['region','mode','sowing_year']).n_sites.nunique().eq(1).all()
independent = site_yields.groupby(group_keys).yield_t_ha.apply(lambda x: np.mean(x.to_numpy()))
assert np.allclose(independent, annual.set_index(group_keys).mean_yield_t_ha)
site_yields.to_csv(HERE/'yield_site_seasons.csv',index=False,encoding='utf-8-sig')
annual.to_csv(HERE/'yield_regional_annual.csv',index=False,encoding='utf-8-sig')

# 3. 两情景分别成图，保持相同布局、颜色和统一 Y 轴，便于直接比较。
plt.rcParams.update({
    'font.family':'Arial',
    'font.size':8,
    'mathtext.fontset':'custom',
    'mathtext.rm':'Arial',
    'mathtext.it':'Arial:italic',
    'pdf.fonttype':42,
    'svg.fonttype':'none',
})
ymax = max(2., float(np.ceil(annual.mean_yield_t_ha.max()/2)*2))
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
            assert values.mean_yield_t_ha.notna().all()
            axis.plot(
                YEARS,
                values.mean_yield_t_ha,
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
    # 5. 在曲线之外的空白位置内嵌透明三线表，保持所有数据可见。
    mean_table = block.groupby(['climate','mode']).mean_yield_t_ha.mean()
    assert block.groupby(['climate','mode']).sowing_year.nunique().eq(12).all()
    table_values = [
        ['', 'BR', 'BI', 'SI'],
        ['Baseline'] + [f'{mean_table.loc[("Baseline", mode)]:.2f}' for mode in MODES],
        ['Tmax+2\u00b0C'] + [f'{mean_table.loc[("Tmax+2", mode)]:.2f}' for mode in MODES],
    ]
    # 高产区使用左下空白，低产区使用左上空白；不覆盖曲线。
    table_bottom = .07 if block.mean_yield_t_ha.max() > ymax * .75 else .76
    table = axis.table(
        cellText=table_values,
        cellLoc='center',
        colWidths=[.34, .22, .22, .22],
        bbox=[.05, table_bottom, .90, .23],
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
    axis.text(.09,1.035,f'{region} | sites={count}',transform=axis.transAxes,
              fontsize=8.5,ha='left',va='bottom')
    axis.set_xlim(2010.6,2022.4)
    axis.set_xticks([2011,2014,2017,2020,2022])
    axis.set_ylim(0,ymax)
    axis.set_yticks(np.arange(0,ymax+.1,2))
    axis.tick_params(labelsize=7,length=3,pad=2)
    axis.grid(axis='y',color='#eeeeee',lw=.5)
    axis.spines[['top','right']].set_visible(False)
fig.supxlabel('Sowing year',y=.055,fontsize=9)
fig.supylabel(r'Mean grain yield (t ha$^{-1}$)',x=.02,fontsize=9)
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
stem = 'Supplementary_yield_annual'
for extension in ['png','pdf','svg']:
    fig.savefig(HERE/f'{stem}.{extension}',dpi=300,facecolor='white')
fig.savefig(HERE/f'{stem}.tif',dpi=300,facecolor='white',
            pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
print(f'{stem}: 12 regions, 6 curves each; solid baseline, dashed warming; Y range 0-{ymax:g}')
print(f'Verified: {len(site_yields)} site-seasons; {len(annual)} regional annual means; valid zeros retained')
