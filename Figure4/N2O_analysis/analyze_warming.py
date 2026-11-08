# 分析入口：沿用原计算方法，直接运行 Python。
import json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
HOME = ROOT/'APSIM_outputs'
OUT = ROOT/'N2O_analysis'
BASE = HOME/'baseline_daily.csv.gz'
WARM = HOME/'tmax_plus2_daily.csv.gz'
SEASONS = ROOT/'N2O_analysis/season_baseline.csv'
YEARS = list(range(2011,2023))
MODES = ['F1_mzI0','F1_mzI2','F2_mzI2']
OUT.mkdir(parents=True, exist_ok=True)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4194304),b''): h.update(block)
    return h.hexdigest()

def write(d,name):
    d.to_csv(OUT/(name+'.csv'),index=False,encoding='utf-8-sig')

hashes={str(p):sha(p) for p in [BASE,WARM,SEASONS]}
baseline=pd.read_csv(SEASONS,parse_dates=['first_record_date','harvest_date'])
excluded=baseline[baseline.xSite.eq('APSIM_012') & baseline.sowing_year.eq(2013)].copy()
assert len(excluded)==3
excluded['reason']='Warming season ends before mature harvest; remove both climate counterparts'
write(excluded,'excluded_pairs')
target=baseline.drop(excluded.index).copy()
keys=['xSite','mode','sowing_year']
assert len(target)==2913 and not target.duplicated(keys).any()
lookup=target.set_index(['xSite','mode','first_record_date'])
records=[]
for climate,path in [('Baseline',BASE),('Tmax+2',WARM)]:
    cols=['xSite','mzManagement','Date()','Maize.stage','mzN2Oatm','xIW','mzNrate']
    d=pd.read_csv(path,usecols=cols)
    d['source_row']=np.arange(len(d))+2
    d['date']=pd.to_datetime(d['Date()'],format='%d/%m/%Y')
    d['mode']=d.mzManagement.str.replace('mzF','F',regex=False)
    assert d.xIW.eq('xIW70').all() and d.mzNrate.eq('mzN240').all()
    for (site,mode),g in d.groupby(['xSite','mode']):
        g=g.sort_values('date')
        start=g['Maize.stage'].eq('sowing') & ~g['Maize.stage'].shift().eq('sowing')
        for _,s in g.groupby(start.cumsum()):
            identity=(site,mode,s.date.iloc[0])
            if identity not in lookup.index: continue
            row=lookup.loc[identity]
            ripe=s[s['Maize.stage'].eq('harvest_ripe')]
            assert len(ripe)==1, (climate,identity,'harvest not unique')
            end=ripe.iloc[0]
            s=s[s.date<=end.date]
            assert s.date.is_unique and len(s)==(end.date-s.date.iloc[0]).days+1
            assert np.isfinite(s.mzN2Oatm).all() and s.mzN2Oatm.ge(0).all()
            assert s.mzN2Oatm.diff().dropna().ge(-1e-8).all()
            if climate=='Baseline':
                assert np.isclose(end.mzN2Oatm,row.mzN2Oatm)
                assert end.date==row.harvest_date
            records.append(dict(xSite=site,region=row.region,mode=mode,sowing_year=int(row.sowing_year),
                climate=climate,first_record_date=s.date.iloc[0],harvest_date=end.date,
                source_first_row=int(s.source_row.iloc[0]),source_harvest_row=int(end.source_row),
                seasonal_N2O=float(end.mzN2Oatm)))
    print(climate,'harvest extraction complete',flush=True)
seasons=pd.DataFrame(records)
assert not seasons.duplicated(keys+['climate']).any()
paired=seasons.pivot(index=keys+['region'],columns='climate',values='seasonal_N2O').reset_index()
assert len(paired)==2913 and paired[['Baseline','Tmax+2']].notna().all().all()
paired['warming_difference']=paired['Tmax+2']-paired['Baseline']
write(seasons,'paired_season_endpoints')
write(paired,'paired_season_differences')

# 2. 每年先按站点等权平均，再累加12年均值；年度SD只保留在表中。
annual=seasons.groupby(['region','mode','sowing_year','climate']).seasonal_N2O.agg(['mean','std','count']).reset_index()
annual=annual.rename(columns={'mean':'annual_site_mean','std':'annual_site_sd','count':'n_sites'})
assert len(annual)==12*3*12*2
assert annual.groupby(['region','mode','sowing_year']).n_sites.nunique().eq(1).all()
totals=annual.groupby(['region','mode','climate']).annual_site_mean.sum().unstack('climate').reset_index()
totals['difference']=totals['Tmax+2']-totals['Baseline']
totals['percent_change']=np.where(totals.Baseline.ne(0),100*totals.difference/totals.Baseline,np.nan)
write(annual,'regional_annual_site_mean_sd')
write(totals,'regional_12year_sum')



assert all(sha(Path(p)) == h for p,h in hashes.items())
print('Verified:', len(paired), 'pairs;', len(annual), 'annual rows;', len(totals), 'regional totals')

