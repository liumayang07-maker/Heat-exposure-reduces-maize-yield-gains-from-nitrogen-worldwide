# 分析入口：沿用原计算方法，直接运行 Python。
import json, hashlib, platform
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
HOME = ROOT / 'APSIM_outputs'
OUT = ROOT / 'N2O_analysis'
# 1. 小时源数据位于归档外，运行时显式指定目录。
parser = argparse.ArgumentParser(description='Rebuild baseline N2O and EDH30 bins.')
parser.add_argument(
    '--hourly-dir',
    type=Path,
    required=True,
    help='Directory containing <xSite>_2010_2025_hourly_utc.csv files.',
)
HOURLY = parser.parse_args().hourly_dir
if not HOURLY.is_dir():
    parser.error(f'Hourly directory does not exist: {HOURLY}')
FIRST_YEAR = 2011
LAST_YEAR = 2022
N_BINS = 100
MODES = ['F1_mzI0', 'F1_mzI2', 'F2_mzI2']
LOSSES = ['mzN2Oatm']
STAGES = ['end_of_juvenil', 'flag_leaf']
SOURCE = HOME / 'baseline_daily.csv.gz'
OUT.mkdir(parents=True, exist_ok=True)

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def write(d, name):
    d.to_csv(OUT / (name + '.csv'), index=False, encoding='utf-8-sig')

def thermal(t):
    g = np.clip(t-10, 0, 20)/24
    e = np.maximum(t-30, 0)
    if not np.isfinite(t).all() or g.sum() <= 0:
        raise ValueError('Invalid temperature or zero seasonal GDD')
    relative = np.minimum(np.cumsum(g)/g.sum(), 1)
    bins = np.clip(np.floor(relative*N_BINS).astype(int), 0, N_BINS-1)
    return g, e, relative, bins

def allocate(cumulative, day_index, bins, n_days):
    # 第0天累计量视为0；每天新增量在该日全部小时均匀分配。
    increments = np.diff(cumulative, prepend=0)
    if not np.isfinite(increments).all() or np.any(increments < -1e-8):
        raise ValueError('Missing or decreasing cumulative loss')
    counts = np.bincount(day_index, minlength=n_days)
    if np.any(counts == 0):
        raise ValueError('Day without temperature hours')
    result = np.bincount(bins, weights=increments[day_index]/counts[day_index], minlength=N_BINS)
    if not np.isclose(result.sum(), cumulative[-1], atol=1e-8):
        raise ValueError('Loss conservation failed')
    return result

# 3. 小型算法测试覆盖跨日分配、首日累计量、零高温及非法累计量。
di = np.repeat(np.arange(3), 24)
bi = np.r_[np.zeros(60, dtype=int), np.ones(12, dtype=int)]
test = allocate(np.array([.1,.3,.6]), di, bi, 3)
assert np.allclose(test[:2], [.45,.15])
assert thermal(np.full(24, 20.))[1].sum() == 0
assert np.isclose(np.std([1,3], ddof=1), np.sqrt(2))
for bad in [np.array([.1,.05,.6]), np.array([.1,np.nan,.6])]:
    try:
        allocate(bad, di, bi, 3)
    except ValueError:
        pass
    else:
        raise AssertionError('Invalid increments accepted')

if True:
    paths = [SOURCE, ROOT/'APSIM_inputs/site_soil_configuration.csv', ROOT/'N2O_analysis/site_metadata.csv']
    hashes = {str(p):sha(p) for p in paths}
    meta = pd.read_csv(paths[2]).set_index('xSite')
    sow = pd.to_datetime(pd.read_csv(paths[1]).set_index('selected_site_id')['SowDate'], format='%d-%b')
    tz = meta['timezone_name']
    columns = ['xSite','xIW','mzNrate','mzManagement','Date()','Maize.stage','Maize.status','mzN2Oatm']
    data = pd.read_csv(SOURCE, usecols=columns)
    data['source_row'] = np.arange(len(data))+2
    data['date'] = pd.to_datetime(data['Date()'], format='%d/%m/%Y')
    assert set(data.mzManagement) == set(MODES)
    bins_out, seasons, stages, issues = [], [], [], []
    conservation = []
    for site, sd in data.groupby('xSite', sort=True):
        weather_path = HOURLY/f'{site}_2010_2025_hourly_utc.csv'
        hashes[str(weather_path)] = sha(weather_path)
        w = pd.read_csv(weather_path, usecols=['time_utc','t2m_C'])
        w.index = pd.to_datetime(w.pop('time_utc'), utc=True).dt.tz_convert(tz[site])
        assert w.index.is_unique and w.index.is_monotonic_increasing
        md = pd.Timestamp(sow[site])
        for mode, d in sd.groupby('mzManagement'):
            d = d.sort_values('date').copy()
            starts = d['Maize.stage'].eq('sowing') & ~d['Maize.stage'].shift().eq('sowing')
            d['season_index'] = starts.cumsum()
            for _, s in d.groupby('season_index'):
                first = s.iloc[0]['date']
                if s.iloc[0]['Maize.stage'] != 'sowing':
                    continue
                configured = pd.Timestamp(first.year, md.month, md.day)
                if configured >= first:
                    configured = configured.replace(year=configured.year-1)
                year = configured.year
                if not FIRST_YEAR <= year <= LAST_YEAR:
                    continue
                key = dict(xSite=site, region=meta.loc[site,'figure2_region_code'], mode=mode, sowing_year=year)
                try:
                    assert (first-configured).days == 1, 'First record is not sowing day + 1'
                    ripe = s[s['Maize.stage'].eq('harvest_ripe')]
                    assert len(ripe) == 1, 'Harvest boundary absent or ambiguous'
                    s = s[s.date <= ripe.iloc[0].date].copy()
                    dates = pd.DatetimeIndex(s.date)
                    assert dates.is_unique and dates.equals(pd.date_range(first, dates[-1])), 'Daily dates incomplete'
                    start = first.tz_localize(tz[site], nonexistent='shift_forward', ambiguous=True)
                    end = (dates[-1]+pd.Timedelta(days=1)).tz_localize(tz[site], nonexistent='shift_forward', ambiguous=True)
                    h = w[(w.index >= start) & (w.index < end)]
                    # 半小时时区的当地记录在半点；完整性仍按UTC整点核查。
                    expected = pd.date_range(start.tz_convert('UTC').ceil('h'), end.tz_convert('UTC'), freq='h', inclusive='left').tz_convert(tz[site])
                    assert h.index.equals(expected), 'Hourly calendar coverage incomplete'
                    day_index = dates.get_indexer(h.index.tz_localize(None).normalize())
                    assert np.all(day_index >= 0)
                    g, e, progress, b = thermal(h.t2m_C.to_numpy())
                    hours = np.bincount(b, minlength=N_BINS)
                    heat = np.bincount(b, weights=e, minlength=N_BINS)
                    frame = pd.DataFrame({'bin':np.arange(1,101), 'hours':hours, 'EDH30_increment':heat,
                        'EDH30_cumulative':heat.cumsum(), 'EDH30_intensity':np.divide(heat,hours,out=np.full(100,np.nan),where=hours>0)})
                    for k,v in key.items(): frame[k] = v
                    for loss in LOSSES:
                        values = allocate(s[loss].to_numpy(), day_index, b, len(s))
                        frame[loss+'_increment'] = values
                        frame[loss+'_cumulative'] = values.cumsum()
                        conservation.append(abs(values.sum()-float(s[loss].iloc[-1])))
                    assert np.isclose(heat.sum(), e.sum()) and hours.sum() == len(h)
                    bins_out.append(frame)
                    seasons.append({**key, 'sowing_date':configured, 'first_record_date':first, 'harvest_date':dates[-1],
                        'source_first_row':int(s.source_row.iloc[0]), 'source_last_row':int(s.source_row.iloc[-1]),
                        'days':len(s), 'hours':len(h), 'GDD':g.sum(), 'EDH30':e.sum(),
                        'dead_days':int(s['Maize.status'].eq('dead').sum()),
                        **{v:float(s[v].iloc[-1]) for v in LOSSES},
                        **{v+'_first_day':float(s[v].iloc[0]) for v in LOSSES}})
                    for stage in STAGES:
                        hit = s[s['Maize.stage'].eq(stage)]
                        if len(hit):
                            day = hit.iloc[0].date
                            pct = progress[day_index == dates.get_loc(day)][-1]*100
                            stages.append({**key,'event':stage,'date':day,'progress':pct})
                except (AssertionError, ValueError) as exc:
                    issues.append({**key,'first_record_date':first,'source_first_row':int(s.source_row.iloc[0]),'reason':str(exc)})
        print(f'Processed {site}: {len(seasons)} seasons', flush=True)
    season = pd.DataFrame(seasons)
    table = pd.concat(bins_out, ignore_index=True)
    assert not season.duplicated(['xSite','mode','sowing_year']).any()
    assert not table.duplicated(['xSite','mode','sowing_year','bin']).any()
    write(season,'season_baseline')
    write(table,'season_bins')
    write(pd.DataFrame(stages),'stage_positions')
    metrics = [c for c in table if c.endswith(('_increment','_cumulative','_intensity'))]
    summary = table.groupby(['region','mode','bin'])[metrics].agg(['mean','std','count'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary = summary.join(table.groupby(['region','mode','bin']).xSite.nunique().rename('n_sites')).reset_index()
    write(summary,'regional_bins')
    coverage = season.groupby(['region','mode']).agg(n_sites=('xSite','nunique'),n_site_years=('sowing_year','size')).reset_index()
    write(coverage,'coverage')

    assert len(season) == 2916 and len(table) == 291600 and not issues
    verify = dict(status='passed', n_sites=int(season.xSite.nunique()),
                  n_seasons=len(season), n_bins=len(table),
                  max_conservation_error=max(conservation),
                  input_hashes=hashes, python=platform.python_version(),
                  numpy=np.__version__, pandas=pd.__version__)
    assert all(sha(Path(p)) == h for p, h in hashes.items())
    (OUT/'baseline_verification.json').write_text(json.dumps(verify, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in verify.items() if k != 'input_hashes'}))

