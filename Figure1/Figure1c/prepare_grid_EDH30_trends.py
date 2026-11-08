"""Rebuild Figure1c from monthly exports, or original maize-masked annual TIFFs.

Trimmed from original scripts 17 (OLS slopes) and 18 (HAC and global BH).
Monthly sums and maize masks reuse Figure1a's original-derived functions.
No annual caches, diagnostic reports, or alternative figures are written.
"""
import argparse
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

LAGS = 3
ALPHA = .05
YEARS = range(1986, 2026)
OUTPUTS = ('EDH30_slope_per_decade.tif', 'EDH30_HAC_BH_adjusted_pvalue.tif',
           'Figure1c_Source_Data.csv.gz')


def hac_test(y):
    n = y.shape[0]
    x = np.arange(n, dtype=float) - (n-1)/2
    sxx = x @ x
    slope = (x[:, None]*y).sum(axis=0)/sxx
    residual = y-y.mean(axis=0)-x[:, None]*slope
    # The slope row of (X'X)^(-1)X' is x / sum(x^2) for centered time.
    score = residual*x[:, None]/sxx
    variance = (score*score).sum(axis=0)
    for lag in range(1, LAGS+1):
        variance += 2*(1-lag/(LAGS+1))*(score[lag:]*score[:-lag]).sum(axis=0)
    se = np.sqrt(np.maximum(variance, 0)*n/(n-2))
    statistic = np.divide(np.abs(slope), se, out=np.full_like(slope,np.inf), where=se>0)
    p = 2*t.sf(statistic, df=n-2)
    constant = np.ptp(y, axis=0)==0
    slope[constant], se[constant], p[constant] = 0, 0, 1
    return slope*10, se*10, p


def read_annual(args):
    if args.tif_dir:
        sys.dont_write_bytecode = True
        path = Path(__file__).resolve().parents[1] / 'Figure1a' / '03_prepare_mean_EDH30.py'
        spec = importlib.util.spec_from_file_location('figure1a_upstream', path)
        shared = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(shared)
    for year in YEARS:
        if args.tif_dir:
            annual, profile = shared.annual_edh30(args.tif_dir / f'era5land_heat_monthly_{year}.tif', year)
            annual[~shared.build_year_maize_mask(year, profile, args.area_dir)] = np.nan
        else:
            path = args.annual_dir / f'gee_era5land_annual_heat_maize_area_{year}.tif'
            with rasterio.open(path) as source:
                band = f'annual_EDH30_{year}'
                if source.descriptions.count(band) != 1 or source.tags().get('maize_mask_applied') != 'True':
                    raise ValueError(f'Annual band or maize-mask provenance missing: {path}')
                annual = source.read(source.descriptions.index(band) + 1, masked=True).filled(np.nan)
                profile = source.profile.copy()
        if profile['crs'] != rasterio.crs.CRS.from_epsg(4326):
            raise ValueError(f'Expected EPSG:4326 in {year}')
        if np.isinf(annual).any() or np.any(annual < 0):
            raise ValueError(f'Invalid EDH30 in {year}')
        yield year, annual, profile


def compute(args):
    profile = None
    for k, (year, annual, current) in enumerate(read_annual(args)):
        if profile is None:
            profile = current
            indices = np.flatnonzero(np.isfinite(annual.ravel()))
            y = np.empty((40, len(indices)), dtype='float64')
        if any(profile[key] != current[key] for key in ('width', 'height', 'transform', 'crs')):
            raise ValueError(f'Grid changed in {year}')
        y[k] = annual.ravel()[indices]
        print(f'Read {year}', flush=True)
    complete = np.isfinite(y).all(axis=0)
    indices, y = indices[complete], y[:, complete]
    if not len(indices):
        raise ValueError('No complete 40-year maize records.')
    # Preserve script 17's chronological sufficient-statistic slope calculation.
    sy, sxy = np.zeros(len(indices)), np.zeros(len(indices))
    for k in range(40):
        sy += y[k]
        sxy += k * y[k]
    slope_raster_values = (10 * (40 * sxy - 780 * sy) / (40 * 20540 - 780**2)).astype('float32')
    slope, se, p = [np.empty(len(indices)) for _ in range(3)]
    for start in range(0, len(indices), 20000):
        stop = min(start + 20000, len(indices))
        slope[start:stop], se[start:stop], p[start:stop] = hac_test(y[:, start:stop])
    np.testing.assert_allclose(slope, slope_raster_values, rtol=1e-6, atol=1e-5)
    _, adjusted, _, _ = multipletests(p, alpha=ALPHA, method='fdr_bh')
    rows, columns = np.unravel_index(indices, (profile['height'], profile['width']))
    longitude, latitude = rasterio.transform.xy(profile['transform'], rows, columns)
    table = pd.DataFrame(dict(
        longitude=longitude, latitude=latitude, slope_per_decade=slope,
        hac_se_per_decade=se, ci_low_per_decade=slope-t.ppf(.975, 38)*se,
        ci_high_per_decade=slope+t.ppf(.975, 38)*se, p_value=p, p_BH_all_grids=adjusted))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, values, dtype in [(OUTPUTS[0], slope_raster_values, 'float32'),
                                (OUTPUTS[1], adjusted, 'float64')]:
        raster = np.full((profile['height'], profile['width']), np.nan, dtype=dtype)
        raster.ravel()[indices] = values
        output_profile = profile.copy()
        output_profile.update(count=1, dtype=dtype, nodata=np.nan, compress='deflate', predictor=3)
        with rasterio.open(args.output_dir / name, 'w', **output_profile) as target:
            target.write(raster, 1)
            target.set_band_description(1, Path(name).stem)
            target.update_tags(start_year='1986', end_year='2025', min_years='40',
                               estimator='OLS with intercept; slope multiplied by 10',
                               test='HAC Bartlett lags=3; finite correction 40/38; two-sided t(38)',
                               correction='BH over all complete grid cells including constants')
    table.to_csv(args.output_dir / OUTPUTS[2], index=False, compression='gzip')
    print(f'Complete cells: {len(table)}; BH p < 0.05: {(adjusted < ALPHA).sum()}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--tif-dir', type=Path)
    source.add_argument('--annual-dir', type=Path)
    parser.add_argument('--area-dir', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.tif_dir and not args.area_dir:
        parser.error('--area-dir is required with --tif-dir')
    if args.annual_dir and args.area_dir:
        parser.error('--area-dir is not used with already masked annual TIFFs')
    if any((args.output_dir / name).exists() for name in OUTPUTS):
        raise FileExistsError('Outputs exist; use a new output directory.')
    for year in YEARS:
        paths = ([args.tif_dir / f'era5land_heat_monthly_{year}.tif',
                  args.area_dir / f'global_maize_area_{min(year, 2020)}.csv.gz'] if args.tif_dir else
                 [args.annual_dir / f'gee_era5land_annual_heat_maize_area_{year}.tif'])
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError(path)
    compute(args)


if __name__ == '__main__':
    main()
