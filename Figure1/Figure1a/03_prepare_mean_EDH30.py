"""Rebuild Figure1a's mean from monthly GEE exports and positive maize areas.

Trimmed from global_maize_edh_1980_2025/scripts/python/
11_plot_multiyear_mean_global_heat_maps.py: yearly_annual_metrics,
build_year_maize_mask, compute_multiyear_mean. No other metrics, caches or plots.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol


def build_year_maize_mask(year, profile, maize_area_dir):
    mask = np.zeros((profile['height'], profile['width']), dtype=bool)
    path = maize_area_dir / f'global_maize_area_{min(year, 2020)}.csv.gz'
    df = pd.read_csv(path, usecols=['maize_lon', 'maize_lat', 'maize_area'])
    df = df.loc[df.maize_area > 0]
    rr, cc = map(np.asarray, rowcol(profile['transform'], df.maize_lon, df.maize_lat))
    inside = (rr >= 0) & (rr < profile['height']) & (cc >= 0) & (cc < profile['width'])
    mask[rr[inside], cc[inside]] = True
    return mask


def annual_edh30(path, year):
    with rasterio.open(path) as src:
        if src.crs != rasterio.crs.CRS.from_epsg(4326):
            raise ValueError(f'Expected EPSG:4326: {path}')
        total = None
        for month in range(1, 13):
            name = f'EDH30_{year}_{month:02d}'
            if src.descriptions.count(name) != 1:
                raise ValueError(f'Missing or duplicate band {name}: {path}')
            arr = src.read(src.descriptions.index(name) + 1, masked=True).astype('float32').filled(np.nan)
            if np.any(arr < 0) or np.isinf(arr).any():
                raise ValueError(f'Invalid EDH30 values: {name}')
            # Original map pipeline sums months in float32, then years in float64.
            total = arr if total is None else total + arr
        return total, src.profile.copy()


def compute_multiyear_mean(tif_dir, area_dir):
    profile = None
    for year in range(1986, 2026):
        annual, current = annual_edh30(tif_dir / f'era5land_heat_monthly_{year}.tif', year)
        if profile is None:
            profile = current
            accum = np.zeros(annual.shape, dtype='float64')
            count = np.zeros(annual.shape, dtype='uint16')
        if any(current[k] != profile[k] for k in ('width', 'height', 'transform', 'crs')):
            raise ValueError(f'Climate grid changed in {year}')
        mask = build_year_maize_mask(year, current, area_dir)
        annual[~mask] = np.nan
        finite = np.isfinite(annual)
        accum[finite] += annual[finite].astype('float64')
        count[finite] += 1
        print(f'{year}: {finite.sum()} valid masked pixels', flush=True)
    mean = np.full(accum.shape, np.nan, dtype='float32')
    valid = count > 0
    mean[valid] = (accum[valid] / count[valid]).astype('float32')
    return mean, profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tif-dir', required=True, type=Path)
    parser.add_argument('--area-dir', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Output exists; use a new output path.')
    for year in range(1986, 2026):
        for p in (args.tif_dir / f'era5land_heat_monthly_{year}.tif',
                  args.area_dir / f'global_maize_area_{min(year, 2020)}.csv.gz'):
            if not p.is_file():
                raise FileNotFoundError(p)
    mean, profile = compute_multiyear_mean(args.tif_dir, args.area_dir)
    profile.update(count=1, dtype='float32', compress='deflate', predictor=3,
                   tiled=True, blockxsize=256, blockysize=256, nodata=np.nan)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(args.output, 'w', **profile) as dst:
        dst.write(mean, 1)
        dst.set_band_description(1, 'mean_annual_EDH30')
    print(f'Created {args.output}', flush=True)


if __name__ == '__main__':
    main()
