"""Regional summer aggregation, trimmed from original scripts 01 and 04.

Uses the preserved two-column spatial crosswalk; writes only 615 input rows.
NH: June-August. SH: previous December-February, labeled by end year.
"""
import argparse
import calendar
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol

ROOT = Path(__file__).resolve().parent

def months(year, hemisphere):
    if hemisphere == "NH":
        return [(year, 6), (year, 7), (year, 8)]
    if hemisphere == "SH":
        return [(year - 1, 12), (year, 1), (year, 2)]
    raise ValueError("Unexpected hemisphere")


def seasonal_stats(values, weights):
    # 三个月均有效才参与加权；真实零值保留，缺失月份不补零。
    values, weights = np.asarray(values, dtype=float), np.asarray(weights, dtype=float)
    if values.shape != (3, len(weights)) or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("Invalid season/area arrays")
    if np.isinf(values).any() or (values[np.isfinite(values)] < 0).any():
        raise ValueError("Invalid EDH30 values")
    total = values.sum(axis=0)
    valid = np.isfinite(total)
    area = weights[valid].sum()
    if area <= 0:
        raise ValueError("No complete-season valid area")
    numerator = np.sum(total[valid] * weights[valid])
    return dict(area_weighted_mean=float(numerator / area), weighted_sum=float(numerator),
                total_maize_area=float(weights.sum()), valid_area=float(area),
                missing_temperature_area=float(weights[~valid].sum()),
                valid_area_share=float(area / weights.sum()), maize_cell_count=len(weights),
                valid_cell_count=int(valid.sum()), zero_season_cell_count=int((total == 0).sum()))



def read_samples(tif_dir, year, month, lon, lat):
    path = tif_dir / f"era5land_heat_monthly_{year}.tif"
    with rasterio.open(path) as src:
        if src.crs != rasterio.crs.CRS.from_epsg(4326) or src.shape != (1800, 3600):
            raise ValueError(f"Unexpected raster grid: {path}")
        with rasterio.open(tif_dir / 'era5land_heat_monthly_1985.tif') as reference:
            np.testing.assert_allclose(tuple(src.transform), tuple(reference.transform), rtol=0, atol=1e-12)
        np.testing.assert_allclose([src.transform.a, src.transform.e], [.1, -.1], rtol=0, atol=1e-6)
        name = f"EDH30_{year}_{month:02d}"
        if src.descriptions.count(name) != 1:
            raise ValueError(f"Missing or duplicate band {name}")
        rr, cc = map(np.asarray, rowcol(src.transform, lon, lat))
        inside = (rr >= 0) & (rr < src.height) & (cc >= 0) & (cc < src.width)
        values = np.full(len(lon), np.nan)
        arr = src.read(src.descriptions.index(name) + 1, masked=True).filled(np.nan)
        values[inside] = arr[rr[inside], cc[inside]]
        return values



def compute_year(year, area_dir, tif_dir, crosswalk, regions):
    # 固定区域归属，面积取季节结束年；2021 年以后复用 2020 年。
    area = pd.read_csv(area_dir / f'global_maize_area_{min(year, 2020)}.csv.gz',
                       usecols=['maize_id', 'maize_lon', 'maize_lat', 'maize_area'])
    if area.maize_id.duplicated().any() or not np.isfinite(area.maize_area).all() or (area.maize_area <= 0).any():
        raise ValueError('Invalid positive-area table')
    area = area.merge(crosswalk, on='maize_id', how='left', validate='one_to_one', indicator=True)
    if not area['_merge'].eq('both').all():
        raise ValueError('Annual grid IDs missing from crosswalk')
    rows = []
    for hemi in ['NH', 'SH']:
        ids = regions.loc[regions.hemisphere.eq(hemi), 'region_id'].to_list()
        subset = area.loc[area.region_id.isin(ids)].reset_index(drop=True)
        ms = months(year, hemi)
        samples = np.vstack([read_samples(tif_dir, y, m, subset.maize_lon, subset.maize_lat) for y, m in ms])
        for rid in ids:
            region = regions.loc[regions.region_id.eq(rid)].iloc[0]
            chosen = subset.region_id.eq(rid).to_numpy()
            record = dict(region_id=rid, region_name=region.region_name, hemisphere=hemi,
                          season_year=year, area_source_year=min(year, 2020),
                          start_month=f'{ms[0][0]}-{ms[0][1]:02d}', end_month=f'{ms[-1][0]}-{ms[-1][1]:02d}',
                          season_days=sum(calendar.monthrange(y, m)[1] for y, m in ms))
            record.update(seasonal_stats(samples[:, chosen], subset.loc[chosen, 'maize_area'].to_numpy()))
            rows.append(record)
    return pd.DataFrame(rows)


def main():
    # 1. 检查输入范围和输出路径，不覆盖已有核心数据。
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tif-dir', required=True, type=Path)
    parser.add_argument('--area-dir', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Output exists; use a new output path.')
    for year in range(1984, 2026):
        path = args.tif_dir / f'era5land_heat_monthly_{year}.tif'
        if not path.is_file():
            raise FileNotFoundError(path)
    for year in range(1985, 2021):
        path = args.area_dir / f'global_maize_area_{year}.csv.gz'
        if not path.is_file():
            raise FileNotFoundError(path)
    crosswalk = pd.read_csv(ROOT / 'maize_region_crosswalk.csv.gz', keep_default_na=False)
    regions = pd.read_csv(ROOT / 'regions.csv').sort_values('order')
    if len(regions) != 15 or regions.region_id.duplicated().any():
        raise ValueError('Expected 15 unique regions')
    if not set(crosswalk.region_id).issubset(set(regions.region_id) | {''}):
        raise ValueError('Unexpected spatial assignment')
    # 2. 按原版季节顺序计算 15 区域的 41 年面积加权暴露。
    tables = []
    for year in range(1985, 2026):
        tables.append(compute_year(year, args.area_dir, args.tif_dir, crosswalk, regions))
        print(f'{year}: 15 summer estimates', flush=True)
    result = pd.concat(tables, ignore_index=True).sort_values(['region_id', 'season_year'])
    # 3. 核对唯一性、加权均值和有效面积后，只写一份输入表。
    assert len(result) == 615 and not result.duplicated(['region_id', 'season_year']).any()
    np.testing.assert_allclose(result.weighted_sum/result.valid_area, result.area_weighted_mean, rtol=1e-12)
    assert result.valid_area_share.between(0, 1, inclusive='right').all()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
