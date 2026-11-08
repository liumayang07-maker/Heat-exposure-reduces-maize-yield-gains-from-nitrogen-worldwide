"""Rebuild Figure1b's 240 annual means, with the original area weighting.

Trimmed from Figure1_v3/scripts/plot_continent_annual_EDH30_moving_average.py,
function compute(). No moving averages, plots, logs or redundant global series.
The supplied crosswalk retains the original continent assignments exactly.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol

CONTINENTS = ('North America', 'South America', 'Europe', 'Asia', 'Africa', 'Oceania')


def compute(tif_dir, area_dir, cross_path):
    cross = pd.read_csv(cross_path, usecols=['maize_id', 'continent'], keep_default_na=False)
    if cross.maize_id.duplicated().any():
        raise ValueError('Duplicate crosswalk IDs')
    rows = []
    grid = None
    for year in range(1986, 2026):
        area_year = min(year, 2020)
        area = pd.read_csv(area_dir / f'global_maize_area_{area_year}.csv.gz',
                           usecols=['maize_id', 'maize_lon', 'maize_lat', 'maize_area'])
        area = area.loc[np.isfinite(area.maize_area) & (area.maize_area > 0)]
        area = area.merge(cross, on='maize_id', how='left', validate='one_to_one', indicator=True)
        if not area['_merge'].eq('both').all():
            raise ValueError(f'Missing crosswalk IDs: {year}')
        if not np.isfinite(area[['maize_lon', 'maize_lat']]).all().all():
            raise ValueError('Nonfinite coordinates')
        with rasterio.open(tif_dir / f'era5land_heat_monthly_{year}.tif') as src:
            if src.crs != rasterio.crs.CRS.from_epsg(4326):
                raise ValueError('Expected EPSG:4326')
            current = (src.width, src.height, src.transform, src.crs)
            if grid is not None and current != grid:
                raise ValueError(f'Climate grid changed in {year}')
            grid = current
            rr, cc = map(np.asarray, rowcol(src.transform, area.maize_lon, area.maize_lat))
            inside = (rr >= 0) & (rr < src.height) & (cc >= 0) & (cc < src.width)
            # Preserve the source continent pipeline's float64 monthly accumulation.
            annual = np.zeros(len(area), dtype=np.float64)
            annual[~inside] = np.nan
            indices = []
            for month in range(1, 13):
                name = f'EDH30_{year}_{month:02d}'
                if src.descriptions.count(name) != 1:
                    raise ValueError(f'Missing or duplicate band: {name}')
                indices.append(src.descriptions.index(name) + 1)
            # Read together to avoid decompressing interleaved TIFF tiles 12 times.
            bands = src.read(indices, masked=True).filled(np.nan)
            for band in bands:
                values = band[rr[inside], cc[inside]]
                if np.any(values < 0) or np.isinf(values).any():
                    raise ValueError('Invalid EDH30')
                annual[inside] += values
        weights = area.maize_area.to_numpy(dtype=float)
        valid = np.isfinite(annual)
        for continent in CONTINENTS:
            member = area.continent.eq(continent).to_numpy()
            use = member & valid
            total, valid_area = weights[member].sum(), weights[use].sum()
            if valid_area <= 0:
                raise ValueError(f'No valid area: {continent}, {year}')
            rows.append(dict(continent=continent, year=year, area_source_year=area_year,
                             total_maize_area=total, valid_area=valid_area,
                             valid_area_share=valid_area/total, valid_cells=int(use.sum()),
                             weighted_sum=np.dot(annual[use], weights[use]),
                             area_weighted_mean=np.dot(annual[use], weights[use])/valid_area))
        print(f'{year}: six continents complete', flush=True)
    return pd.DataFrame(rows).sort_values(['continent', 'year'])


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
    table = compute(args.tif_dir, args.area_dir,
                    Path(__file__).resolve().with_name('maize_continent_crosswalk.csv.gz'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    print(f'Created {args.output}: {len(table)} rows', flush=True)


if __name__ == '__main__':
    main()
