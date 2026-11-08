"""Extract the annual maize area inputs needed by Figure1a/b.

Trimmed from global_maize_edh_1980_2025/scripts/python/
01_prepare_global_maize_area.py and edh_config.py. Preserves HDF5 indexing,
grid centers, positive-area selection and maize IDs. No union grid or logs.
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def write_year_area(dset, year, out_path):
    # Use multiplication as in edh_config.py to preserve original coordinates.
    lons = -180.0 + (np.arange(4320, dtype=np.float64) + 0.5) * (1.0 / 12.0)
    lats = 90.0 - (np.arange(2160, dtype=np.float64) + 0.5) * (1.0 / 12.0)
    header = True
    positive_cells = 0
    for lon0 in range(0, dset.shape[1], 80):
        lon1 = min(dset.shape[1], lon0 + 80)
        block = dset[year - 1961, lon0:lon1, :].astype(np.float64)
        positive = np.isfinite(block) & (block > 0)
        if not positive.any():
            continue
        lon_pos, lat_idx = np.where(positive)
        lon_idx = lon0 + lon_pos
        df = pd.DataFrame({
            'maize_id': lon_idx.astype(np.int64) * 2160 + lat_idx.astype(np.int64) + 1,
            'maize_lon': lons[lon_idx], 'maize_lat': lats[lat_idx],
            'maize_area': block[positive].astype(np.float64)})
        df.to_csv(out_path, index=False, mode='w' if header else 'a',
                  header=header, compression='gzip')
        header = False
        positive_cells += len(df)
    if header:
        raise ValueError(f'No positive maize area: {year}')
    print(f'{year}: {positive_cells} positive maize cells', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--area-h5', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--start-year', type=int, default=1986, choices=(1985, 1986),
                        help='Use 1985 to include Figure1d heatmap inputs.')
    args = parser.parse_args()
    paths = [args.output_dir / f'global_maize_area_{y}.csv.gz' for y in range(args.start_year, 2021)]
    if any(p.exists() for p in paths):
        raise FileExistsError('Use a new output directory; existing annual tables are not overwritten.')
    with h5py.File(args.area_h5, 'r') as handle:
        dset = handle['Maize']
        if dset.shape != (60, 4320, 2160):
            raise ValueError(f'Unexpected Maize dimensions: {dset.shape}')
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for year, path in zip(range(args.start_year, 2021), paths):
            write_year_area(dset, year, path)


if __name__ == '__main__':
    main()
