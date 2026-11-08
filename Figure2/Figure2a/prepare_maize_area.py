"""Rebuild the 2011-2020 complete-decade maize harvested-area raster."""
import argparse
from pathlib import Path

import h5py
import numpy as np
import rasterio
from rasterio.transform import from_origin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--area-h5", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # 1. 累计十年有效面积；原始零面积写为 NaN，不能填零。
    total = np.zeros((4320, 2160), dtype=np.float64)
    counts = np.zeros(total.shape, dtype=np.uint8)
    with h5py.File(args.area_h5, "r") as handle:
        ds = handle["Maize"]
        assert ds.shape == (60, 4320, 2160)
        for year in range(2011, 2021):
            values = ds[year - 1961]
            valid = np.isfinite(values)
            assert not np.isinf(values).any() and not (values[valid] < 0).any()
            total[valid] += values[valid]
            counts += valid
    # 2. 仅保留十年均有效像元；由经纬度数组转为栅格行列。
    complete = counts == 10
    mean = np.full(total.shape, np.nan, dtype=np.float32)
    mean[complete] = total[complete] / 10
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        args.output,
        "w",
        driver="GTiff",
        width=4320,
        height=2160,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=from_origin(-180, 90, 1 / 12, 1 / 12),
        nodata=np.nan,
        compress="deflate",
    ) as dst:
        dst.write(mean.T, 1)
    print(f"Saved {args.output}; {complete.sum()} complete-decade cells.")


if __name__ == "__main__":
    main()
