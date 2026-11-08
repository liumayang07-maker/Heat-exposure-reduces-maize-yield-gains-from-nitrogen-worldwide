"""Verify/rebuild regional median parameters and the shared 0.25-degree soil."""
import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import Resampling, reproject

from prepare_inputs import HERE, SHARED, STAGES, FRACTIONS, read_table


def regional_parameters():
    # 1. 复用 RF 原始输入，每个环境只计一次；不重复保存整张主表。
    source = read_table(HERE.parent / "RF_model/data/training_source.csv")
    columns = [
        "env_ID", "lat", "lon", "GDD_total", "EDH30_total", "total_rainfall",
        "annual_DD10_30_C_d", "annual_EDH30_C_h", "annual_precip_mm",
    ]
    assert source.groupby("env_ID")[columns].nunique(dropna=False).le(1).all().all()
    env = source[columns].drop_duplicates("env_ID")
    env = env.merge(read_table(SHARED / "environment_timing_fractions.csv"), on="env_ID", validate="one_to_one")
    polygons = gpd.read_file(HERE.parents[1] / "Figure2/Figure2a/study_regions.geojson").to_crs(4326)
    points = gpd.GeoDataFrame(env, geometry=gpd.points_from_xy(env.lon, env.lat), crs=4326)
    joined = gpd.sjoin(points, polygons[["region_id", "geometry"]], how="left", predicate="intersects")
    assert not joined.env_ID.duplicated().any()
    lookup = read_table(SHARED / "regional_parameters_historical.csv")[["region_code", "region_id", "region_name"]]
    rows = []
    for region in lookup.itertuples(index=False):
        current = joined.loc[joined.region_id.eq(region.region_id)]
        heat = current.loc[current.EDH30_total.gt(0)]
        medians = heat[STAGES].median()
        medians /= medians.sum()
        annual_positive = current.loc[current.annual_EDH30_C_h.gt(0)]
        rows.append(dict(
            region_code=region.region_code,
            region_id=region.region_id,
            region_name=region.region_name,
            growing_season_GDD_fraction=float((current.GDD_total / current.annual_DD10_30_C_d).median()),
            growing_season_precip_fraction=float((current.total_rainfall / current.annual_precip_mm).median()),
            growing_season_EDH30_fraction=float((annual_positive.EDH30_total / annual_positive.annual_EDH30_C_h).quantile(.5, interpolation="linear")),
            **medians.to_dict(),
        ))
    parameters = pd.DataFrame(rows).sort_values("region_code")
    parameters.loc[parameters.growing_season_EDH30_fraction.eq(0), STAGES] = 0.
    density = read_table(SHARED / "density_settings.csv")
    result = {}
    for kind, dataset in [
        ("historical", "historical_2011_2022"),
        ("future", "SSP245_2030_2050_fixed2020"),
    ]:
        frame = parameters.copy()
        settings = density.loc[density.dataset.eq(dataset)].set_index("region_code")
        frame["plant_density"] = frame.region_code.map(settings.plant_density)
        expected = read_table(SHARED / f"regional_parameters_{kind}.csv")
        np.testing.assert_array_equal(frame.region_code, expected.region_code)
        np.testing.assert_allclose(frame[[*FRACTIONS, "plant_density"]], expected[[*FRACTIONS, "plant_density"]], atol=1e-12, rtol=1e-12)
        result[kind] = frame[expected.columns]
    return result


def soil_layer():
    # 2. 原版有效像元平均重采样，不插值补缺、不再次换算 SOC 单位。
    with rasterio.open(SHARED / "static_soil_0p25.tif") as ref:
        shape, profile = ref.shape, ref.profile.copy()
        target_transform, target_crs = ref.transform, ref.crs
        expected = ref.read()
    with rasterio.open(SOURCE_ROOT / "soil_0p1.tif") as src:
        assert src.descriptions[:3] == ("SOC_g_kg", "soil_total_N_g_kg", "soil_pH_water")
        raw = src.read([1, 2, 3], masked=True).filled(np.nan)
        valid = (src.read(4) == 1) & np.isfinite(raw).all(axis=0)
        kwargs = dict(
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=target_transform,
            dst_crs=target_crs,
            resampling=Resampling.average,
        )
        values = np.full((4, *shape), np.nan, dtype="float32")
        for band in range(3):
            reproject(np.where(valid, raw[band], np.nan), values[band], src_nodata=np.nan, dst_nodata=np.nan, **kwargs)
        values[3] = 0
        reproject(valid.astype("float32"), values[3], **kwargs)
    target = np.zeros(shape, dtype=bool)
    target.ravel()[read_table(SHARED / "future_target_grid.csv").pixel_id.to_numpy(int)] = True
    values = np.where(target[None, :, :] & np.isfinite(values), values, -9999.).astype("float32")
    np.testing.assert_array_equal(values, expected)
    return values, profile


def main():
    global SOURCE_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Optional new directory; default: verify only")
    parser.add_argument("--source-dir", type=Path, required=True, help="External sources directory containing soil_0p1.tif.")
    args = parser.parse_args()
    SOURCE_ROOT = args.source_dir
    if args.output_dir is not None and args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    parameters = regional_parameters()
    values, profile = soil_layer()
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True)
        for kind, frame in parameters.items():
            frame.to_csv(args.output_dir / f"regional_parameters_{kind}.csv", index=False, encoding="utf-8-sig")
        with rasterio.open(args.output_dir / "static_soil_0p25.tif", "w", **profile) as dst:
            dst.write(values)
            for band, name in enumerate(["SOC", "soil_total_N", "soil_pH", "soil_valid_overlap_fraction"], 1):
                dst.set_band_description(band, name)
    print("PASS: 12-region parameters reconstructed; 4-band soil values exactly reproduced; no RF inference")


if __name__ == "__main__":
    main()
