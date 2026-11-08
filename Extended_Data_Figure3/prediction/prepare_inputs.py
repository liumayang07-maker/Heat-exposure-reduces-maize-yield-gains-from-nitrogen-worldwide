"""Rebuild the 96 archived RF matrices from local climate and shared inputs."""
import argparse
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import rasterio
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
SHARED = HERE / "shared"
CLIMATE = ["annual_DD10_30_C_d", "annual_EDH30_C_h", "annual_precip_mm"]
SOILS = ["SOC", "soil_total_N", "soil_pH"]
STAGES = [
    "EDH30_0_40_growing_season_fraction",
    "EDH30_40_70_growing_season_fraction",
    "EDH30_70_100_growing_season_fraction",
]
HEAT = ["growing_season_EDH30_fraction", *STAGES]
FRACTIONS = ["growing_season_GDD_fraction", "growing_season_precip_fraction", *HEAT]


def read_table(path):
    return pd.read_csv(path, float_precision="round_trip")


def overlap(ns, nt, source_origin, source_step, latitude=False):
    # 1. 沿用原版球面重叠权重，裁掉导出栅格超出地球边界的微小部分。
    origin, sign = (90., -1.) if latitude else (-180., 1.)
    bound = 90. if latitude else 180.
    source_edges = np.clip(source_origin + np.arange(ns + 1) * source_step, -bound, bound)
    target_edges = origin + sign * np.arange(nt + 1) * .25
    low = np.minimum(source_edges[:-1], source_edges[1:])
    high = np.maximum(source_edges[:-1], source_edges[1:])
    rows, columns, values = [], [], []
    for target in range(nt):
        target_low, target_high = sorted(target_edges[target:target + 2])
        for source in np.flatnonzero((high > target_low) & (low < target_high)):
            source_low, source_high = sorted(source_edges[source:source + 2])
            lo, hi = max(source_low, target_low), min(source_high, target_high)
            if hi <= lo:
                continue
            fraction = (
                (np.sin(np.deg2rad(hi)) - np.sin(np.deg2rad(lo)))
                / (np.sin(np.deg2rad(source_high)) - np.sin(np.deg2rad(source_low)))
            ) if latitude else (hi - lo) / (source_high - source_low)
            rows.append(target)
            columns.append(source)
            values.append(fraction)
    matrix = csr_matrix((values, (rows, columns)), shape=(nt, ns))
    np.testing.assert_allclose(np.asarray(matrix.sum(axis=0)).ravel(), 1., atol=1e-10)
    return matrix


def historical_climate():
    grid = read_table(SHARED / "historical_target_grid.csv")
    ids = grid.pixel_id.to_numpy(int)
    fixed = np.zeros((1800, 3600), dtype="float64")
    reference_transform = None
    # 2. 用 2011–2020 年均玉米面积作为固定权重，不把缺测气候填成零。
    for year in range(2011, 2021):
        path = SOURCE_ROOT / f"historical/MaizeCropYearClimate_{year}.tif"
        with rasterio.open(path) as src:
            assert src.shape == fixed.shape and src.crs.to_epsg() == 4326
            assert src.descriptions[3] == "maize_area"
            if reference_transform is None:
                reference_transform = src.transform
            assert src.transform == reference_transform
            area = src.read(4, masked=True).filled(0).astype(float)
        fixed += np.where(np.isfinite(area) & (area > 0), area, 0) / 10
    lat = overlap(1800, 720, reference_transform.f, reference_transform.e, True)
    lon = overlap(3600, 1440, reference_transform.c, reference_transform.a)

    def aggregate(values):
        return lat @ (lon @ values.T).T

    area = aggregate(fixed)
    np.testing.assert_allclose(area.sum(), fixed.sum(), rtol=1e-12)
    np.testing.assert_array_equal(area.ravel()[ids], grid.fixed_area_weight)
    for year in range(2011, 2023):
        with rasterio.open(SOURCE_ROOT / f"historical/MaizeCropYearClimate_{year}.tif") as src:
            assert src.transform == reference_transform and src.shape == fixed.shape
            raw = src.read([1, 2, 3], masked=True).filled(np.nan).astype(float)
        valid = np.isfinite(raw).all(axis=0) & (raw >= 0).all(axis=0)
        available = aggregate(np.where(valid, fixed, 0))
        complete = (area > 0) & np.isclose(available, area, rtol=1e-10, atol=0)
        values = np.full((3, 720, 1440), np.nan)
        for band in range(3):
            numerator = aggregate(np.where(valid, raw[band], 0) * fixed)
            np.divide(numerator, area, out=values[band], where=complete)
        yield year, values.reshape(3, -1)[:, ids].T


def prior_exclusions(scenario):
    paths = [SHARED / f"prior_exclusions_{scenario}.csv"]
    if scenario == "ssp245":
        paths.append(SHARED / "additional_exclusions_ssp245.csv")
    frames = []
    for path in paths:
        frame = read_table(path)
        if "pixel_id" not in frame:
            frame["pixel_id"] = frame.row.astype(int) * 1440 + frame.col.astype(int)
        frames.append(frame[["year", "pixel_id"]])
    return pd.concat(frames).drop_duplicates()


def future_climate(scenario, years, ids):
    excluded = prior_exclusions(scenario)
    for year in years:
        path = SOURCE_ROOT / f"{scenario}/MaizeFutureClimate_v2_ACCESS-CM2_{scenario}_{year}.tif"
        with rasterio.open(path) as src:
            assert src.shape == (720, 1440) and src.crs.to_epsg() == 4326
            assert tuple(src.transform)[:6] == (.25, 0., -180., 0., -.25, 90.)
            assert src.descriptions == ("DD1030", "EDH30", "Rain", "maize_area", "crop_year_type")
            values = src.read(masked=True).filled(np.nan).reshape(5, -1)[:, ids].T
        # 3. 只执行实际用于预测的剔除规则；z=3.5 诊断不触发剔除。
        removed = np.isin(ids, excluded.loc[excluded.year.eq(year), "pixel_id"])
        removed |= values[:, 2] > 8000
        values[removed] = np.nan
        yield year, values


def reconstruct(package, names, soil_raster, historical_values):
    historical = package["name"].startswith("historical")
    kind = "historical" if historical else "future"
    target = read_table(SHARED / f"{kind}_target_grid.csv")
    ids = target.pixel_id.to_numpy(int)
    parameters = read_table(SHARED / f"regional_parameters_{kind}.csv").set_index("region_code")
    dynamic = read_table(SHARED / "country_year_parameters.csv")
    soil = soil_raster.reshape(3, -1)[:, ids].T
    static = np.isfinite(soil).all(axis=1) & target.country_code.gt(0).to_numpy()
    years = package["years"]
    if historical:
        climate = historical_values
        logs = {
            year: target.country_code.map(dynamic.loc[dynamic.year.eq(year)].set_index("country_code").log_gdp_pc_ppp_prior3).to_numpy()
            for year in years
        }
        common = static.copy()
    else:
        logs2020 = target.country_code.map(dynamic.loc[dynamic.year.eq(2020)].set_index("country_code").log_gdp_pc_ppp_prior3).to_numpy()
        static &= np.isfinite(logs2020)
        logs = {year: logs2020 for year in years}
        common = static.copy()
        climate = None
        # 同一时期两个 SSP 使用相同的 21 年共同有效空间范围。
        for scenario in ["ssp245", "ssp585"]:
            current = {}
            for year, values in future_climate(scenario, years, ids):
                common &= np.isfinite(values).all(axis=1) & (values[:, :3] >= 0).all(axis=1)
                current[year] = values[:, :3]
                good = np.isfinite(values).all(axis=1)
                np.testing.assert_array_equal(values[good, 3], target.fixed_area_weight.to_numpy()[good])
                np.testing.assert_array_equal(values[good, 4], target.crop_year_type.to_numpy()[good])
            if package["name"].startswith(scenario):
                climate = current
    for year in years:
        common &= np.isfinite(climate[year]).all(axis=1) & np.isfinite(logs[year])
    grid = read_table(HERE / package["name"] / "tables/prediction_grid.csv")
    np.testing.assert_array_equal(ids[common], grid.pixel_id)
    assert len(grid) == package["rows"] and grid.region_code.nunique() == 12
    np.testing.assert_array_equal(target.region_code.to_numpy()[common], grid.region_code)
    np.testing.assert_array_equal(target.fixed_area_weight.to_numpy()[common], grid.fixed_area_weight)
    region = grid.region_code.to_numpy()
    for year in years:
        matrix = pd.DataFrame(0., index=np.arange(len(grid)), columns=names)
        matrix[CLIMATE] = climate[year][common]
        matrix[SOILS] = soil[common]
        matrix["Exp_year"] = year if historical else 2020.
        matrix["log_gdp_pc_ppp_prior3"] = logs[year][common]
        for feature in ["plant_density", *FRACTIONS]:
            matrix[feature] = pd.Series(region).map(parameters[feature]).to_numpy()
        zero = matrix.annual_EDH30_C_h.eq(0) | matrix.growing_season_EDH30_fraction.eq(0)
        matrix.loc[zero, HEAT] = 0.
        assert np.isfinite(matrix).all().all()
        np.testing.assert_allclose(matrix[STAGES].sum(axis=1), (~zero).astype(float), atol=1e-10)
        yield year, dict(
            X=matrix.to_numpy(),
            features=np.array(names),
            pixel_id=grid.pixel_id.to_numpy(int),
            region_code=region,
            area_weight=grid.fixed_area_weight.to_numpy(),
        )


def digest_source(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def run(output_dir=None, source_dir=None):
    global SOURCE_ROOT
    SOURCE_ROOT = Path(source_dir)
    # 4. 默认只验证；如需生成文件，必须指定全新输出目录，绝不覆盖归档。
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    if output_dir is not None:
        output_dir = output_dir.resolve()
        if output_dir.exists():
            raise FileExistsError(output_dir)
    for record in manifest.get("external_source_files", []):
        path = SOURCE_ROOT / Path(record["path"]).relative_to("sources")
        assert digest_source(path) == record["sha256"], path
    for record in manifest["files"]:
        with (HERE / record["path"]).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == record["sha256"], record["path"]
    for name, expected in manifest["dependencies"].items():
        with (HERE / name).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == expected, name
    names = joblib.load(HERE.parent / "RF_model/rf_model.joblib")["features"]
    with rasterio.open(SHARED / "static_soil_0p25.tif") as src:
        assert src.descriptions[:3] == tuple(SOILS)
        soil = src.read([1, 2, 3], masked=True).filled(np.nan)
    historical_values = dict(historical_climate())
    checked = 0
    maximum_difference = 0.
    for package in manifest["packages"]:
        for year, values in reconstruct(package, names, soil, historical_values):
            relative = Path(package["name"]) / "data" / f"inputs_{year}_median_P50.npz"
            with np.load(HERE / relative, allow_pickle=False) as archived:
                for key in values:
                    if key == "X":
                        difference = float(np.max(np.abs(values[key] - archived[key])))
                        maximum_difference = max(maximum_difference, difference)
                        np.testing.assert_allclose(values[key], archived[key], rtol=1e-13, atol=1e-10)
                    else:
                        np.testing.assert_array_equal(values[key], archived[key])
            if output_dir is not None:
                path = output_dir / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(path, **values)
            checked += 1
        print(f"PASS {package['name']}: {len(package['years'])} matrices, {package['rows']} cells", flush=True)
    result = dict(
        status="PASS_INPUT_RECONSTRUCTION",
        files_hashed=len(manifest["files"]),
        matrices_checked=checked,
        maximum_absolute_feature_difference=maximum_difference,
        model_inference_executed=False,
    )
    assert checked == 96
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Optional new directory; default: verify in memory")
    parser.add_argument("--source-dir", type=Path, required=True, help="External sources directory with historical/ and SSP climate rasters.")
    args = parser.parse_args()
    run(args.output_dir, args.source_dir)
