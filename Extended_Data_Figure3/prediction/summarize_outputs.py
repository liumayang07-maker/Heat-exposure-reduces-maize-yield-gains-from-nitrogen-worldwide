"""Verify archived RF outputs and reconstruct summaries without RF inference."""
import argparse
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import rasterio

HERE = Path(__file__).resolve().parent
KEYS = ["year", "region_code", "N_rate", "irrigated"]


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def table(path):
    return pd.read_csv(path, float_precision="round_trip")


def raster_name(package, year, nitrogen, water):
    prefix = "" if package.startswith("historical") else package.split("_")[0] + "_"
    middle = "" if package.startswith("historical") else "_median_P50"
    return f"{prefix}{year}{middle}_N{nitrogen}_{water}.tif"


def annual_gains(annual):
    index = ["year", "region_code", "region_id", "region_name", "irrigated"]
    wide = annual.pivot(index=index, columns="N_rate", values="yield_Mg_ha")
    result = []
    for nitrogen in [120, 240]:
        frame = (wide[nitrogen] - wide[0]).rename("delta_yield_Mg_ha").reset_index()
        frame["N_rate"] = nitrogen
        result.append(frame)
    return pd.concat(result, ignore_index=True)


def compare(left, right, keys, columns, tolerance=1e-12):
    assert not left.duplicated(keys).any() and not right.duplicated(keys).any()
    left = left.set_index(keys).sort_index()
    right = right.set_index(keys).sort_index()
    pd.testing.assert_index_equal(left.index, right.index)
    np.testing.assert_allclose(left[columns], right[columns], atol=tolerance, rtol=0, equal_nan=True)


def summarize_package(package, bundle, output_dir=None, annual_root=None):
    name, years = package["name"], package["years"]
    folder = HERE / name
    grid = table(folder / "tables/prediction_grid.csv")
    annual = table(folder / "tables/regional_annual_yield.csv")
    expected = {(y, r, n, i) for y in years for r in grid.region_code.unique() for n in [0, 120, 240] for i in [0, 1]}
    assert set(map(tuple, annual[KEYS].to_numpy())) == expected
    assert len(annual) == len(expected) and annual.heat_scenario.eq("median_P50").all()
    indexed = annual.set_index(KEYS)
    ids = grid.pixel_id.to_numpy(int)
    regions = grid.region_code.to_numpy()
    weights = grid.fixed_area_weight.to_numpy()
    names = bundle["features"]
    matrices = {}
    for year in years:
        with np.load(folder / f"data/inputs_{year}_median_P50.npz", allow_pickle=False) as data:
            np.testing.assert_array_equal(data["pixel_id"], ids)
            np.testing.assert_array_equal(data["region_code"], regions)
            np.testing.assert_array_equal(data["area_weight"], weights)
            np.testing.assert_array_equal(data["features"], names)
            matrices[year] = data["X"].copy()
    max_weighted_error = 0.
    rasters_checked = 0
    for nitrogen in [0, 120, 240]:
        for irrigated, water in [(0, "rainfed"), (1, "irrigated")]:
            total = np.zeros(len(ids))
            for year in years:
                # 1. 核查完整掩膜与权重；TIFF 为 float32，原区域表由 float64 预测汇总。
                path = annual_root / name / "outputs" / raster_name(name, year, nitrogen, water)
                with rasterio.open(path) as src:
                    assert src.shape == (720, 1440) and src.crs.to_epsg() == 4326
                    assert src.count == 1 and src.dtypes == ("float32",) and src.nodata == -9999.
                    assert tuple(src.transform)[:6] == (.25, 0., -180., 0., -.25, 90.)
                    raw = src.read(1).ravel()
                    profile = src.profile.copy()
                    np.testing.assert_array_equal(np.flatnonzero(np.isfinite(raw) & (raw != src.nodata)), ids)
                    values = raw[ids].astype(float)
                assert np.isfinite(values).all() and (values >= 0).all()
                total += values
                matrix = matrices[year].copy()
                matrix[:, names.index("N_rate")] = nitrogen
                matrix[:, names.index("irrigated")] = irrigated
                outside = np.zeros(len(ids), dtype=bool)
                for column, feature in enumerate(names):
                    low, high = bundle["ranges"][feature]
                    outside |= (matrix[:, column] < low) | (matrix[:, column] > high)
                for region in np.unique(regions):
                    mask = regions == region
                    row = indexed.loc[(year, region, nitrogen, irrigated)]
                    assert row.pixels == mask.sum()
                    np.testing.assert_allclose(row.area_weight_sum, weights[mask].sum(), rtol=1e-13)
                    mean = float(np.average(values[mask], weights=weights[mask]))
                    max_weighted_error = max(max_weighted_error, abs(mean - row.yield_Mg_ha))
                    np.testing.assert_allclose(
                        [mean, values[mask].min(), values[mask].max()],
                        [row.yield_Mg_ha, row.min_yield_Mg_ha, row.max_yield_Mg_ha],
                        atol=1e-5,
                        rtol=0,
                    )
                    share = weights[mask & outside].sum() / weights[mask].sum()
                    np.testing.assert_allclose(share, row.outside_training_range_area_share, atol=1e-12, rtol=0)
                    within = mask & ~outside
                    within_yield = np.average(values[within], weights=weights[within]) if within.any() else np.nan
                    np.testing.assert_allclose(within_yield, row.within_range_yield_Mg_ha, atol=1e-5, rtol=0, equal_nan=True)
                rasters_checked += 1
            # 2. 按原代码对每年 float32 栅格等权平均，核对时期均值的全部有效像元。
            filename = f"mean_yield_{years[0]}_{years[-1]}_N{nitrogen}_{water}.tif"
            mean_values = (total / len(years)).astype("float32")
            with rasterio.open(folder / "summary" / filename) as src:
                assert src.shape == (720, 1440) and src.transform == rasterio.Affine(.25, 0, -180, 0, -.25, 90)
                raw = src.read(1).ravel()
                np.testing.assert_array_equal(np.flatnonzero(raw != src.nodata), ids)
                np.testing.assert_array_equal(raw[ids], mean_values)
            if output_dir is not None:
                target = output_dir / name / "summary" / filename
                target.parent.mkdir(parents=True, exist_ok=True)
                raster = np.full((720, 1440), -9999., dtype="float32")
                raster.ravel()[ids] = mean_values
                with rasterio.open(target, "w", **profile) as dst:
                    dst.write(raster, 1)
                    dst.set_band_description(1, f"{years[0]}_{years[-1]}_mean_yield_Mg_ha")
                    dst.update_tags(units="Mg/ha", temporal_weighting="equal years; fixed area support")
    keys = ["region_code", "N_rate", "irrigated"]
    means = annual.groupby(keys, as_index=False).agg(
        mean_yield_Mg_ha=("yield_Mg_ha", "mean"),
        sd_across_years_Mg_ha=("yield_Mg_ha", "std"),
        n_years=("year", "nunique"),
        mean_outside_training_range_area_share=("outside_training_range_area_share", "mean"),
    )
    compare(means, table(folder / "tables/regional_mean_yield.csv"), keys, list(means.columns[3:]))
    # 3. 同一年、同区域、同水分管理下先配对减 N0，再计算跨年均值与标准差。
    gains = annual_gains(annual)
    compare(gains, table(folder / "tables/regional_annual_N_gain.csv"), KEYS, ["delta_yield_Mg_ha"])
    heatmap = gains.groupby(["region_id", "irrigated", "N_rate"], as_index=False).agg(
        mean_Mg_ha=("delta_yield_Mg_ha", "mean"),
        interannual_sd_Mg_ha=("delta_yield_Mg_ha", "std"),
        years=("delta_yield_Mg_ha", "count"),
    )
    heatmap["dataset"] = name
    if output_dir is not None:
        target = output_dir / name / "tables"
        target.mkdir(parents=True, exist_ok=True)
        gains.to_csv(target / "regional_annual_N_gain.csv", index=False, encoding="utf-8-sig")
        means.to_csv(target / "regional_mean_yield.csv", index=False, encoding="utf-8-sig")
    print(f"PASS {name}: {rasters_checked} annual rasters, 6 mean rasters, {len(annual)} regional rows", flush=True)
    return heatmap, max_weighted_error, len(annual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Optional new directory for reproduced summaries")
    parser.add_argument("--annual-root", type=Path, help="External root containing <package>/outputs annual TIFFs; omit to check deposited files only.")
    args = parser.parse_args()
    if args.output_dir is not None and args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    inputs = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    outputs = json.loads((HERE / "output_manifest.json").read_text(encoding="utf-8"))
    assert digest(HERE / "manifest.json") == outputs["input_manifest_sha256"]
    model = HERE.parent / "RF_model/rf_model.joblib"
    assert digest(model) == outputs["model_sha256"]
    for record in [*inputs["files"], *outputs["files"]]:
        assert digest(HERE / record["path"]) == record["sha256"], record["path"]
    # 精简归档默认只核对现存结果；不要求未上传的年度栅格。
    if args.annual_root is None:
        if args.output_dir is not None:
            parser.error('--output-dir requires --annual-root for reconstructing summaries.')
        print('PASS: deposited NPZ inputs, period means and regional tables match their hashes; no inference.')
        return
    bundle = joblib.load(model)
    result = [summarize_package(p, bundle, args.output_dir, args.annual_root) for p in inputs["packages"]]
    heatmap = pd.concat([x[0] for x in result], ignore_index=True)
    assert len(heatmap) == 240
    compare(
        heatmap, table(HERE / "summary/regional_N_gain_source.csv"),
        ["dataset", "region_id", "irrigated", "N_rate"],
        ["mean_Mg_ha", "interannual_sd_Mg_ha", "years"],
    )
    if args.output_dir is not None:
        target = args.output_dir / "summary"
        target.mkdir(parents=True, exist_ok=True)
        heatmap.to_csv(target / "regional_N_gain_source.csv", index=False)
    print(json.dumps(dict(
        status="PASS_OUTPUT_AND_SUMMARY_INTEGRITY",
        annual_rasters=576,
        period_mean_rasters=30,
        annual_regional_rows=sum(x[2] for x in result),
        heatmap_rows=len(heatmap),
        maximum_raster_vs_float64_regional_mean_difference=max(x[1] for x in result),
        model_inference_executed=False,
    ), indent=2))


if __name__ == "__main__":
    main()
