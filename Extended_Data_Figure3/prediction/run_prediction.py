"""Check or explicitly reproduce one archived RF prediction package."""
import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import rasterio

from summarize_outputs import HERE, annual_gains, digest, raster_name, table


def predict_batches(model, matrix, names, batch_size):
    values = np.empty(len(matrix), dtype="float64")
    for start in range(0, len(matrix), batch_size):
        stop = min(len(matrix), start + batch_size)
        values[start:stop] = model.predict(pd.DataFrame(matrix[start:stop], columns=names))
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Nonfinite or negative predicted yield")
    return values


def regional_summary(values, matrix, data, bundle, year, nitrogen, irrigated):
    outside = np.zeros(len(values), dtype=bool)
    for column, name in enumerate(bundle["features"]):
        low, high = bundle["ranges"][name]
        outside |= (matrix[:, column] < low) | (matrix[:, column] > high)
    management = f"N{nitrogen}_{'irrigated' if irrigated else 'rainfed'}"
    records = []
    for region in np.unique(data["region_code"]):
        mask = data["region_code"] == region
        weights = data["area_weight"][mask]
        within = mask & ~outside
        records.append(dict(
            year=year,
            heat_scenario="median_P50",
            management_id=management,
            scenario_id=f"median_P50_{management}",
            N_rate=nitrogen,
            irrigated=irrigated,
            region_code=int(region),
            pixels=int(mask.sum()),
            area_weight_sum=float(weights.sum()),
            yield_Mg_ha=float(np.average(values[mask], weights=weights)),
            min_yield_Mg_ha=float(values[mask].min()),
            max_yield_Mg_ha=float(values[mask].max()),
            outside_training_range_area_share=float(data["area_weight"][mask & outside].sum() / weights.sum()),
            within_range_yield_Mg_ha=float(np.average(values[within], weights=data["area_weight"][within])) if within.any() else None,
        ))
    return records


def save_raster(path, values, ids, profile, description):
    raster = np.full((profile["height"], profile["width"]), -9999., dtype="float32")
    raster.ravel()[ids] = values
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(raster, 1)
        dst.set_band_description(1, description)
        dst.update_tags(units="Mg/ha")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", required=True)
    parser.add_argument("--output-dir", type=Path, help="Required for --execute; must not already exist")
    parser.add_argument("--execute", action="store_true", help="Explicitly run RF inference; default checks inputs only")
    parser.add_argument("--accept-provisional-climate", action="store_true")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=20000)
    args = parser.parse_args()
    if args.workers < 1 or args.batch_size < 1:
        parser.error("workers and batch-size must be positive")
    inputs = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    packages = {package["name"]: package for package in inputs["packages"]}
    if args.package not in packages:
        parser.error("package must be one of: " + ", ".join(packages))
    package = packages[args.package]
    historical = args.package.startswith("historical")
    if args.execute:
        if args.output_dir is None:
            parser.error("--execute requires --output-dir")
        if args.output_dir.exists():
            raise FileExistsError(args.output_dir)
        if not historical and not args.accept_provisional_climate:
            parser.error("Future inference requires --accept-provisional-climate")
    # 1. 默认只核对同一归档模型、年度输入与管理组合，不进行任何预测。
    model_path = HERE.parent / "RF_model/rf_model.joblib"
    assert digest(model_path) == inputs["dependencies"]["../RF_model/rf_model.joblib"]
    package_files = [record for record in inputs["files"] if record["path"].startswith(args.package + "/")]
    for record in package_files:
        assert digest(HERE / record["path"]) == record["sha256"], record["path"]
    for name in ["shared/management_scenarios.csv", "shared/regional_parameters_historical.csv"]:
        record = next(record for record in inputs["files"] if record["path"] == name)
        assert digest(HERE / name) == record["sha256"], name
    management = table(HERE / "shared/management_scenarios.csv")
    assert set(zip(management.N_rate, management.irrigated)) == {(n, i) for n in [0, 120, 240] for i in [0, 1]}
    if not args.execute:
        print(f"PASS {args.package}: {len(package['years'])} inputs, {len(package['years']) * 6} management runs; no inference")
        return
    bundle = joblib.load(model_path)
    bundle["model"].set_params(n_jobs=args.workers)
    names = bundle["features"]
    folder = HERE / args.package
    destination = args.output_dir / args.package
    for name in ["outputs", "tables", "summary"]:
        (destination / name).mkdir(parents=True, exist_ok=True)
    with rasterio.open(folder / "data/prediction_region_mask.tif") as src:
        profile = src.profile.copy()
    profile.update(driver="GTiff", count=1, dtype="float32", nodata=-9999., compress="deflate", predictor=3)
    sums = {(n, i): np.zeros(package["rows"]) for n in [0, 120, 240] for i in [0, 1]}
    records, artifacts = [], []
    # 2. 沿用原版分批预测与面积加权汇总，原始区域表保留 float64 精度。
    for year in package["years"]:
        with np.load(folder / f"data/inputs_{year}_median_P50.npz", allow_pickle=False) as source:
            data = {key: source[key].copy() for key in source.files}
        np.testing.assert_array_equal(data["features"], names)
        for row in management.itertuples(index=False):
            nitrogen, irrigated = int(row.N_rate), int(row.irrigated)
            matrix = data["X"].copy()
            matrix[:, names.index("N_rate")] = nitrogen
            matrix[:, names.index("irrigated")] = irrigated
            values = predict_batches(bundle["model"], matrix, names, args.batch_size)
            records.extend(regional_summary(values, matrix, data, bundle, year, nitrogen, irrigated))
            path = destination / "outputs" / raster_name(args.package, year, nitrogen, row.water_regime)
            save_raster(path, values, data["pixel_id"], profile, "predicted_yield_Mg_ha")
            sums[(nitrogen, irrigated)] += values.astype("float32").astype(float)
            artifacts.append(dict(path=path.relative_to(args.output_dir).as_posix(), sha256=digest(path)))
        print(f"Predicted {args.package} {year}: six management runs", flush=True)
    lookup = table(HERE / "shared/regional_parameters_historical.csv")[["region_code", "region_id", "region_name"]]
    annual = pd.DataFrame(records).merge(lookup, on="region_code", validate="many_to_one")
    assert len(annual) == len(package["years"]) * 72
    annual.to_csv(destination / "tables/regional_annual_yield.csv", index=False, encoding="utf-8-sig")
    means = annual.groupby(["region_code", "region_id", "region_name", "N_rate", "irrigated"], as_index=False).agg(
        mean_yield_Mg_ha=("yield_Mg_ha", "mean"),
        sd_across_years_Mg_ha=("yield_Mg_ha", "std"),
        n_years=("year", "nunique"),
        mean_outside_training_range_area_share=("outside_training_range_area_share", "mean"),
    )
    means.to_csv(destination / "tables/regional_mean_yield.csv", index=False, encoding="utf-8-sig")
    gains = annual_gains(annual)
    gains.to_csv(destination / "tables/regional_annual_N_gain.csv", index=False, encoding="utf-8-sig")
    heatmap = gains.groupby(["region_id", "irrigated", "N_rate"], as_index=False).agg(
        mean_Mg_ha=("delta_yield_Mg_ha", "mean"),
        interannual_sd_Mg_ha=("delta_yield_Mg_ha", "std"),
        years=("delta_yield_Mg_ha", "count"),
    )
    heatmap["dataset"] = args.package
    heatmap.to_csv(destination / "tables/regional_N_gain_source.csv", index=False)
    for (nitrogen, irrigated), total in sums.items():
        water = "irrigated" if irrigated else "rainfed"
        path = destination / "summary" / f"mean_yield_{package['years'][0]}_{package['years'][-1]}_N{nitrogen}_{water}.tif"
        save_raster(path, total / len(package["years"]), data["pixel_id"], profile, "period_mean_yield_Mg_ha")
    # 3. 仅写一份总清单，不生成每个栅格的重复 JSON。
    result = dict(
        status="PASS_CONDITIONAL_RF_PREDICTION" if historical else "PASS_TRIAL_PREDICTION",
        model_sha256=digest(model_path),
        input_manifest_sha256=digest(HERE / "manifest.json"),
        annual_rasters=len(artifacts),
        annual_regional_rows=len(annual),
        artifacts=artifacts,
        interpretation="Conditional RF predictions; nitrogen contrasts are not causal effects",
        climate_limitations="Future daily-sine heat is not calibrated to hourly ERA5; one climate model; input screening does not validate extrapolation",
    )
    (destination / "prediction_manifest.json").write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
