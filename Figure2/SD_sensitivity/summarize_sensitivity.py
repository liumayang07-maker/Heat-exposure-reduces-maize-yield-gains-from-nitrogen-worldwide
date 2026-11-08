"""Summarize existing SD sensitivity results into core CSV files only.

Input root contains run_status.csv and results/<model>/<scenario>/*.csv.
Use --input-root for temporary runner results and --output-dir for retained CSVs.
No models are fitted and no input files are deleted.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SCENARIOS = ("Baseline", "CV05", "CV10", "CV20")
MODELS = {
    "Figure2": ("Figure2a", "Figure2b", "Figure2c", "Figure2d"),
    "Figure3": ("Figure3b",),
}
# 1. 按已核对的原始 CSV 和绘图坐标定义映射，不重算预测值。
CURVES = {
    "Figure2b": ("N_group_EDH30_percentile_curves_REML.csv", "percentile",
                 ["N_group"], "ci_lb_pct", "ci_ub_pct",
                 "log1p(EDH30) percentile within N-rate group"),
    "Figure2c": ("adjusted_N_response_curves_REML.csv", "N_rate",
                 ["Env_Type"], "effect_size_ci_lb", "effect_size_ci_ub",
                 "N rate (kg N ha-1)"),
    "Figure2d": ("fixed_N_percentile_curves_REML.csv", "percentile",
                 ["Env_Type", "N_control_kg_ha"], "effect_size_ci_lb",
                 "effect_size_ci_ub", "log1p(EDH30) percentile within heat environment"),
    "Figure3b": ("regional_N_EDH30_predictions_REML.csv", "log1p_EDH30",
                 ["region", "N_scenario"], "ci_lb_pct", "ci_ub_pct",
                 "Within-region ln(1+EDH30)"),
}
SUCCESS = {"success", "ok", "completed", "complete", "passed", "done"}


def read_csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def tagged(frame, model, scenario, source):
    frame = frame.copy()
    if "model" in frame:
        frame = frame.rename(columns={"model": "source_model"})
    frame.insert(0, "model", model)
    frame.insert(1, "scenario", scenario)
    frame.insert(2, "source_file", source)
    return frame


def combine(frames):
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def direction(values):
    values = np.asarray(values)
    positive = np.any(values > 1e-10)
    negative = np.any(values < -1e-10)
    if positive and negative:
        return "mixed"
    return "positive" if positive else "negative" if negative else "flat"


def summarize_curves(predictions, model):
    _, x, groups, lower, upper, _ = CURVES[model]
    summaries = []
    for keys, curve in predictions.groupby(["scenario"] + groups, sort=False):
        scenario, *values = keys
        curve = curve.sort_values(x)
        metadata = dict(zip(groups, values))
        row = {"model": model, "scenario": scenario, **metadata, "x_variable": x}
        row["n"] = curve["n"].iloc[0]
        y = curve["effect_size_pct"].to_numpy()
        row.update(
            points=len(curve),
            x_min=curve[x].iloc[0],
            x_max=curve[x].iloc[-1],
            effect_at_start=y[0],
            effect_at_end=y[-1],
            ci_lower_at_start=curve[lower].iloc[0],
            ci_upper_at_start=curve[upper].iloc[0],
            ci_lower_at_end=curve[lower].iloc[-1],
            ci_upper_at_end=curve[upper].iloc[-1],
            effect_min=y.min(),
            effect_max=y.max(),
            effect_range_pp=y.max() - y.min(),
            endpoint_change_pp=y[-1] - y[0],
            effect_sign=direction(y),
            curve_direction=direction(np.diff(y)),
        )
        baseline = predictions.loc[predictions["scenario"].eq("Baseline")]
        for name, value in metadata.items():
            baseline = baseline.loc[baseline[name].eq(value)]
        # 相同百分位并不保证相同暴露；同时匹配所有实际预测坐标。
        grid = list(dict.fromkeys(groups + [x] + [
            name for name in ("log1p_EDH30", "EDH30_total", "N_rate", "N_control_kg_ha")
            if name in curve
        ]))
        measures = ["effect_size_pct", lower, upper]
        matched = curve[grid + measures].merge(
            baseline[grid + measures],
            on=grid,
            how="outer",
            suffixes=("", "_baseline"),
            indicator=True,
            validate="one_to_one",
        )
        same_grid = not baseline.empty and matched["_merge"].eq("both").all()
        row["baseline_grid_status"] = "same grid" if same_grid else "missing baseline" if baseline.empty else "grid mismatch"
        row["matched_points"] = int(matched["_merge"].eq("both").sum())
        if same_grid:
            matched = matched.drop(columns="_merge").sort_values(x)
            for column, name in zip(measures, ("effect_shift_pp", "ci_lower_shift_pp", "ci_upper_shift_pp")):
                matched[name] = matched[column] - matched[column + "_baseline"]
            matched["ci_width_change_pp"] = (
                matched[upper] - matched[lower]
                - matched[upper + "_baseline"] + matched[lower + "_baseline"]
            )
            for name in ("effect_shift_pp", "ci_lower_shift_pp", "ci_upper_shift_pp", "ci_width_change_pp"):
                row[name + "_mean"] = matched[name].mean()
                row[name + "_max_abs"] = matched[name].abs().max()
                peak = matched.loc[matched[name].abs().idxmax()]
                row[name + "_at_max_abs"] = peak[name]
                row[name + "_max_abs_x"] = peak[x]
            row["endpoint_change_shift_pp"] = matched["effect_shift_pp"].iloc[-1] - matched["effect_shift_pp"].iloc[0]

        summaries.append(row)
    return pd.DataFrame(summaries)


def check_predictions(frame, model):
    _, x, groups, lower, upper, _ = CURVES[model]
    required = groups + [x, "effect_size_pct", lower, upper]
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing prediction columns: {sorted(missing)}")
    if frame.empty or frame[required].isna().any().any():
        raise ValueError("Empty predictions or missing curve values")
    if not np.isfinite(frame[[x, "effect_size_pct", lower, upper]].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite prediction values")
    if frame.duplicated(groups + [x]).any():
        raise ValueError("Duplicate curve coordinates")
    if (frame[lower] > frame[upper]).any():
        raise ValueError("Reversed confidence interval bounds")


def build_report(root, figure, output_dir=None):
    root = Path(root)
    output_dir = Path(output_dir) if output_dir is not None else root
    # 2. 只消费 runner 的完成状态；运行中/失败目录不作为成功结果。
    status_path = root / "run_status.csv"
    status = read_csv(status_path)
    if not {"model", "scenario", "status"}.issubset(status):
        raise ValueError("run_status.csv requires model, scenario, status")
    if status.duplicated(["model", "scenario"]).any():
        raise ValueError("Duplicate model/scenario in run_status.csv")
    expected = pd.MultiIndex.from_product([MODELS[figure], SCENARIOS], names=["model", "scenario"]).to_frame(index=False)
    status = expected.merge(status, on=["model", "scenario"], how="left")
    status["status"] = status["status"].fillna("not recorded")
    buckets = {name: [] for name in ("Overall", "Key_coefficients", "Regional_heat", "Tests", "Variance_components", "Model_parameters")}
    issues = []
    predictions = {model: [] for model in MODELS[figure] if model in CURVES}
    for index, record in status.iterrows():
        model, scenario = record["model"], record["scenario"]
        folder = root / "results" / model / scenario
        if str(record["status"]).lower() not in SUCCESS:
            status.loc[index, "report_status"] = "not consumed"
            continue
        required = (["overall_REML_results.csv"] if model == "Figure2a" else [
            "coefficients_REML.csv", "variance_components_REML.csv", "model_parameters.csv", CURVES[model][0]
        ])
        if model == "Figure3b":
            required += ["regional_heat_coefficients.csv", "regional_heat_N_joint_tests_REML.csv"]
        required += {
            "Figure2b": ["interaction_joint_Wald.csv"],
            "Figure2c": ["interaction_Wald_tests.csv"],
            "Figure2d": ["three_way_interaction_Wald_tests.csv", "fixed_N_curve_Wald_tests_REML.csv"],
            "Figure3b": ["omnibus_Wald_tests_REML.csv"],
        }.get(model, [])
        missing = [name for name in required if not (folder / name).is_file()]
        status.loc[index, "report_status"] = "incomplete files" if missing else "consumed"
        for name in missing:
            issues.append({"model": model, "scenario": scenario, "file": name, "issue": "Missing required output"})
        # n 为该模型全部效应量行数，不冒充某区域/曲线的独立样本数。
        sample_size = np.nan
        if model != "Figure2a":
            try:
                parameters = read_csv(folder / "model_parameters.csv")
                if len(parameters) != 1 or "effects" not in parameters:
                    raise ValueError("Expected one model_parameters row with effects")
                sample_size = pd.to_numeric(parameters["effects"], errors="raise").iloc[0]
                if not np.isfinite(sample_size) or sample_size <= 0 or sample_size != int(sample_size):
                    raise ValueError("effects must be a positive integer")
            except (ValueError, OSError, pd.errors.ParserError) as error:
                status.loc[index, "report_status"] = "invalid output"
                issues.append({"model": model, "scenario": scenario, "file": "model_parameters.csv", "issue": str(error)})
        for path in sorted(folder.glob("*.csv")):
            name = path.name
            sheet = {
                "overall_REML_results.csv": "Overall",
                "coefficients_REML.csv": "Key_coefficients",
                "regional_heat_coefficients.csv": "Regional_heat",
                "variance_components_REML.csv": "Variance_components",
                "model_parameters.csv": "Model_parameters",
            }.get(name)
            if "wald" in name.lower() or name == "regional_heat_N_joint_tests_REML.csv":
                sheet = "Tests"
            is_prediction = model in CURVES and name == CURVES[model][0]
            if sheet is None and not is_prediction:
                continue
            try:
                frame = read_csv(path)
                if frame.empty:
                    raise ValueError("Empty output table")
                if model != "Figure2a":
                    frame["n"] = sample_size
                if is_prediction:
                    check_predictions(frame, model)
                    predictions[model].append(tagged(frame, model, scenario, name))
                else:
                    buckets[sheet].append(tagged(frame, model, scenario, name))
            except (ValueError, KeyError, OSError, pd.errors.ParserError) as error:
                status.loc[index, "report_status"] = "invalid output"
                issues.append({"model": model, "scenario": scenario, "file": name, "issue": str(error)})

    # 3. 仅保留曲线摘要；CI 边界变化不等于情景间差异的置信区间。
    curve_summaries = [
        summarize_curves(combine(frames), model)
        for model, frames in predictions.items()
        if frames
    ]
    tables = {name: combine(frames) for name, frames in buckets.items()}
    # 运行状态仅保留核心字段，不复制耗时或其他运行元数据。
    core_status = status.reindex(columns=[
        "model",
        "scenario",
        "status",
        "warnings",
        "error",
        "report_status",
    ])
    tables.update({
        "Model_status": core_status,
        "Curve_summaries": combine(curve_summaries),
        "Reporting_issues": pd.DataFrame(
            issues,
            columns=["model", "scenario", "file", "issue"],
        ),
    })
    filenames = {
        "Overall": "overall.csv",
        "Key_coefficients": "coefficients.csv",
        "Regional_heat": "regional_heat_coefficients.csv",
        "Tests": "joint_tests.csv",
        "Variance_components": "variance_components.csv",
        "Model_parameters": "model_parameters.csv",
        "Model_status": "model_status.csv",
        "Curve_summaries": "curve_summaries.csv",
        "Reporting_issues": "reporting_issues.csv",
    }
    # 4. 仅输出有内容的核心表；检查失败时不覆盖已经核验的结果。
    if issues or not status["report_status"].eq("consumed").all():
        raise ValueError(f"Incomplete results: {issues}")
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in tables.items():
        if frame.empty:
            continue
        frame = frame.drop(columns=["caution", "source_file", "source_model"], errors="ignore")
        if name == "Key_coefficients":
            terms = frame["term"].astype(str)
            focus = terms.str.contains(r"heat|EDH|N_c|N_group|N_rate|Env_Type|region", case=False)
            frame = frame.loc[(focus | terms.eq("intrcpt")) & ~terms.str.contains("country", case=False)]
        elif name == "Overall":
            fields = ["model", "scenario", "k", "estimate_lnRR", "se", "ci_lb_lnRR", "ci_ub_lnRR",
                      "pval", "response_change_percent", "ci_lb_response_change_percent",
                      "ci_ub_response_change_percent", "sigma2_reference_ID", "sigma2_ref_site_year_ID",
                      "sigma2_effect_id", "QE", "QE_df", "QE_pval", "unique_references"]
            frame = frame[fields].rename(columns={"k": "n"})
        elif name == "Model_parameters":
            fields = ["model", "scenario", "n", "references", "environments", "regions",
                      "design_rank", "method", "optimizer", "logLik_REML", "AIC_REML", "BIC_REML",
                      "QE", "QE_p", "R_version", "metafor_version"]
            frame = frame[[field for field in fields if field in frame]]
            frame["QE_df"] = frame["n"] - frame["design_rank"]
        elif name == "Curve_summaries":
            fields = ["model", "scenario", "N_group", "Env_Type", "N_control_kg_ha", "region",
                      "N_scenario", "n", "x_variable", "x_min", "x_max", "effect_at_start",
                      "ci_lower_at_start", "ci_upper_at_start", "effect_at_end", "ci_lower_at_end",
                      "ci_upper_at_end", "endpoint_change_pp", "effect_sign", "curve_direction",
                      "effect_shift_pp_max_abs", "effect_shift_pp_at_max_abs", "effect_shift_pp_max_abs_x",
                      "ci_lower_shift_pp_max_abs", "ci_upper_shift_pp_max_abs", "ci_width_change_pp_max_abs"]
            frame = frame[[field for field in fields if field in frame]]
        if name == "Model_status":
            frame = frame[["model", "scenario", "status", "warnings", "error"]]
        frame.to_csv(output_dir / filenames[name], index=False)
    print(f"Core CSVs written: {output_dir}")
    print(f"Consumed: {(status['report_status'] == 'consumed').sum()}/{len(status)}; reporting issues: {len(issues)}")
    return tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", choices=tuple(MODELS), default="Figure2")
    parser.add_argument(
        "--input-root",
        type=Path,
        help="Directory containing run_status.csv and results/ (default: figure SD_sensitivity)",
    )
    parser.add_argument(
        "--output-dir",
        "--output-root",
        dest="output_dir",
        type=Path,
        help="Core CSV destination (default: input root)",
    )
    args = parser.parse_args()
    manuscripts = Path(__file__).resolve().parents[2]
    root = args.input_root or manuscripts / args.figure / "SD_sensitivity"
    tables = build_report(root, args.figure, args.output_dir)
    # 不完整迁移返回非零状态，调用者不得据此清理原始结果。
    if not tables["Model_status"]["report_status"].eq("consumed").all():
        raise SystemExit(1)


if __name__ == "__main__":
    main()
