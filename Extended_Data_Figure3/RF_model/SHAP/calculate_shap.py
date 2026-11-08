from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import shap


HERE = Path(__file__).resolve().parent
MODEL = HERE.parent


def explain_batch(model, frame):
    explainer = shap.TreeExplainer(
        model,
        feature_perturbation="tree_path_dependent",
        model_output="raw",
    )
    return np.asarray(explainer.shap_values(frame, check_additivity=True))


def main():
    parser = argparse.ArgumentParser(description="Recalculate complete RF SHAP values without refitting the RF")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    # 1. 共用正式模型及完整训练集，保留已有解释结果。
    if args.output_dir.exists():
        raise FileExistsError(f"Use a new output directory: {args.output_dir}")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    assert hashlib.sha256((MODEL / "rf_model.joblib").read_bytes()).hexdigest() == manifest["model_sha256"]
    assert hashlib.sha256((MODEL / "tables/training_arms_v7.csv").read_bytes()).hexdigest() == manifest["training_table_sha256"]
    if shap.__version__ != manifest["shap_version"]:
        raise RuntimeError(f"Reproduction requires shap {manifest['shap_version']}")
    bundle = joblib.load(MODEL / "rf_model.joblib")
    data = pd.read_csv(MODEL / "tables/training_arms_v7.csv")
    X = data[bundle["features"]].apply(pd.to_numeric, errors="raise")
    assert X.shape == (6040, 17) and np.isfinite(X.to_numpy()).all()

    # 2. 沿用原 TreeExplainer 口径、500 行批次及完整样本顺序，不抽样。
    model = bundle["model"]
    values = np.concatenate(joblib.Parallel(n_jobs=args.jobs)(
        joblib.delayed(explain_batch)(model, X.iloc[start:start + 500])
        for start in range(0, len(X), 500)
    ))
    explainer = shap.TreeExplainer(model, feature_perturbation="tree_path_dependent", model_output="raw")
    base = float(np.asarray(explainer.expected_value).reshape(-1)[0])
    prediction = model.predict(X)
    error = float(np.max(np.abs(base + values.sum(axis=1) - prediction)))
    assert error < 1e-6

    # 3. 保留完整逐行贡献、预测分解和特征重要性，不生成弃用单图。
    tables = args.output_dir / "tables"
    tables.mkdir(parents=True)
    ids = ["candidate_arm_id", "reference_ID", "env_ID", "arm", "country", "Exp_year"]
    # Exp_year 元数据与同名 SHAP 列沿用原格式，绘图代码按位置读取贡献列。
    pd.concat([data[ids].reset_index(drop=True), pd.DataFrame(values, columns=X.columns)], axis=1).to_csv(
        tables / "shap_values.csv", index=False, encoding="utf-8-sig")
    data[ids].assign(base_value=base, prediction_Mg_ha=prediction, shap_sum=values.sum(axis=1),
                     reconstructed_prediction_Mg_ha=base + values.sum(axis=1)).to_csv(
        tables / "prediction_decomposition.csv", index=False, encoding="utf-8-sig")
    absolute = pd.DataFrame(np.abs(values), columns=X.columns)
    study_equal = absolute.groupby(data.reference_ID).mean().mean()
    metadata = manifest["feature_metadata"]
    importance = pd.DataFrame(dict(
        feature=X.columns,
        feature_label=[metadata[name]["feature_label"] for name in X],
        feature_block=[metadata[name]["feature_block"] for name in X],
        mean_abs_SHAP=absolute.mean().to_numpy(),
        study_equal_mean_abs_SHAP=study_equal.reindex(X.columns).to_numpy(),
    )).sort_values("mean_abs_SHAP", ascending=False)
    importance["rank"] = np.arange(1, len(importance) + 1)
    importance["contribution_share"] = importance.mean_abs_SHAP / importance.mean_abs_SHAP.sum()
    importance.to_csv(tables / "importance.csv", index=False, encoding="utf-8-sig")
    manifest.update(base_value=base, max_additivity_error=error,
                    input_files={f"tables/{path.name}": hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in tables.glob("*.csv")})
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Completed {len(X)} rows; additivity error {error:.3g}; RF not refitted")


if __name__ == "__main__":
    main()
