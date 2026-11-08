from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold


HERE = Path(__file__).resolve().parent


def metrics(y, prediction):
    return dict(
        n=len(y),
        R2=float(r2_score(y, prediction)),
        RMSE=float(mean_squared_error(y, prediction) ** 0.5),
        MAE=float(mean_absolute_error(y, prediction)),
        bias=float(np.mean(prediction - y)),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce nested grouped RF training in a new output directory")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    # 1. 只读完整精确训练快照，不覆盖归档模型或既有结果。
    if args.output_dir.exists():
        raise FileExistsError(f"Use a new output directory: {args.output_dir}")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    config = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    if sklearn.__version__ != config["sklearn"]:
        raise RuntimeError(f"Reproduction requires scikit-learn {config['sklearn']}")
    data = pd.read_csv(HERE / "tables/training_arms_v7.csv", float_precision="round_trip")
    X = data[config["features"]]
    y = data.yield_Mg_ha.to_numpy(dtype=float)
    groups = data.validation_component
    folds = data.fold.to_numpy()
    assert len(data) == 6040 and np.isfinite(X.to_numpy()).all() and np.isfinite(y).all()
    output = args.output_dir.resolve()
    tables = output / "tables"
    tables.mkdir(parents=True)

    def save(frame, name):
        frame.to_csv(tables / name, index=False, encoding="utf-8-sig")

    def tune(indices, label):
        # 2. 使用归档的三折内层划分，48 组参数、RMSE 选择准则不变。
        membership = pd.read_csv(HERE / f"tables/inner_folds_{label}.csv")
        assert set(membership.row_index) == set(indices) and membership.row_index.is_unique
        labels = membership.set_index("row_index").inner_fold.reindex(indices).to_numpy()
        splits = [(np.flatnonzero(labels != fold), np.flatnonzero(labels == fold))
                  for fold in sorted(np.unique(labels))]
        for train, test in splits:
            assert not set(groups.iloc[indices[train]]) & set(groups.iloc[indices[test]])
        search = GridSearchCV(
            RandomForestRegressor(random_state=config["seed"], n_jobs=1),
            param_grid=config["grid"],
            cv=splits,
            scoring="neg_root_mean_squared_error",
            n_jobs=args.jobs,
            pre_dispatch=args.jobs,
            error_score="raise",
            return_train_score=False,
        )
        search.fit(X.iloc[indices], y[indices])
        save(pd.DataFrame(search.cv_results_), f"tuning_{label}.csv")
        save(membership, f"inner_folds_{label}.csv")
        return search.best_estimator_, search.best_params_

    # 3. 五折外层仅在训练折内调参，保存每条样本的留出预测。
    oof = np.full(len(data), np.nan)
    fold_metrics, parameters = [], []
    for fold in sorted(data.fold.unique()):
        train = np.flatnonzero(folds != fold)
        test = np.flatnonzero(folds == fold)
        assert not set(groups.iloc[train]) & set(groups.iloc[test])
        assert not set(data.reference_ID.iloc[train]) & set(data.reference_ID.iloc[test])
        model, best = tune(train, f"outer_{fold}")
        oof[test] = model.predict(X.iloc[test])
        fold_metrics.append(dict(fold=int(fold), **metrics(y[test], oof[test])))
        parameters.append(dict(fold=int(fold), **best))
    assert np.isfinite(oof).all()
    save(pd.DataFrame(fold_metrics), "fold_metrics.csv")
    save(pd.DataFrame(parameters), "fold_best_parameters.csv")
    save(data[["candidate_arm_id", "reference_ID", "country", "Exp_year", "yield_Mg_ha", "fold"]].assign(
        prediction_Mg_ha=oof, residual_Mg_ha=y - oof), "out_of_fold_predictions.csv")
    save(pd.DataFrame([metrics(y, oof)]), "overall_metrics.csv")

    # 4. 全数据内层调参后拟合最终模型；等权、不裁剪产量。
    model, _ = tune(np.arange(len(data)), "full_data")
    model.set_params(n_jobs=args.jobs)
    bundle = joblib.load(HERE / "rf_model.joblib")
    bundle["model"] = model
    bundle["training_csv_sha256"] = hashlib.sha256((HERE / "tables/training_arms_v7.csv").read_bytes()).hexdigest()
    joblib.dump(bundle, output / "rf_model.joblib", compress=3)

    # 5. 时间和国家留出沿用固定最终参数，仅作为迁移诊断。
    temporal = []
    for cutoff in [2010, 2015]:
        years = data.groupby("reference_ID").Exp_year.agg(["min", "max"])
        train_refs = set(years.index[years["max"] <= cutoff])
        test_refs = set(years.index[years["min"] > cutoff])
        train = data.reference_ID.isin(train_refs).to_numpy()
        test = data.reference_ID.isin(test_refs).to_numpy()
        fitted = RandomForestRegressor(**model.get_params()).fit(X.loc[train], y[train])
        temporal.append(dict(cutoff=cutoff, train_n=int(train.sum()), test_n=int(test.sum()),
                             train_studies=len(train_refs), test_studies=len(test_refs),
                             **metrics(y[test], fitted.predict(X.loc[test]))))
    save(pd.DataFrame(temporal), "temporal_forward_validation.csv")
    country_prediction = np.full(len(data), np.nan)
    country_fold = np.full(len(data), -1, dtype=int)
    country_metrics = []
    for fold, (train, test) in enumerate(GroupKFold(n_splits=5).split(X, y, data.country), 1):
        fitted = RandomForestRegressor(**model.get_params()).fit(X.iloc[train], y[train])
        country_prediction[test] = fitted.predict(X.iloc[test])
        country_fold[test] = fold
        country_metrics.append(dict(fold=fold, **metrics(y[test], country_prediction[test])))
    country_metrics.append(dict(fold="overall", **metrics(y, country_prediction)))
    save(pd.DataFrame(country_metrics), "country_holdout_metrics.csv")
    save(data[["candidate_arm_id", "reference_ID", "country", "yield_Mg_ha"]].assign(
        country_fold=country_fold, prediction_Mg_ha=country_prediction), "country_holdout_predictions.csv")
    save(data, "training_arms_v7.csv")
    save(data[["candidate_arm_id", "reference_ID", "validation_component", "country", "Exp_year", "fold"]], "outer_folds.csv")
    config.update(parameters=model.get_params(), overall=metrics(y, oof),
                  model_sha256=hashlib.sha256((output / "rf_model.joblib").read_bytes()).hexdigest())
    (output / "manifest.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    print(f"Completed training in {output}")


if __name__ == "__main__":
    main()
