from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict maize yield with the 17-feature RF v7 model")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError("Output must differ from input")

    # 1. 读取完整特征并检查有效范围；零热暴露不能作为缺失值。
    bundle = joblib.load(Path(__file__).resolve().parent / "rf_model.joblib")
    data = pd.read_csv(args.input)
    X = data[bundle["features"]].apply(pd.to_numeric, errors="raise").copy()
    heat_timing = [
        "EDH30_0_40_growing_season_fraction",
        "EDH30_40_70_growing_season_fraction",
        "EDH30_70_100_growing_season_fraction",
    ]
    if not np.isfinite(X.to_numpy()).all():
        raise ValueError("Missing or nonfinite RF inputs")
    fractions = [name for name in X if "fraction" in name]
    if ((X[fractions] < 0) | (X[fractions] > 1 + 1e-6)).any().any():
        raise ValueError("Fraction outside [0, 1]")
    if not X.irrigated.isin([0, 1]).all():
        raise ValueError("Irrigation must be encoded as 0 or 1")
    if X.loc[X.growing_season_EDH30_fraction.eq(0), heat_timing].ne(0).any().any():
        raise ValueError("Nonzero heat-timing fraction with zero growing-season EDH30 fraction")
    positive = X.growing_season_EDH30_fraction.gt(0)
    if not np.allclose(X.loc[positive, heat_timing].sum(axis=1), 1.0, rtol=0, atol=1e-6):
        raise ValueError("Positive-exposure heat-timing fractions must sum to one")

    # 2. 使用归档模型预测，同时标识训练范围以外的输入。
    data["predicted_yield_Mg_ha"] = bundle["model"].predict(X)
    outside = [
        (X[feature] < lower) | (X[feature] > upper)
        for feature, (lower, upper) in bundle["ranges"].items()
    ]
    data["outside_training_range"] = np.logical_or.reduce(outside)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"Predicted {len(data)} rows; outside ranges: {int(data.outside_training_range.sum())}")


if __name__ == "__main__":
    main()
