from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams.update({
    "font.weight": "normal",
    "font.style": "normal",
    "axes.titleweight": "normal",
    "axes.labelweight": "normal",
    "mathtext.fontset": "custom",
    "mathtext.default": "regular",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial",
    "mathtext.bf": "Arial",
})
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
OOF = TABLES / "out_of_fold_predictions.csv"
TRAINING = TABLES / "training_arms_v7.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    # 1. 从全部分组留出预测读取作图数据，不重新训练模型。
    oof = pd.read_csv(OOF)
    training = pd.read_csv(TRAINING, usecols=["candidate_arm_id", "validation_component"])
    data = oof.merge(training, on="candidate_arm_id", how="left", validate="one_to_one")
    if data.validation_component.isna().any():
        raise ValueError("Every OOF row must map to a validation component")

    observed = data.yield_Mg_ha.to_numpy(dtype=float)
    predicted = data.prediction_Mg_ha.to_numpy(dtype=float)
    if not np.isfinite(observed).all() or not np.isfinite(predicted).all():
        raise ValueError("Observed and OOF-predicted yields must be finite")

    metrics = {
        "n": len(data),
        "R2": float(r2_score(observed, predicted)),
        "RMSE": float(mean_squared_error(observed, predicted) ** 0.5),
        "MAE": float(mean_absolute_error(observed, predicted)),
        "bias": float(np.mean(predicted - observed)),
    }

    # 2. 校准趋势的置信区间按验证连通分组计算聚类稳健协方差。
    design = sm.add_constant(observed)
    fit = sm.OLS(predicted, design).fit(
        cov_type="cluster", cov_kwds={"groups": data.validation_component}
    )
    grid = np.linspace(float(observed.min()), float(observed.max()), 300)
    grid_design = np.column_stack([np.ones_like(grid), grid])
    fitted = grid_design @ fit.params
    covariance = np.asarray(fit.cov_params())
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", grid_design, covariance, grid_design), 0))

    limit_low = float(min(observed.min(), predicted.min()))
    limit_high = float(max(observed.max(), predicted.max()))
    margin = (limit_high - limit_low) * 0.035
    limits = (max(0.0, limit_low - margin), limit_high + margin)

    # 3. 保留原六边形计数与校准拟合，仅精简标题和色阶表达。
    fig, ax = plt.subplots(figsize=(5.7, 5.1))
    density = ax.hexbin(
        observed,
        predicted,
        gridsize=48,
        mincnt=1,
        cmap="cividis",
        linewidths=0,
        bins="log",
    )
    ax.fill_between(grid, fitted - 1.96 * se, fitted + 1.96 * se, color="#0072B2", alpha=0.18, linewidth=0)
    ax.plot(grid, fitted, color="#0072B2", linewidth=1.8, label="Calibration trend")
    ax.plot(limits, limits, color="#222222", linestyle="--", linewidth=1.1, label="1:1 line")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Observed yield (Mg ha$^{-1}$)", fontsize=10)
    ax.set_ylabel("Out-of-fold predicted yield (Mg ha$^{-1}$)", fontsize=10)
    ax.text(
        0.04,
        0.96,
        "\n".join(
            [
                f"n = {metrics['n']:,}",
                f"R$^2$ = {metrics['R2']:.3f}",
                f"RMSE = {metrics['RMSE']:.3f} Mg ha$^{{-1}}$",
                f"MAE = {metrics['MAE']:.3f} Mg ha$^{{-1}}$",
                f"Bias = {metrics['bias']:.3f} Mg ha$^{{-1}}$",
            ]
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
    )
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", top=False, right=False, labelsize=8.5)
    colorbar = fig.colorbar(density, ax=ax, pad=0.025, fraction=0.047)
    maximum_count = int(density.get_array().max())
    count_ticks = [value for value in (1, 2, 5, 10, 20) if value <= maximum_count]
    if count_ticks[-1] != maximum_count:
        count_ticks.append(maximum_count)
    colorbar.set_ticks(count_ticks, labels=[str(value) for value in count_ticks])
    colorbar.minorticks_off()
    colorbar.set_label("Observations per hexagon", fontsize=9)
    colorbar.ax.tick_params(labelsize=8, direction="out", length=3, width=.7, pad=3)
    fig.tight_layout()

    stem = FIGURES / "validation_oof_observed_vs_predicted"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    confidence = fit.conf_int(alpha=0.05)
    summary = pd.DataFrame(
        [
            {
                **metrics,
                "calibration_intercept": float(fit.params[0]),
                "calibration_intercept_ci_low": float(confidence[0, 0]),
                "calibration_intercept_ci_high": float(confidence[0, 1]),
                "calibration_slope": float(fit.params[1]),
                "calibration_slope_ci_low": float(confidence[1, 0]),
                "calibration_slope_ci_high": float(confidence[1, 1]),
                "calibration_covariance": "cluster robust by validation_component",
            }
        ]
    )
    summary.to_csv(TABLES / "validation_plot_metrics.csv", index=False, encoding="utf-8-sig")

    image = mpimg.imread(stem.with_suffix(".png"))
    if image.shape[0] < 1200 or image.shape[1] < 1200 or float(np.nanstd(image)) < 0.02:
        raise AssertionError("Validation PNG is blank or too small")
    result = {
        "scope": "Nested grouped out-of-fold predictions; not full-data fitted predictions",
        "rows": len(data),
        "validation_components": int(data.validation_component.nunique()),
        "oof_sha256": sha256(OOF),
        "training_sha256": sha256(TRAINING),
        "metrics": metrics,
        "calibration_covariance": "cluster robust by validation_component",
        "display": {
            "title": None,
            "y_axis": "Out-of-fold predicted yield (Mg ha^-1)",
            "colorbar_label": "Observations per hexagon",
            "colorbar_scale": "logarithmic",
            "colorbar_ticks": count_ticks,
            "maximum_hexagon_count": maximum_count,
        },
        "outputs": [
            str(stem.with_suffix(".png").relative_to(HERE)),
            str(stem.with_suffix(".pdf").relative_to(HERE)),
            str((TABLES / "validation_plot_metrics.csv").relative_to(HERE)),
        ],
    }
    (HERE / "validation_figure_manifest.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
