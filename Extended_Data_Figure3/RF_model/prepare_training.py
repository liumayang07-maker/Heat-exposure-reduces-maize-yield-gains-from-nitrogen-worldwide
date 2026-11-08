from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SHARED = [
    "reference_ID", "site_year_ID", "env_ID", "Cultivars name",
    "Irrigation_decision", "Irrigation_method", "plant_density", "Tillage",
    "SOC", "soil_total_N", "soil_pH",
]
CONTROL_KEY = SHARED + ["Label_of_control_treatment", "mean_c", "sd_c", "n_c"]
HEAT = [
    "EDH30_0_40_growing_season_fraction",
    "EDH30_40_70_growing_season_fraction",
    "EDH30_70_100_growing_season_fraction",
]


def build_training() -> pd.DataFrame:
    # 1. 用当前主表的必要字段重建操作性共享对照，不删除真实处理记录。
    source = pd.read_csv(HERE / "data/training_source.csv", float_precision="round_trip")
    assert source.data_ID.is_unique and source.Label_of_control_treatment.eq("N0").all()
    # 分期比例从完整 100 分箱重算，不以预计算比例替代其上游输入。
    bins = pd.read_csv(HERE / "data/environment_EDH30_100bins.csv", float_precision="round_trip")
    assert not bins.duplicated(["env_ID", "bin"]).any()
    assert bins[["bin_hours", "EDH30_sum"]].notna().all().all()
    assert bins[["bin_hours", "EDH30_sum"]].ge(0).all().all()
    coverage = bins.groupby("env_ID").bin.agg(["size", "nunique", "min", "max"])
    assert (coverage["size"].eq(100) & coverage["nunique"].eq(100)
            & coverage["min"].eq(1) & coverage["max"].eq(100)).all()
    bins["stage"] = np.select([bins.bin.le(40), bins.bin.le(70)], [0, 1], default=2)
    sums = bins.groupby(["env_ID", "stage"]).EDH30_sum.sum().unstack("stage")
    fractions = sums.div(sums.sum(axis=1).where(sums.sum(axis=1).gt(0)), axis=0).fillna(0)
    for stage, column in enumerate(HEAT):
        source[column] = source.env_ID.map(fractions[stage])
    assert source[HEAT].notna().all().all()
    source["candidate_control_id"] = "C" + (
        source.groupby(CONTROL_KEY, dropna=False, sort=True).ngroup() + 1
    ).astype(str).str.zfill(5)
    rows = []
    for row in source.to_dict("records"):
        rows.append(dict(
            candidate_arm_id=f"T{row['data_ID']}",
            candidate_control_id=row["candidate_control_id"],
            arm="N_treatment",
            source_data_IDs=str(row["data_ID"]),
            yield_Mg_ha=row["mean_t"],
            N_rate=row["N_rate"],
            **{key: row[key] for key in SHARED},
        ))
    for control, group in source.groupby("candidate_control_id", sort=True):
        rows.append(dict(
            candidate_arm_id=control,
            candidate_control_id=control,
            arm="control_N0_assumed",
            source_data_IDs=";".join(group.data_ID.astype(str)),
            yield_Mg_ha=group.mean_c.iloc[0],
            N_rate=0,
            **{key: group[key].iloc[0] for key in SHARED},
        ))
    training = pd.DataFrame(rows)

    # 2. 保持原有文献与完全相同坐标的连通分组，不引入空间距离阈值。
    parent = {str(ref): str(ref) for ref in source.reference_ID.unique()}

    def find(ref):
        while parent[ref] != ref:
            parent[ref] = parent[parent[ref]]
            ref = parent[ref]
        return ref

    for _, group in source.groupby(["lat", "lon"], dropna=False):
        refs = sorted(group.reference_ID.astype(str).unique())
        for ref in refs[1:]:
            parent[find(ref)] = find(refs[0])
    training["validation_component"] = training.reference_ID.astype(str).map(find)

    # 3. 每个处理或共享对照必须对应唯一的年度、气候、土壤与管理输入。
    indexed = source.set_index("data_ID")
    direct = [
        "annual_DD10_30_C_d", "annual_EDH30_C_h", "annual_precip_mm",
        "Exp_year", "country", "EDH30_total", "GDD_total", "total_rainfall",
        *HEAT,
    ]
    for feature in direct:
        values = []
        for ids in training.source_data_IDs:
            unique = indexed.loc[[int(x) for x in ids.split(";")], feature].drop_duplicates()
            if len(unique) != 1 or unique.isna().any():
                raise ValueError(f"Non-unique or missing {feature} for source IDs {ids}")
            values.append(unique.iloc[0])
        training[feature] = values
    training["irrigated"] = training.Irrigation_decision.map({"Rainfed": 0, "Irrigation": 1})
    for numerator, denominator, output in [
        ("GDD_total", "annual_DD10_30_C_d", "growing_season_GDD_fraction"),
        ("EDH30_total", "annual_EDH30_C_h", "growing_season_EDH30_fraction"),
        ("total_rainfall", "annual_precip_mm", "growing_season_precip_fraction"),
    ]:
        training[output] = training[numerator] / training[denominator].where(training[denominator] > 0)
        if numerator == "EDH30_total":
            valid_zero = training.annual_EDH30_C_h.eq(0) & training.annual_DD10_30_C_d.ge(0)
            assert training.loc[valid_zero, "EDH30_total"].eq(0).all()
            training.loc[valid_zero, output] = 0.0

    # 4. GDP 取试验前三年均值的自然对数，边界处理沿用原计算规则。
    codes = pd.read_csv(HERE / "data/country_iso3.csv").set_index("country").iso3
    gdp = pd.read_csv(HERE / "data/world_bank_gdp.csv").dropna(subset=["gdp_per_capita_ppp_constant_2021"])
    lookup = gdp.set_index(["iso3", "year"]).gdp_per_capita_ppp_constant_2021
    first_year = gdp.groupby("iso3").year.min()
    values = []
    for row in training.itertuples(index=False):
        code = codes.loc[row.country]
        prior = [lookup.get((code, max(year, first_year.loc[code])), np.nan)
                 for year in range(int(row.Exp_year) - 3, int(row.Exp_year))]
        valid = [value for value in prior if np.isfinite(value)]
        values.append(float(np.log(np.mean(valid))) if valid else np.nan)
    training["log_gdp_pc_ppp_prior3"] = values

    # 5. 按处理标识读取保存的折分，避免重新分配验证样本。
    folds = pd.read_csv(HERE / "tables/outer_folds.csv")
    training = training.merge(folds[["candidate_arm_id", "fold"]], on="candidate_arm_id", validate="one_to_one")
    columns = pd.read_csv(HERE / "tables/training_arms_v7.csv", nrows=0).columns
    training = training[list(columns)]
    assert len(training) == 6040 and training.candidate_arm_id.is_unique
    return training


def main() -> None:
    parser = argparse.ArgumentParser(description="Rebuild all RF training arms from archived source fields")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Preserve existing data: {args.output}")
    training = build_training()
    features = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))["features"]
    assert np.isfinite(training[features + ["yield_Mg_ha"]].to_numpy()).all()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    training.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"Prepared {len(training)} arms: {args.output}")


if __name__ == "__main__":
    main()
