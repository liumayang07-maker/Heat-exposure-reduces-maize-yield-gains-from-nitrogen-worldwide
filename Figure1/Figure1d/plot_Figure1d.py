"""Reproduce both Figure1d plots and an English-only Source Data workbook.

Original sources: scripts 07 (row Z-scores/Ward.D2) and 09 (OLS/HAC/BH).
"""
import argparse
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
OUTPUTS = ('Figure1d_heatmap.pdf', 'Figure1d_trend.pdf', 'Figure1d_Source_Data.xlsx')


def load_module(path, name):
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fit_models(data, styles, cluster_order):
    # 柱图只用 1986–2025 年原始暴露，不拟合热图 Z-score。
    annual = data.loc[data.season_year.between(1986, 2025)]
    values = annual.pivot(index='season_year', columns='region_id', values='area_weighted_mean')
    y = values.loc[list(range(1986, 2026)), styles.region_id].to_numpy()
    assert len(annual) == 600 and y.shape == (40, 15) and np.isfinite(y).all() and (y >= 0).all()
    assert (annual.area_source_year == np.minimum(annual.season_year, 2020)).all()
    shared = load_module(ROOT.parent / 'Figure1c/prepare_grid_EDH30_trends.py', 'grid_hac')
    slope, se, p = shared.hac_test(y)
    design = sm.add_constant(np.arange(40) - 19.5)
    parameters, covariance = [], []
    for i, rid in enumerate(styles.region_id):
        fit = sm.OLS(y[:, i], design).fit(
            cov_type='HAC',
            cov_kwds={
                'maxlags': 3,
                'kernel': 'bartlett',
                'use_correction': True,
            },
            use_t=True,
        )
        np.testing.assert_allclose([slope[i], se[i], p[i]],
            [fit.params[1]*10, fit.bse[1]*10, fit.pvalues[1]], rtol=1e-7, atol=1e-9)
        parameters.append(dict(region_id=rid, year_center=2005.5, intercept=fit.params[0],
                               slope_per_year=fit.params[1], residual_sigma=np.sqrt(fit.scale),
                               n_years=40, df=38, hac_lags=3))
        terms = ['intercept', 'year_c']
        for row in range(2):
            for col in range(2):
                covariance.append(dict(region_id=rid, term_row=terms[row], term_col=terms[col],
                                       covariance=fit.cov_params()[row, col]))
    q = multipletests(p, method='fdr_bh')[1]
    rename = {'CHN_NE': 'Northeast China', 'CHN_NCP': 'North China Plain',
              'CHN_SW': 'Southwest China', 'CHN_NW': 'Northwest China'}
    names = [rename.get(rid, name) for rid, name in zip(styles.region_id, styles.region_name)]
    ci = t.ppf(.975, 38) * se
    result = pd.DataFrame({'region_id': styles.region_id.to_numpy(), 'region_name': names,
        'start_year': 1986, 'end_year': 2025, 'n_years': 40,
        'slope_per_decade': slope, 'HAC_se_per_decade': se,
        'ci95_lower': slope-ci, 'ci95_upper': slope+ci, 'p_raw': p, 'p_BH': q,
        'significant_BH_005': q < .05})
    result['plot_order'] = result.region_id.map({rid: i+1 for i, rid in enumerate(cluster_order)})
    result = result.sort_values('plot_order').reset_index(drop=True)
    return result, pd.DataFrame(parameters), pd.DataFrame(covariance)


def plot_trend(result, output):
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': .6,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(135/25.4, 78/25.4))
    x = np.arange(15)
    plot_slope = result['slope_per_decade'].to_numpy()
    plot_ci = (result['ci95_upper'] - result['slope_per_decade']).to_numpy()
    plot_q = result['p_BH'].to_numpy()
    plot_names = result['region_name'].to_list()
    bars = ax.bar(x, plot_slope, width=.7, color=np.where(plot_slope >= 0, '#C74A43', '#3679AC'),
                  edgecolor='#333333', linewidth=.35)
    ax.errorbar(x, plot_slope, yerr=plot_ci, fmt='none', ecolor='#303030', elinewidth=.65,
                capsize=2, capthick=.65)
    ax.axhline(0, color='#333333', linewidth=.6)
    span = max((plot_slope+plot_ci).max(), 0) - min((plot_slope-plot_ci).min(), 0)
    for i in range(15):
        stars = '***' if plot_q[i] < .001 else '**' if plot_q[i] < .01 else '*' if plot_q[i] < .05 else ''
        pos = plot_slope[i]+plot_ci[i]+span*.025 if plot_slope[i] >= 0 else plot_slope[i]-plot_ci[i]-span*.025
        ax.text(i, pos, stars, ha='center', va='bottom' if plot_slope[i] >= 0 else 'top', fontsize=8)
    ax.set_ylim(min((plot_slope-plot_ci).min(), 0)-span*.13, max((plot_slope+plot_ci).max(), 0)+span*.13)
    ax.set_xticks(x, plot_names, rotation=65, ha='right', fontsize=6)
    ax.set_ylabel('Summer EDH30 trend (\N{DEGREE SIGN}C\N{MIDDLE DOT}h per decade)', fontsize=7)
    ax.tick_params(width=.5, length=2)
    ax.set_xlim(-.6, 14.6)
    fig.subplots_adjust(left=.14, right=.98, top=.97, bottom=.43)
    assert len(bars) == 15
    fig.savefig(output, dpi=600, facecolor='white')
    plt.close(fig)


def main():
    # 1. 检查输入及覆盖选项，原始数据始终只读。
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'region_summer_1985_2025.csv')
    parser.add_argument('--output-dir', type=Path, default=ROOT)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not args.overwrite and any((args.output_dir / name).exists() for name in OUTPUTS):
        raise FileExistsError('Outputs exist; use --overwrite explicitly.')
    data = pd.read_csv(args.input)
    styles = pd.read_csv(ROOT / 'regions.csv').sort_values('order')
    assert len(data) == 615 and not data.duplicated(['region_id', 'season_year']).any()
    assert len(styles) == 15 and not styles.region_id.duplicated().any()
    with tempfile.TemporaryDirectory(prefix='figure1d_plot_') as temporary:
        temp = Path(temporary)
        # 2. 先生成热图及聚类顺序，再按同一顺序拟合和绘制趋势。
        subprocess.run(
            [
                args.rscript,
                str(ROOT / 'plot_heatmap.R'),
                str(args.input.resolve()),
                str(ROOT / 'regions.csv'),
                str(temp),
            ],
            check=True,
        )
        heatmap = pd.read_csv(temp / 'Heatmap_data.csv')
        cluster_order = heatmap.sort_values('row_order').region_id.drop_duplicates().to_list()
        assert len(cluster_order) == 15 and set(cluster_order) == set(styles.region_id)
        result, parameters, covariance = fit_models(data, styles, cluster_order)
        plot_trend(result, temp / 'Figure1d_trend.pdf')
        writer = load_module(ROOT.parent / 'build_source_data.py', 'source_workbook')
        # 3. 只导出五类英文数据，临时表不进入投稿目录。
        writer.write_workbook(temp / 'Figure1d_Source_Data.xlsx', {
            'Annual_data': data, 'Heatmap_data': heatmap, 'Trend_models': result,
            'Model_parameters': parameters, 'Coefficient_covariance': covariance})
        args.output_dir.mkdir(parents=True, exist_ok=True)
        # 所有计算和工作簿回读核对成功后，再写入最终目录。
        for name in OUTPUTS:
            shutil.copy2(temp / name, args.output_dir / name)
    print('Verified: 615 heatmap cells; 15 OLS/HAC models; five Source Data sheets.')


if __name__ == '__main__':
    main()
