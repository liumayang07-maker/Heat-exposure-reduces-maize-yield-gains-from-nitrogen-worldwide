"""Compact V7 SHAP figure with restricted cubic spline (RCS) trends."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import dmatrix, build_design_matrices

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT.parent
FIGURE_SIZE_INCHES = (297 / 25.4, 210 / 25.4)
IMPORTANCE_LABELS = {
    'growing_season_GDD_fraction': 'GDD ratio (season/annual)',
    'growing_season_EDH30_fraction': 'EDH30 ratio (season/annual)',
    'growing_season_precip_fraction': 'Rainfall ratio (season/annual)',
    'EDH30_0_40_growing_season_fraction': 'Early-season EDH30 proportion',
    'EDH30_40_70_growing_season_fraction': 'Mid-season EDH30 proportion',
    'EDH30_70_100_growing_season_fraction': 'Late-season EDH30 proportion',
}
LABELS = {
    'N_rate': ('N rate', r'N rate (kg N ha$^{-1}$)'),
    'plant_density': ('Plant density', r'Density (plants m$^{-2}$)'),
    'annual_DD10_30_C_d': ('Annual DD', 'Annual DD$_{10-30}$ (\N{DEGREE SIGN}C d)'),
    'annual_EDH30_C_h': ('Annual EDH30', 'Annual EDH$_{30}$ (\N{DEGREE SIGN}C h)'),
    'annual_precip_mm': ('Annual precipitation', 'Annual precipitation (mm)'),
    'SOC': ('Soil organic C', r'SOC (g kg$^{-1}$)'),
    'soil_total_N': ('Soil total N', r'Soil total N (g kg$^{-1}$)'),
    'soil_pH': ('Soil pH', 'Soil pH'),
    'irrigated': ('Irrigation', 'Water management'),
    'growing_season_GDD_fraction': ('Growing-season-to-\nannual GDD ratio', 'Growing-season GDD\n/ annual DD'),
    'growing_season_EDH30_fraction': ('Growing-season-to-\nannual EDH30 ratio', 'Growing-season EDH30\n/ annual EDH30'),
    'growing_season_precip_fraction': ('Growing-season-to-\nannual precipitation ratio', 'Growing-season rainfall\n/ annual rainfall'),
    'EDH30_0_40_growing_season_fraction': ('Early-season EDH30\nproportion', 'EDH30 (0-40% GDD)\n/ growing-season EDH30'),
    'EDH30_40_70_growing_season_fraction': ('Mid-season EDH30\nproportion', 'EDH30 (40-70% GDD)\n/ growing-season EDH30'),
    'EDH30_70_100_growing_season_fraction': ('Late-season EDH30\nproportion', 'EDH30 (70-100% GDD)\n/ growing-season EDH30'),
    'Exp_year': ('Experiment year', 'Experiment year'),
    'log_gdp_pc_ppp_prior3': ('Prior-3-year GDPpc', 'ln prior-3-year GDPpc (PPP)'),
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inputs():
    # 1. 共用 RF_model 中的模型和完整训练表，不生成重复输入快照。
    manifest = json.loads((ROOT / 'manifest.json').read_text(encoding='utf-8'))
    assert sha(MODEL / 'rf_model.joblib') == manifest['model_sha256']
    assert sha(MODEL / 'tables/training_arms_v7.csv') == manifest['training_table_sha256']
    for name, expected in manifest['input_files'].items():
        assert sha(ROOT / name) == expected, name
    return manifest


def load(manifest):
    directory = ROOT / 'tables'
    train = pd.read_csv(MODEL / 'tables/training_arms_v7.csv')
    raw = pd.read_csv(directory / 'shap_values.csv')
    features = manifest['features']
    # Exp_year appears twice: metadata first, SHAP contribution in the feature block.
    with (directory / 'shap_values.csv').open(encoding='utf-8-sig', newline='') as stream:
        header = next(csv.reader(stream))
    assert header[6:] == features
    values = raw.iloc[:, 6:].copy()
    values.columns = features
    assert train.candidate_arm_id.is_unique and raw.candidate_arm_id.is_unique
    assert raw.candidate_arm_id.equals(train.candidate_arm_id)
    assert raw.reference_ID.equals(train.reference_ID)
    assert raw.iloc[:, 5].equals(train.Exp_year)
    x = train[features].astype(float)
    assert len(x) == manifest['rows'] and np.isfinite(x.to_numpy()).all()
    assert np.isfinite(values.to_numpy()).all()
    importance = pd.read_csv(directory / 'importance.csv').sort_values('rank')
    recalculated = values.abs().mean().reindex(importance.feature).to_numpy()
    assert np.allclose(recalculated, importance.mean_abs_SHAP, rtol=1e-10)
    assert np.allclose(recalculated / recalculated.sum(), importance.contribution_share)
    prediction = pd.read_csv(directory / 'prediction_decomposition.csv')
    assert prediction.candidate_arm_id.equals(train.candidate_arm_id)
    error = np.max(np.abs(values.sum(axis=1) + manifest['base_value'] - prediction.prediction_Mg_ha))
    assert error < 1e-6
    return x, values, importance, float(error), train.reference_ID


def fit_rcs_curves(x, values, groups):
    curves, methods = [], []
    quantiles = [.05, .275, .50, .725, .95]
    for feature in x.columns:
        xx = x[feature].to_numpy()
        yy = values[feature].to_numpy()
        if feature == 'irrigated':
            assert np.array_equal(np.unique(xx), [0, 1])
            knots = []
            design = dmatrix('1 + x', {'x': xx}, return_type='dataframe')
            grid = np.array([0., 1.])
            method = 'Binary category means'
        else:
            knots = np.unique(np.quantile(xx, quantiles)).tolist()
            if len(knots) < 3:
                raise ValueError(f'Insufficient distinct RCS knots: {feature}')
            # Natural cubic regression spline spans the RCS space, linear in both tails.
            design = dmatrix(
                '1 + cr(x, knots=inner, lower_bound=lower, upper_bound=upper, constraints="center")',
                {'x': xx, 'inner': knots[1:-1], 'lower': knots[0], 'upper': knots[-1]},
                return_type='dataframe')
            grid = np.unique(np.r_[np.linspace(xx.min(), xx.max(), 240), knots])
            method = 'Gaussian GLM restricted cubic spline (natural cubic basis)'
        assert np.linalg.matrix_rank(design) == design.shape[1], feature
        fit = sm.GLM(yy, design, family=sm.families.Gaussian()).fit(
            cov_type='cluster', cov_kwds={'groups': groups.to_numpy()})
        assert fit.converged, feature
        parameters = {'inner': knots[1:-1], 'lower': knots[0], 'upper': knots[-1]} if knots else {}
        new = build_design_matrices([design.design_info], {'x': grid, **parameters})[0]
        predicted = fit.get_prediction(new).summary_frame(alpha=.05)
        curve = pd.DataFrame({'feature': feature, 'x': grid,
                              'mean': predicted['mean'].to_numpy(),
                              'lower_95': predicted.mean_ci_lower.to_numpy(),
                              'upper_95': predicted.mean_ci_upper.to_numpy()})
        assert np.isfinite(curve.iloc[:, 1:].to_numpy()).all(), feature
        assert curve.lower_95.le(curve['mean']).all() and curve.upper_95.ge(curve['mean']).all()
        tail_errors = []
        if knots:
            step = (knots[-1] - knots[0]) * .1
            for points in (knots[0] - step * np.array([3, 2, 1]),
                           knots[-1] + step * np.array([1, 2, 3])):
                tail_design = build_design_matrices([design.design_info], {'x': points, **parameters})[0]
                # Verify the basis itself, not just one set of fitted coefficients.
                error = float(np.max(np.abs(np.diff(np.asarray(tail_design), n=2, axis=0))))
                assert error < 1e-8, (feature, error)
                tail_errors.append(error)
        else:
            assert np.allclose(curve['mean'], [yy[xx == level].mean() for level in (0, 1)])
        curves.append(curve)
        methods.append({'feature': feature, 'method': method, 'knots': knots,
                        'design_columns': design.shape[1], 'n': len(xx),
                        'clusters': int(groups.nunique()), 'converged': bool(fit.converged),
                        'tail_basis_second_difference': tail_errors})
    combined = pd.concat(curves, ignore_index=True)
    combined.to_csv(ROOT / 'outputs/tables/dependence_curves_RCS.csv', index=False)
    metadata = {'family': 'Gaussian', 'link': 'identity', 'knot_quantiles': quantiles,
                'duplicate_knots': 'Removed; no sample removal',
                'covariance': 'cluster robust by reference_ID',
                'interval': 'Pointwise 95% confidence interval for the mean trend',
                'panels': methods}
    (ROOT / 'outputs/tables/dependence_method_RCS.json').write_text(
        json.dumps(metadata, indent=2), encoding='utf-8')
    return combined


def clean(ax):
    ax.tick_params(labelsize=7.7, length=2, width=.6, pad=2)
    for spine in ax.spines.values():
        spine.set_color('#9b9b9b')
        spine.set_linewidth(.6)


def importance_panel(fig, slot, x, values, importance):
    ax = fig.add_subplot(slot)
    bars = ax.twiny()
    bars.barh(np.arange(17), importance.contribution_share * 100,
              color='#9a9a9a', alpha=.23, height=.68, zorder=0)
    bars.set_xlim(0, 40)
    bars.set_xticks([0, 10, 20, 30])
    bars.set_xlabel('Mean |SHAP| share (%)', fontsize=8.5, labelpad=4)
    bars.tick_params(labelsize=7.7, length=2)
    ax.set_zorder(2)
    ax.patch.set_visible(False)
    rng = np.random.default_rng(20260914)
    cmap = LinearSegmentedColormap.from_list('feature_value', ['#176caa', '#bebad3', '#e27433'])
    for row, feature in enumerate(importance.feature):
        y = values[feature].to_numpy()
        xx = x[feature].to_numpy()
        counts, edges = np.histogram(y, bins=55)
        bins = np.minimum(np.searchsorted(edges, y, side='right') - 1, len(counts) - 1)
        jitter = rng.uniform(-1, 1, len(y)) * .31 * np.sqrt(counts[bins] / counts.max())
        colors = (xx - xx.min()) / (xx.max() - xx.min())
        order = rng.permutation(len(y))
        ax.scatter(y[order], row + jitter[order], c=colors[order], cmap=cmap, vmin=0, vmax=1,
                   s=2.2, alpha=.5, linewidths=0, rasterized=True)
        ax.text(.985, row, f'{importance.iloc[row].contribution_share * 100:.1f}%',
                transform=ax.get_yaxis_transform(), va='center', ha='right', fontsize=7.5,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .68, 'pad': .3})
    ax.axvline(0, color='#888888', lw=.6, zorder=0)
    ax.set_yticks(range(17), [IMPORTANCE_LABELS.get(f, LABELS[f][0]) for f in importance.feature])
    ax.set_ylim(16.7, -.7)
    ax.set_xlim(-4, 3.7)
    ax.set_xticks([-4, -2, 0, 2])
    ax.set_xlabel('SHAP value (Mg/ha)', fontsize=8.5)
    ax.set_title('a', loc='left', fontsize=11, fontweight='bold', pad=32)
    clean(ax)
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=cmap), ax=ax,
                      orientation='horizontal', fraction=.035, pad=.12, aspect=24)
    cb.set_ticks([0, 1], labels=['Low', 'High'])
    cb.ax.tick_params(labelsize=7.5, length=0)
    cb.set_label('Feature value (within variable)', fontsize=7.5, labelpad=1)


def dependence(fig, slot, feature, x, y, curve, letter, histogram_rows):
    gs = slot.subgridspec(2, 1, height_ratios=[4.5, 1], hspace=.025)
    ax = fig.add_subplot(gs[0])
    hist = fig.add_subplot(gs[1], sharex=ax)
    rng = np.random.default_rng(20260914)
    binary = feature == 'irrigated'
    xp = x + rng.uniform(-.06, .06, len(x)) if binary else x
    ax.scatter(xp, y, s=3.5, color='#0b5c8f', alpha=.15, linewidths=0, rasterized=True)
    if binary:
        center = curve['mean'].to_numpy()
        ax.errorbar(curve.x, center, yerr=np.vstack([center - curve.lower_95, curve.upper_95 - center]),
                    fmt='o', color='#d38c24', ms=4, lw=1.3, capsize=3, zorder=5)
        bins = np.array([-.5, .5, 1.5])
        hist.set_xticks([0, 1], ['Rainfed', 'Irrigated'])
    else:
        ax.fill_between(curve.x, curve.lower_95, curve.upper_95, color='#e1ac49', alpha=.35, linewidth=0, zorder=3)
        ax.plot(curve.x, curve['mean'], color='#d08a1e', lw=1.5, zorder=4)
        bins = np.linspace(x.min(), x.max(), 33)
        hist.xaxis.set_major_locator(MaxNLocator(3, integer=feature == 'Exp_year'))
        hist.ticklabel_format(axis='x', style='plain', useOffset=False)
    counts, edges = np.histogram(x, bins=bins)
    hist.bar(edges[:-1], counts, width=np.diff(edges), align='edge', color='#707070', edgecolor='white', linewidth=.18)
    for i, count in enumerate(counts):
        histogram_rows.append({'feature': feature, 'left': edges[i], 'right': edges[i+1], 'count': int(count)})
    ax.axhline(0, color='#969696', lw=.65)
    ax.grid(color='#eeeeee', lw=.4)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_ylabel('SHAP (Mg/ha)', fontsize=8.2, labelpad=2)
    ax.tick_params(labelbottom=False)
    ax.set_title(letter, loc='left', fontsize=9.5, fontweight='bold', pad=4)
    hist.set_xlabel(LABELS[feature][1], fontsize=8.2, labelpad=2)
    hist.set_yticks([])
    hist.set_ylabel('Count', fontsize=7.2, labelpad=2)
    hist.yaxis.set_label_coords(-.18, .5)
    hist.set_xlim((-.4, 1.4) if binary else (float(x.min()), float(x.max())))
    clean(ax)
    clean(hist)


def main():
    global ROOT
    parser = argparse.ArgumentParser(description='Plot the archived full-training SHAP analysis')
    parser.add_argument('--data-dir', type=Path, default=ROOT)
    args = parser.parse_args()
    ROOT = args.data_dir.resolve()
    # 2. 保持原拟合与面板顺序，以 A4 横向紧凑排版。
    for folder in ('outputs/figures', 'outputs/tables'):
        (ROOT / folder).mkdir(parents=True, exist_ok=True)
    manifest = inputs()
    x, values, importance, error, groups = load(manifest)
    curves = fit_rcs_curves(x, values, groups)
    studies = groups.nunique()
    # 数学标注也使用 Arial，避免上下标及单位回退到 DejaVu Sans。
    plt.rcParams.update({
        'font.family': 'Arial',
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'font.size': 8,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
        'mathtext.default': 'regular',
        'mathtext.fallback': None,
    })
    fig = plt.figure(figsize=FIGURE_SIZE_INCHES)
    grid = fig.add_gridspec(
        4, 5,
        left=.14,
        right=.99,
        bottom=.085,
        top=.90,
        wspace=.30,
        hspace=.48,
    )
    importance_panel(fig, grid[:2, 0], x, values, importance)
    slots = [(r, c) for r in range(4) for c in range(5) if (r, c) not in [(0, 0), (1, 0), (3, 4)]]
    histogram_rows = []
    panels = []
    for idx, (feature, (r, c)) in enumerate(zip(importance.feature, slots)):
        curve = curves[curves.feature.eq(feature)]
        letter = chr(98 + idx)
        dependence(fig, grid[r, c], feature, x[feature].to_numpy(), values[feature].to_numpy(), curve, letter, histogram_rows)
        panels.append({'panel': letter, 'feature': feature, 'rank': idx + 1})
    legend = fig.add_subplot(grid[3, 4])
    legend.set_axis_off()
    legend.legend(handles=[Line2D([], [], marker='o', linestyle='none', ms=4, color='#0b5c8f', label='Sample SHAP'),
                           Line2D([], [], lw=1.5, color='#d08a1e', label='RCS fit'),
                           Patch(facecolor='#e1ac49', alpha=.35, label='Pointwise 95% CI'),
                           Patch(facecolor='#707070', label='Sample distribution')],
                  loc='upper left', frameon=False, fontsize=8.5, handlelength=1.8, borderaxespad=0)
    files = []
    for extension in ('png', 'pdf', 'svg'):
        path = ROOT / f'outputs/figures/Extended_Figure_RF_V7_SHAP_compact.{extension}'
        fig.savefig(path, dpi=300, facecolor='white')
        files.append({'path': path.relative_to(ROOT).as_posix(), 'sha256': sha(path)})
    plt.close(fig)
    pd.DataFrame(histogram_rows).to_csv(ROOT / 'outputs/tables/histogram_counts.csv', index=False)
    pd.DataFrame(panels).to_csv(ROOT / 'outputs/tables/panel_features.csv', index=False)
    report = {'status': 'PASS', 'rows': len(x), 'studies': int(studies), 'features': 17,
              'maximum_additivity_error': error, 'importance_recalculation': 'PASS',
              'duplicate_Exp_year_handling': 'Metadata first; feature block read positionally and verified against header',
              'model_sha256': manifest['model_sha256'], 'script_sha256': sha(Path(__file__)),
              'input_manifest_sha256': sha(ROOT / 'manifest.json'),
              'trend': 'Gaussian GLM RCS; study-clustered pointwise mean CI',
              'RCS_tail_linearity_check': 'PASS', 'trend_refitted': True,
              'RF_refitted': False, 'SHAP_recalculated': False,
              'rcs_curves_sha256': sha(ROOT / 'outputs/tables/dependence_curves_RCS.csv'),
              'rcs_method_sha256': sha(ROOT / 'outputs/tables/dependence_method_RCS.json'),
              'no_axis_trimming': True, 'figure_outputs': files,
              'layout': {'size_mm': [297, 210], 'orientation': 'landscape',
                         'title_font_pt': 9.5, 'axis_label_font_pt': 8.2,
                         'tick_font_pt': 7.7, 'wspace': .30, 'hspace': .48}}
    (ROOT / 'outputs/verification.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(f'PASS: {len(x)} samples, 17 predictors, SHAP reconstruction error {error:.3g}; PNG/PDF/SVG exported.')


if __name__ == '__main__':
    main()
