#!/usr/bin/env python3
"""Source Data tables for every panel of Figures 2, 4 and 5.

One CSV per panel, written to <Figure>/source_data/, holding the values that are
actually drawn -- plus the attributes applied at plotting time (tissue colour,
group membership, per-box n), which otherwise live only in the plotting code.

Panel-level annotations that appear as text on a figure rather than as a point
(rho, P, n, test name) go to <Figure>/source_data/panel_statistics.csv, one row
per panel, so the per-panel tables stay tidy.

    python3 analysis/make_source_data.py            # all figures
    python3 analysis/make_source_data.py Figure4    # one figure

R panels (Fig. 2c and the Figure 5 RDS panels) are handled by
make_source_data.R, which this script calls if Rscript is available.
"""
import os
import pickle
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))


def _load(fig, name):
    with open(os.path.join(BASE, fig, 'figure_data', name + '.pickle'), 'rb') as fh:
        return pickle.load(fh)


def _out(fig):
    d = os.path.join(BASE, fig, 'source_data')
    os.makedirs(d, exist_ok=True)
    return d


def _write(df, fig, panel, index=False):
    p = os.path.join(_out(fig), f'{panel}.csv')
    df.to_csv(p, index=index)
    print(f'   {panel + ".csv":<44} {df.shape[0]:>6} rows x {df.shape[1]:>2} cols')
    return p


# --------------------------------------------------------------------------- #
# Figure 2
# --------------------------------------------------------------------------- #

def figure2():
    print('Figure2')
    stats = []

    # 2a -- one row per GTEx sample
    panels = _load('Figure2', 'fig2a_panels')
    names = _load('Figure2', 'fig2a_tissue_names')
    rows = []
    for p, display in zip(panels, names):
        rows.append(pd.DataFrame({
            'tissue': p['tissue'], 'tissue_display': display, 'color': p['color'],
            'unproductive_pct': p['pct'], 'UPF3A_tpm': p['upf3a'],
        }))
    _write(pd.concat(rows, ignore_index=True), 'Figure2', 'fig2a')

    # 2b -- one row per tissue pair, with the series it is drawn in
    sd = _load('Figure2', 'fig2b_source_data').copy()
    sd['series'] = (sd.splicing_class + '_'
                    + np.where(sd['significant_at_fdr_0.1'], 'significant', 'ns'))
    sd['plotted_x_spearman_rho'] = sd.spearman_rho
    sd['plotted_y_minus_log10_p'] = -np.log10(sd.p_value_two_sided)
    _write(sd, 'Figure2', 'fig2b')

    # 2e -- one row per sample, with the group each box is pooled into
    bx = _load('Figure2', 'fig2e_boxplot_df').copy()
    st = _load('Figure2', 'fig2e_boxplot_stats')
    n_by_tissue = dict(zip(st['tissues'], st['n_per_tissue']))
    bx['group'] = np.where(bx.tissue.isin(st['group_a']), 'brain', 'non-brain')
    bx['n_in_box'] = bx.tissue.map(n_by_tissue)
    _write(bx, 'Figure2', 'fig2e')

    t = st[st['primary']]
    stats.append({'panel': 'fig2e', 'comparison': t['label'], 'test': t['test'],
                  'statistic_name': t['stat_name'], 'statistic': t['statistic'],
                  'p_value': t['pvalue'], 'n': t.get('n_pairs')})
    pd.DataFrame(stats).to_csv(os.path.join(_out('Figure2'), 'panel_statistics.csv'), index=False)
    print(f'   {"panel_statistics.csv":<44} {len(stats):>6} rows')


# --------------------------------------------------------------------------- #
# Figure 4
# --------------------------------------------------------------------------- #

def figure4():
    print('Figure4')
    stats = []

    d = _load('Figure4', 'fig4A_data')
    _write(pd.DataFrame({
        'tissue': d['sorted_tissues'], 'tissue_display': d['tissue_names'],
        'color': d['colors'], 'n_psQTL_clusters': d['pr_counts'],
        'n_usQTL_clusters': d['up_counts'],
        'n_samples': d['n_samples'],
        'n_donors': [d['donor_counts'].get(t) for t in d['sorted_tissues']],
    }), 'Figure4', 'fig4A')

    for panel, key in (('fig4B_usQTL', 'fig4B_usQTL'), ('fig4B_psQTL', 'fig4B_psQTL'),
                       ('fig4_colocs_sc2', 'fig4_colocs_sc2')):
        fit = _load('Figure4', key + '_fit')
        _write(pd.DataFrame({'sQTL_beta': fit['x'], 'eQTL_beta': fit['y']}),
               'Figure4', panel)
        try:
            s = _load('Figure4', key + '_stats')
        except FileNotFoundError:
            continue
        stats.append({'panel': panel, 'test': s['test'], 'statistic_name': 'rho',
                      'statistic': s['rho'], 'p_value': s['pvalue'], 'n': s['n'],
                      'df': s.get('df'), 'ci95_low': s['ci95'][0], 'ci95_high': s['ci95'][1],
                      'n_donors': s.get('n_donors')})

    c = _load('Figure4', 'fig4C_data')
    sd = _load('Figure4', 'fig4C_source_data').copy()
    sd['color'] = sd.tissue.map(dict(zip(c['sorted_tissues'], c['colors'])))
    sd['tissue_display'] = sd.tissue.map(dict(zip(c['sorted_tissues'], c['tissue_names'])))
    _write(sd, 'Figure4', 'fig4C')

    # 4d -- one row per (trait, tissue) box observation
    d = _load('Figure4', 'fig4D_data')
    boxes = d['boxes']
    long = boxes.reset_index().melt(id_vars=boxes.index.name or 'index',
                                    var_name='trait', value_name='pct_loci_colocalizing')
    long = long.rename(columns={long.columns[0]: 'tissue'}).dropna()
    disp = dict(zip(boxes.columns, d['trait_names'])) if len(d['trait_names']) == len(boxes.columns) else {}
    if disp:
        long['trait_display'] = long.trait.map(disp)
    _write(long, 'Figure4', 'fig4D')

    for panel in ('fig4E_sQTL', 'fig4E_eQTL'):
        bx = _load('Figure4', panel + '_boxplot_df').copy()
        bs = _load('Figure4', panel + '_box_stats')
        grp = bs['groups']
        counts = dict(zip(grp['genotype'], grp['n']))
        bx['n_in_box'] = bx.genotype.map(counts)
        bx['color'] = 'tab:red' if panel.endswith('sQTL') else 'tab:blue'
        _write(bx, 'Figure4', panel)
        s = _load('Figure4', panel + '_stats')
        stats.append({'panel': panel, 'test': s['test'], 'statistic_name': 'beta',
                      'statistic': s['beta'], 'p_value': s['pvalue'], 'n': bs['n_total'],
                      'df': s.get('df'), 'ci95_low': s['ci95'][0], 'ci95_high': s['ci95'][1],
                      'phenotype': s['phenotype'], 'variant': s['variant']})

    _write(_load('Figure4', 'fig4E_ASB16_annot'), 'Figure4', 'fig4E_ASB16')

    f = _load('Figure4', 'fig4F_data')
    pts = f['points'].copy()
    pts.index.name = 'variant'
    pts = pts.reset_index()
    # The lead is stored by rsID and coordinates, not by the index key, so mark
    # it by matching the three plotted values it was recorded with.
    lead = f['lead']
    pts['is_lead_variant'] = (np.isclose(pts.GWAS, lead['GWAS'])
                              & np.isclose(pts.sQTL, lead['sQTL'])
                              & np.isclose(pts.eQTL, lead['eQTL']))
    if pts.is_lead_variant.sum() != 1:
        print(f"      note: lead variant matched {pts.is_lead_variant.sum()} rows")
    _write(pts, 'Figure4', 'fig4F')
    stats.append({'panel': 'fig4F', 'test': 'locuscompare (no test)', 'n': f['n_snps'],
                  'lead_variant': f['lead']['label']})

    pd.DataFrame(stats).to_csv(os.path.join(_out('Figure4'), 'panel_statistics.csv'), index=False)
    print(f'   {"panel_statistics.csv":<44} {len(stats):>6} rows')


# --------------------------------------------------------------------------- #
# Figure 5
# --------------------------------------------------------------------------- #

def figure5():
    print('Figure5')
    d = _load('Figure5', 'Fig5D_right_data')
    rows = []
    for panel, values in d['panels'].items():
        v = pd.DataFrame(values) if not isinstance(values, pd.DataFrame) else values.copy()
        v['panel'] = panel
        rows.append(v)
    out = pd.concat(rows, ignore_index=True)
    # 'geno' is the genotype the box is drawn for; counts are the per-box n
    out['n_in_box'] = out['geno'].map(d['counts'])
    out['genotype_order'] = out['geno'].map({g: i for i, g in enumerate(d['geno_order'])})
    out['color'] = np.where(out.panel.str.contains('sQTL'), 'tab:red', 'tab:blue')
    out = out.sort_values(['panel', 'genotype_order']).reset_index(drop=True)
    _write(out, 'Figure5', 'Fig5D-right')

    stats = []
    for panel, s in d['stats_by_panel'].items():
        stats.append({'panel': f'Fig5D-right ({panel})', **{k: v for k, v in s.items()
                                                            if not isinstance(v, (list, dict, np.ndarray))}})
    pd.DataFrame(stats).to_csv(os.path.join(_out('Figure5'), 'panel_statistics.csv'), index=False)
    print(f'   {"panel_statistics.csv":<44} {len(stats):>6} rows')



# --------------------------------------------------------------------------- #
# Supplementary figures
# --------------------------------------------------------------------------- #

SUP = 'SupFigures'


def supfigures():
    """Every supplementary panel that has a pickle behind it."""
    print('SupFigures')
    stats = []

    # --- panels driven by the shared per-sample table (SupFig_Fig2.ipynb) ---
    upf = _load(SUP, 'fig2group_upf_df')
    sys.path.insert(0, os.path.join(BASE, SUP))
    from SupFig_Fig2_helpers import gtex_colors

    # sup_fig4A / 4B / 4UPF1 / 4UPF2 / 4UPF3B: one point per tissue, the median
    for panel, xvar in (('sup_fig4A', 'RIN'), ('sup_fig4B', 'logUPF3A'),
                        ('sup_fig4UPF1', 'logUPF1'), ('sup_fig4UPF2', 'logUPF2'),
                        ('sup_fig4UPF3B', 'logUPF3B')):
        if xvar not in upf.columns:
            print(f'   {panel:<44} skipped: {xvar} not in the table')
            continue
        d = upf.groupby('tissue')[[xvar, 'pct']].median().reset_index()
        d.columns = ['tissue', xvar, 'unproductive_pct']
        d['color'] = ['#' + gtex_colors[t]['tissue_color_hex'] for t in d.tissue]
        d['n_samples'] = d.tissue.map(upf.groupby('tissue').size())
        _write(d, SUP, panel)

    # sup_fig7A/B/C: gene TPM per sample, by tissue
    from Figure2_helpers import TEN_TISSUES
    for panel, gene in (('sup_fig7A', 'MYOM2'), ('sup_fig7B', 'SRSF3'), ('sup_fig7C', 'DLG4')):
        if gene not in upf.columns:
            print(f'   {panel:<44} skipped: {gene} not in the table')
            continue
        d = upf.loc[upf.tissue.isin(TEN_TISSUES), ['tissue', 'entity:sample_id', gene]].copy()
        d = d.rename(columns={'entity:sample_id': 'sample_id', gene: f'{gene}_tpm'})
        d['n_in_box'] = d.tissue.map(d.groupby('tissue').size())
        _write(d.sort_values('tissue'), SUP, panel)

    # --- sup_fig_lambda: the QQ p-values, one row per test ---
    qq = _load(SUP, 'sqtl_qq')
    rows = []
    for tissue, series in qq.items():
        for name, vals in series.items():
            v = np.sort(np.asarray(vals, dtype=float))
            rows.append(pd.DataFrame({
                'tissue': tissue, 'series': name, 'nominal_p': v,
                'observed_minus_log10_p': -np.log10(v),
                'expected_minus_log10_p': -np.log10(
                    (np.arange(1, len(v) + 1) - 0.5) / len(v)),
            }))
    _write(pd.concat(rows, ignore_index=True), SUP, 'sup_fig_lambda')

    # --- the polished grids ---
    for panel, key, xname in (('sup_fig5A_polished', 'fig5A_data', 'logUPF3A'),
                              ('sup_fig5B_polished', 'fig5B_data', 'RIN')):
        d = _load(SUP, key)
        pts = pd.concat([pd.DataFrame({'tissue': e['tissue'], 'color': e['color'],
                                       xname: e['x'], 'unproductive_pct': e['y']})
                         for e in d], ignore_index=True)
        _write(pts, SUP, panel)
        for e in d:
            stats.append({'panel': panel, 'tissue': e['tissue'], 'n': e['n'],
                          'pearson_rho': e['pearson']['rho'], 'pearson_p': e['pearson']['pvalue'],
                          'spearman_rho': e['spearman']['rho'], 'spearman_p': e['spearman']['pvalue']})

    # cumulative curves, one row per drawn point
    d = _load(SUP, 'up_by_expression_data')
    rows = []
    for e in d:
        for q, (x, y) in enumerate(e['curves'], start=1):
            rows.append(pd.DataFrame({'tissue': e['tissue'], 'color': e['color'],
                                      'expression_quintile': q,
                                      'log10_unproductive_pct': np.asarray(x),
                                      'cumulative_fraction': np.asarray(y)}))
        stats.append({'panel': 'sup_fig_UP_by_expression_polished', 'tissue': e['tissue'],
                      'n': e['n'], 'pearson_rho': e['pearson']['rho'],
                      'pearson_p': e['pearson']['pvalue'],
                      'spearman_rho': e['spearman']['rho'], 'spearman_p': e['spearman']['pvalue']})
    _write(pd.concat(rows, ignore_index=True), SUP, 'sup_fig_UP_by_expression_polished')

    # 49 x 49 correlation matrix, long, both metrics with their clustering order
    c = _load(SUP, 'corr_matrix_data')
    tissues = c['tissues']
    rows = []
    for metric in ('pearson', 'spearman'):
        M, order = c[metric]['matrix'], list(c[metric]['order'])
        pos = {t: order.index(i) for i, t in enumerate(tissues)}
        for i, t1 in enumerate(tissues):
            for j, t2 in enumerate(tissues):
                rows.append((metric, t1, t2, M[i, j], pos[t1], pos[t2]))
    _write(pd.DataFrame(rows, columns=['metric', 'splicing_tissue_row',
                                       'expression_tissue_column', 'rho',
                                       'row_position_in_figure',
                                       'column_position_in_figure']),
           SUP, 'sup_fig_corr_across_tissues_polished')

    pd.DataFrame(stats).to_csv(os.path.join(_out(SUP), 'panel_statistics.csv'), index=False)
    print(f'   {"panel_statistics.csv":<44} {len(stats):>6} rows')



# --------------------------------------------------------------------------- #
# Figure 1
# --------------------------------------------------------------------------- #

def figure1():
    print('Figure1')
    # 1c: the BED12 junction tracks behind the browser screenshot, one row per
    # junction drawn, with the class the track is coloured by
    tracks = _load('Figure1', 'fig1c_tracks')
    rows = []
    for q, df in tracks.items():
        d = df.copy()
        d.insert(0, 'usage_quartile', q)
        rows.append(d)
    _write(pd.concat(rows, ignore_index=True), 'Figure1', 'fig1c_bed12_tracks')

    _write(_load('Figure1', 'fig1d_data'), 'Figure1', 'fig1d')
    _write(_load('Figure1', 'fig1e_data'), 'Figure1', 'fig1e')

    # 1g: the drawn ECDF, plus the per-junction log2FCs behind it
    g = _load('Figure1', 'fig1g_data')
    _write(pd.concat([pd.DataFrame({
        'comparison': s['comparison_label'], 'category': s['category_label'],
        'color': s['color'], 'log2fc_x': s['ecdf_x'], 'cumulative_fraction': s['ecdf_y'],
    }) for s in g], ignore_index=True), 'Figure1', 'fig1g')
    _write(pd.concat([pd.DataFrame({
        'comparison': s['comparison_label'], 'category': s['category_label'],
        'log2fc': s['log2fc'],
    }) for s in g], ignore_index=True), 'Figure1', 'fig1g_log2fc_values')

    _write(_load('Figure1', 'fig1h_data'), 'Figure1', 'fig1h')

    stats = [{'panel': 'fig1g', 'comparison': s['comparison_label'],
              'category': s['category_label'], 'n': s['n'],
              'median_log2fc': s['median_log2fc']} for s in g]
    h = _load('Figure1', 'fig1h_data')
    for _, r in h.iterrows():
        stats.append({'panel': 'fig1h', 'comparison': r['rule'], 'category': r['class'],
                      'n': (r['n_low'] + r['n_high']) if pd.notna(r['n_low']) else None,
                      'test': 'Mann-Whitney U, one-sided',
                      'statistic': r.get('U'), 'p_value': r['p_value'],
                      'delta_log2fd': r['delta_log2fd']})
    pd.DataFrame(stats).to_csv(os.path.join(_out('Figure1'), 'panel_statistics.csv'), index=False)
    print(f'   {"panel_statistics.csv":<44} {len(stats):>6} rows')


def main(argv):
    todo = argv or ['Figure1', 'Figure2', 'Figure4', 'Figure5', 'SupFigures']
    fns = {'Figure1': figure1, 'Figure2': figure2, 'Figure4': figure4,
           'Figure5': figure5, 'SupFigures': supfigures}
    for f in todo:
        fns[f]()
    r = os.path.join(BASE, 'make_source_data.R')
    rscript = shutil.which('Rscript') or '/software/R-4.1.0-el8-x86_64/bin/Rscript'
    if os.path.exists(r) and os.path.exists(rscript):
        print('\nR panels:')
        env = dict(os.environ)
        # the cluster R needs these on the library path to start at all
        env['LD_LIBRARY_PATH'] = ':'.join(filter(None, [
            '/software/gcc-13.2.0-el8-x86_64/lib64',
            '/software/openblas-0.3.13-el8-x86_64/lib',
            env.get('LD_LIBRARY_PATH')]))
        subprocess.run([rscript, r], check=False, env=env)
    elif os.path.exists(r):
        print('\nR panels: skipped, Rscript not found')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
