"""Data for the three polished supplementary panels.

  sup_fig5A_polished            per-tissue eQTL-style scatter, logUPF3A vs %UP
  sup_fig5B_polished            the same against RIN
  sup_fig_UP_by_expression_polished
                                per-tissue cumulative %UP by expression quintile,
                                each tissue ranked by its OWN expression

All three are 7x7 grids over the 49 GTEx tissues. Note that 7 x 7 = 49 exactly,
so every panel holds a tissue and none is spare -- see plot helpers.

run_all() does the expensive reading once and pickles plot-ready objects into
figure_data/; the notebook plots from those.
"""
import os
import pickle

import numpy as np
import pandas as pd

BASE = '/project/yangili1/cfbuenabadn/leafcutter2_paper'
sys_path_added = os.path.join(BASE, 'analysis', 'Figure2')

import sys
if sys_path_added not in sys.path:
    sys.path.insert(0, sys_path_added)

# gtex_colors and the per-sample unproductive-percentage / TPM chain are already
# implemented for Figure 2; reuse them rather than keeping a second copy.
from Figure2_helpers import (gtex_colors, load_unproductive_pct, load_gene_tpm,
                             build_upf_df)

JUNCS_DIR = f'{BASE}/code/results/juncs_with_gene'
PHENO_DIR = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx'
RIN_TABLE = (f'{BASE}/code/resources/GTEx/'
             'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
MEDIAN_TPM = ('/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex_tables/'
              'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz')

PLOTS_DIR = f'{BASE}/analysis/SupFigures/plots'


def all_tissues():
    """The 49 GTEx tissues, sorted -- the panel order of every grid here."""
    return sorted(os.listdir(PHENO_DIR))


# --------------------------------------------------------------------------- #
# sup_fig5A / sup_fig5B -- per-sample UPF3A, RIN and unproductive percentage
# --------------------------------------------------------------------------- #

def load_rin(rin_table=RIN_TABLE):
    m = pd.read_csv(rin_table, sep='\t')[['SAMPID', 'SMRIN']]
    m.columns = ['entity:sample_id', 'RIN']
    return m


def make_upf_df():
    """One row per GTEx sample: tissue, %unproductive reads, UPF3A TPM, RIN."""
    pct_df = load_unproductive_pct()
    gene_df = load_gene_tpm(genes=['UPF3A'])
    upf_df = build_upf_df(gene_df, pct_df)
    upf_df = upf_df.merge(load_rin(), on='entity:sample_id')
    upf_df['logUPF3A'] = np.log10(upf_df.UPF3A + 1)
    return upf_df


def make_scatter_grid_data(upf_df, xvar):
    """Per-tissue points, RLM fit line and Spearman stats for one 7x7 grid."""
    import statsmodels.api as sm
    from scipy.stats import spearmanr

    out = []
    for tissue, df in upf_df.sort_values('tissue').groupby('tissue'):
        X = sm.add_constant(df[[xvar]])
        res = sm.RLM(df.pct, X).fit()
        slope = res.params.loc[xvar]
        const = res.params.loc['const']
        rho, pval = spearmanr(df[xvar], df.pct)
        x = np.linspace(df[xvar].min(), df[xvar].max(), 100)
        out.append({
            'tissue': tissue,
            'color': '#' + gtex_colors[tissue]['tissue_color_hex'],
            'x': np.array(df[xvar]), 'y': np.array(df.pct),
            'fit_x': x, 'fit_y': slope * x + const,
            'rho': float(rho), 'pvalue': float(pval), 'n': int(df.shape[0]),
        })
    return out


# --------------------------------------------------------------------------- #
# sup_fig_UP_by_expression -- cumulative %UP by expression quintile
# --------------------------------------------------------------------------- #

def _transform(s):
    import re
    return re.sub(r'\)', '', re.sub(r'\(', '_', s.replace(' ', '')))


def load_median_expression(path=MEDIAN_TPM):
    me = pd.read_csv(path, sep='\t', skiprows=2)
    me.columns = [_transform(c) for c in me.columns]
    me['gene'] = me.Name.apply(lambda x: x.split('.')[0])
    return me


def get_UP_table(tissue, juncs_dir=JUNCS_DIR):
    """Per-gene unproductive / productive junction-read percentages."""
    df = pd.read_csv(f'{juncs_dir}/{tissue}.tab.gz', sep='\t')
    s = df.groupby(['ensembl_id', 'annot']).counts.sum().reset_index()

    total = s.groupby('ensembl_id').counts.sum().reset_index()
    total.columns = ['ensembl_id', 'total_counts']
    pr = s.loc[s.annot == 'PR', ['ensembl_id', 'counts']]
    pr.columns = ['ensembl_id', 'PR_counts']
    up = s.loc[s.annot == 'UP', ['ensembl_id', 'counts']]
    up.columns = ['ensembl_id', 'UP_counts']
    nonpr = s.loc[s.annot != 'PR'].groupby('ensembl_id').counts.sum().reset_index()
    nonpr.columns = ['ensembl_id', 'nonPR_counts']

    m = (total.merge(pr, on='ensembl_id')
              .merge(up, on='ensembl_id')
              .merge(nonpr, on='ensembl_id'))
    m['UP_percentage'] = 100 * m.UP_counts / m.total_counts
    m['nonPR_percentage'] = 100 * m.nonPR_counts / m.total_counts
    return m


def get_splicing_expression(up_table, median_expression, tissue2):
    shared = pd.Index(up_table.ensembl_id).intersection(pd.Index(median_expression.gene))
    med = median_expression.loc[median_expression.gene.isin(shared), ['gene', tissue2]]
    med.columns = ['gene', 'expression']
    se = up_table.merge(med, left_on='ensembl_id', right_on='gene').dropna()
    se['logUP'] = np.log10(np.maximum(0.001, se.UP_percentage))
    return se


def cumulative(x, steps=100):
    x = np.array(x)
    X = np.linspace(np.min(x), np.max(x), steps)
    return X, [np.mean(x <= z) for z in X]


def quintile_curves(se, plot_by='logUP'):
    """Five cumulative curves, genes split by expression quintile."""
    qs = [se.expression.quantile(q) for q in (0.2, 0.4, 0.6, 0.8)]
    bounds = [(-np.inf, qs[0]), (qs[0], qs[1]), (qs[1], qs[2]),
              (qs[2], qs[3]), (qs[3], np.inf)]
    curves = []
    for lo, hi in bounds:
        genes = se.loc[(se.expression > lo) & (se.expression <= hi), 'gene']
        sub = se.loc[se.ensembl_id.isin(genes), plot_by]
        curves.append(cumulative(sub))
    return curves


def make_selfsorted_grid_data(tissues=None, median_expression=None):
    """Per tissue: quintile curves ranked by that tissue's OWN expression.

    This is the leftmost panel of UP_splicing_by_exppression.BA24_v_selected
    ('BA24 itself'), computed for every tissue. Slow: one junction table per
    tissue.
    """
    from scipy.stats import pearsonr
    from tqdm import tqdm

    if tissues is None:
        tissues = all_tissues()
    if median_expression is None:
        median_expression = load_median_expression()

    out = []
    for tissue in tqdm(tissues, position=0, leave=True):
        up_table = get_UP_table(tissue)
        se = get_splicing_expression(up_table, median_expression, tissue)
        rho, pval = pearsonr(np.log1p(se.expression), se.logUP)
        out.append({
            'tissue': tissue,
            'color': '#' + gtex_colors[tissue]['tissue_color_hex'],
            'curves': quintile_curves(se),
            'rho': float(rho), 'pvalue': float(pval), 'n': int(len(se)),
        })
    return out


# --------------------------------------------------------------------------- #

PLOT_READY_VARS = ['fig5A_data', 'fig5B_data', 'up_by_expression_data']


def run_all(data_dir='figure_data'):
    os.makedirs(data_dir, exist_ok=True)
    upf_df = make_upf_df()
    data = {
        'fig5A_data': make_scatter_grid_data(upf_df, 'logUPF3A'),
        'fig5B_data': make_scatter_grid_data(upf_df, 'RIN'),
        'up_by_expression_data': make_selfsorted_grid_data(),
    }
    for k, v in data.items():
        with open(os.path.join(data_dir, f'{k}.pickle'), 'wb') as fh:
            pickle.dump(v, fh)
    return data


def load_plot_data(data_dir='figure_data', vars=PLOT_READY_VARS):
    data = {}
    for k in vars:
        with open(os.path.join(data_dir, f'{k}.pickle'), 'rb') as fh:
            data[k] = pickle.load(fh)
    return data
