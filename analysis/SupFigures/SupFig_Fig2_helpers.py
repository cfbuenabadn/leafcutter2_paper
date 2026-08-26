"""
Data loading and processing for the Fig2.ipynb supplementary panels.

Created with Figure code cleaner.
Source notebook: ../Fig2.ipynb

Serves sup_fig4A, sup_fig4B, sup_fig4UPF1/2/3B, sup_fig5A, sup_fig5B and
sup_fig7A/B/C. The three GSEA panels from the same notebook are a separate
chain (gseapy + MSigDB) and live in SupFig_GSEA.ipynb.

Panels are written to plots/, mirroring ../../code/plots/.

Only standard packages are imported at module level; heavy or optional ones
(tqdm, scipy, statsmodels, gseapy, tabix, rpy2) are imported inside the
functions that use them.
"""

import os
import gzip
import pickle
import sys

import numpy as np
import pandas as pd

BASE = '/project/yangili1/cfbuenabadn/leafcutter2_paper'

# Two destinations, matching the source notebooks.
PLOTS_DIR = 'plots'
REVISION_PLOTS_DIR = 'revision_plots'

# --------------------------------------------------------------------------- #
# Reuse of the Figure 2 helpers
#
# The Fig2.ipynb supplementary panels share their whole data chain with the
# main Figure 2 panels, so rather than duplicating the GTEx TPM reader and the
# unproductive-percentage loader, this module imports them from ../Figure2.
# Nothing is copied; if those loaders change, these panels follow.
# --------------------------------------------------------------------------- #

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Figure2'))

from Figure2_helpers import (            # noqa: E402
    load_unproductive_pct,
    load_gene_tpm,
    build_upf_df,
    make_boxplot_df,
    gtex_colors,
    TEN_TISSUES,
    TEN_TISSUES_CLEAN,
    BOXPLOT_PALETTE,
)

# GTEx sample attributes, for the RIN column (Fig2.ipynb cell 39).
GTEX_SAMPLE_ATTRIBUTES = f'{BASE}/code/resources/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'

# Genes the supplementary panels need on top of the two Figure 2 uses.
# Fig2.ipynb log-transforms all of these; sup_fig7A/B/C plot MYOM2/SRSF3/DLG4.
SUP_GENES = ['UPF1', 'UPF2', 'UPF3A', 'UPF3B', 'GABBR1',
             'ABCA5', 'AKAP8L', 'CNNM3', 'DLG4', 'EIF4A2',
             'MRTO4', 'MYOM2', 'NOC2L', 'SFT2D1', 'SRSF3']

LOG_GENES = ['UPF1', 'UPF2', 'UPF3A', 'UPF3B', 'GABBR1',
             'ABCA5', 'AKAP8L', 'CNNM3', 'DLG4', 'EIF4A2',
             'MRTO4', 'MYOM2', 'NOC2L', 'SFT2D1', 'SRSF3']


# =========================================================================== #
# Fig2.ipynb group -- sup_fig4*, sup_fig5*, sup_fig7*
# =========================================================================== #

def load_sample_rin(attributes_table=GTEX_SAMPLE_ATTRIBUTES):
    """Per-sample RIN from the GTEx sample attributes file (Fig2.ipynb cell 39)."""
    metadata = pd.read_csv(attributes_table, sep='\t')
    metadata_rin = metadata[['SAMPID', 'SMRIN']]
    metadata_rin.columns = ['entity:sample_id', 'RIN']
    return metadata_rin


def build_sup_upf_df(genes=SUP_GENES, log_genes=LOG_GENES):
    """`upf_df` as Fig2.ipynb has it by the time the supplementary panels run.

    Per-sample gene TPM joined to the unproductive-read percentage, with the
    log10 gene columns and the RIN column added (Fig2.ipynb cells 24, 25, 40).
    """
    pct_df = load_unproductive_pct()
    gene_df = load_gene_tpm(genes=genes)
    upf_df = build_upf_df(gene_df, pct_df)

    for gene in log_genes:
        upf_df[f'log{gene}'] = np.log10(upf_df[gene])
    upf_df['ratio'] = upf_df.UPF3A / upf_df.UPF3B

    upf_df = upf_df.merge(load_sample_rin(), on='entity:sample_id')
    return upf_df


def make_median_by_tissue(upf_df, xvar):
    """Per-tissue medians of `xvar` and pct, with the GTEx colour of each tissue.

    This is the `upf_df.groupby('tissue')[[xvar, 'pct']].median()` that
    sup_fig4A / sup_fig4B / sup_fig4UPF* all start from.
    """
    data = upf_df.groupby('tissue')[[xvar, 'pct']].median()
    colors = ['#' + gtex_colors[t]['tissue_color_hex'] for t in data.index]
    return {'data': data, 'xvar': xvar, 'colors': colors}


def fit_rlm_with_band(data, xvar, n_draws=100, draw_size=30, seed=None):
    """Robust linear fit of pct on `xvar`, with a subsample band.

    Reproduces the block shared by sup_fig4A/4B/4UPF*: statsmodels RLM on all
    49 tissue medians for the line, and the 10th-90th percentile of the RLM
    slope over `n_draws` random subsamples of `draw_size` tissues for the band.

    `draw_size` is 30 for the gene panels and 15 for sup_fig4A (RIN) -- the
    source notebook uses different values and that changes the band width.

    The originals left the RNG unseeded; pass `seed` to make the band
    reproducible.
    """
    import statsmodels.api as sm

    Y = data.pct
    X = sm.add_constant(data[xvar])

    rng = np.random.default_rng(seed) if seed is not None else np.random
    rho_list = []
    for _ in range(n_draws):
        idx = rng.choice(np.arange(len(data)), draw_size, replace=False)
        res_ = sm.RLM(Y.iloc[idx], sm.add_constant(X.iloc[idx])).fit()
        rho_list.append(res_.params.iloc[1])

    results = sm.RLM(Y, X).fit()
    return {
        'slope': results.params.loc[xvar],
        'const': results.params.loc['const'],
        'slope_90': np.quantile(rho_list, 0.9),
        'slope_10': np.quantile(rho_list, 0.1),
        'mean_x': float(np.mean(X[xvar])),
        'x_min': float(np.min(data[xvar])),
        'x_max': float(np.max(data[xvar])),
    }


def make_per_tissue_series(upf_df, xvar, n_draws=100, draw_size=50, seed=None):
    """Per-tissue scatter + RLM fit for the 7x7 grids (sup_fig5A / sup_fig5B).

    One entry per tissue, in `sort_values('tissue').groupby('tissue')` order,
    each carrying the points, the tissue colour and the same RLM-with-band fit
    the single-panel version uses.
    """
    import statsmodels.api as sm
    from scipy.stats import spearmanr

    rng = np.random.default_rng(seed) if seed is not None else np.random
    series = []
    for tissue, df in upf_df.sort_values('tissue').groupby('tissue'):
        Y = df.pct
        X = sm.add_constant(df[xvar])

        rho_list = []
        for _ in range(n_draws):
            idx = rng.choice(np.arange(df.shape[0]), draw_size, replace=False)
            res_ = sm.RLM(Y.iloc[idx], sm.add_constant(X.iloc[idx])).fit()
            rho_list.append(res_.params.iloc[1])

        results = sm.RLM(Y, X).fit()
        spear = spearmanr(df[xvar], df.pct)   # annotated on each subplot
        series.append({
            'tissue': tissue,
            'spearman_rho': spear[0],
            'spearman_pval': spear[1],
            'color': '#' + gtex_colors[tissue]['tissue_color_hex'],
            'x': np.array(df[xvar]),
            'y': np.array(df.pct),
            'slope': results.params.loc[xvar],
            'const': results.params.loc['const'],
            'pvalue': results.pvalues.loc[xvar],
            'slope_90': np.quantile(rho_list, 0.9),
            'slope_10': np.quantile(rho_list, 0.1),
            'mean_x': float(np.mean(df[xvar])),
            'x_min': float(np.min(df[xvar])),
            'x_max': float(np.max(df[xvar])),
        })
    return series
