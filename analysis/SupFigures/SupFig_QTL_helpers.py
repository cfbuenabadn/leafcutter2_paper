"""
Data loading for the QTL_analysis.ipynb supplementary panels.

Created with Figure code cleaner.
Source notebook: ../QTL_analysis.ipynb

Serves fig4_supfig1 (lambda inflation per tissue) and sup_fig_lambda (per-tissue
QQ plots). Both share their data chain with main Figure 4, so this module reuses
../Figure4/Figure4_helpers.py rather than duplicating the permutation-pass and
eQTL nominal-pass readers.

Heavy packages are imported inside the functions that use them.
"""

import os
import gzip
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Figure4'))

from Figure4_helpers import (          # noqa: E402
    load_sqtl_stats,
    get_perm_counts,
    get_sqtl_expression_effects,
    tissue_color,
    gtex_colors,
    SQTL_DIR,
    EQTL_DIR,
)

PLOTS_DIR = 'plots'
REVISION_PLOTS_DIR = 'revision_plots'


# =========================================================================== #
# fig4_supfig1 -- lambda inflation per tissue (cheap: reads one .tsv.gz)
# =========================================================================== #

def make_fig4_supfig1_data():
    """Plot-ready arrays for fig4_supfig1, ordered by descending u-sQTL lambda."""
    sqtl_stats = load_sqtl_stats()

    u = sqtl_stats.loc[sqtl_stats.sqtl_type == 'u_sqtl'].set_index('tissue')
    sorted_lambda_tissues = list(u.sort_values('lambda_inflation').index[::-1])

    return {
        'sorted_tissues': sorted_lambda_tissues,
        'colors': [tissue_color(t) for t in sorted_lambda_tissues],
        'lambda_u': np.array(u.loc[sorted_lambda_tissues].lambda_inflation),
        'lambda_pp': np.array(
            sqtl_stats.loc[sqtl_stats.sqtl_type == 'pp_sqtl']
            .set_index('tissue').loc[sorted_lambda_tissues].lambda_inflation),
    }


# =========================================================================== #
# sup_fig_lambda -- per-tissue QQ plots (expensive)
# =========================================================================== #

def get_var_eqtls(tissue, max_pvals=10000, use_source_tissue_bug=True):
    """A null distribution of eQTL nominal p-values.

    !! REPRODUCES A BUG IN THE SOURCE NOTEBOOK !!

    QTL_analysis.ipynb defines this as `get_var_eqtls(tissue, ...)` but hard-codes
    Testis in the path:

        nom_file = f'../code/results/eqtl/GTEx/Testis/cis_100000/nom/{chrom}.txt.gz'

    so the `tissue` argument is ignored and every panel of sup_fig_lambda shows
    the SAME Testis null (the grey points), not that panel's own tissue. That is
    what the published figure contains, so it is the default here.

    Pass use_source_tissue_bug=False to read the requested tissue instead, which
    is what the code was evidently meant to do. That changes the figure.
    """
    read_tissue = 'Testis' if use_source_tissue_bug else tissue

    chroms = ['chr' + str(i) for i in range(1, 23)]
    epvals = []
    for chrom in chroms:
        nom_file = f'{EQTL_DIR}/{read_tissue}/cis_100000/nom/{chrom}.txt.gz'
        with gzip.open(nom_file, 'rb') as fh:
            fh.readline()
            for line in fh:
                epvals.append(float(line.decode().rstrip().split('\t')[11]))

    return np.random.choice(epvals, max_pvals, replace=False)


def make_sqtl_qq(tissues=None, use_source_tissue_bug=True, seed=None):
    """Per-tissue p-value vectors for the sup_fig_lambda QQ grid.

    For each tissue: nominal eQTL p-values of the variants behind its u-sQTLs
    and its p-sQTLs, plus a null drawn from the eQTL nominal pass (see the bug
    note in get_var_eqtls).

    Very slow -- one tabix query per significant sQTL, and the null reads whole
    eQTL nominal chromosomes. The originals left the sampling unseeded.
    """
    from tqdm import tqdm

    if tissues is None:
        tissues = sorted(os.listdir(SQTL_DIR))
    if seed is not None:
        np.random.seed(seed)

    sqtl_qq = {}
    for tissue in tqdm(tissues, position=0, leave=True):
        perm_counts = get_perm_counts(tissue)
        u_sqtl_effects = get_sqtl_expression_effects(
            perm_counts, tissue, ctype='PR,UP', itype='UP', qmax=1e-1, clu_psi_min=0.1)
        pp_sqtl_effects = get_sqtl_expression_effects(
            perm_counts, tissue, ctype='PR', itype='PR', qmax=1e-1, clu_psi_min=0.1)

        # The source sizes the null to the u-sQTL count of this tissue.
        max_sqtls = len(u_sqtl_effects[2])
        epvals = get_var_eqtls(tissue, max_pvals=max_sqtls,
                               use_source_tissue_bug=use_source_tissue_bug)

        sqtl_qq[tissue] = {'usQTLs': u_sqtl_effects[2],
                           'ppsQTLs': pp_sqtl_effects[2],
                           'eQTLs': epvals}
    return sqtl_qq


def prepare_qq(y_):
    """Expected/observed -log10 p pairs for one QQ series (QTL_analysis.ipynb)."""
    y = sorted(-np.log10(np.array(y_)))
    x = sorted(-np.log10(np.linspace(1, len(y) + 1, len(y)) / len(y)))
    return x, y
