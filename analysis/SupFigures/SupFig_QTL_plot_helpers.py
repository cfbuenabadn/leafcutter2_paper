"""
Plotting for the QTL_analysis.ipynb supplementary panels.

Created with Figure code cleaner.
Source notebook: ../QTL_analysis.ipynb
"""

import numpy as np
from matplotlib import pyplot as plt


def plot_fig4_supfig1(data, figsize=(10, 3.5), dpi=300):
    """fig4_supfig1 -- genomic-control lambda in eQTL effect, per tissue.

    u-sQTLs (red) against p-sQTLs (blue), 49 tissues ordered by descending
    u-sQTL lambda, with the GTEx tissue colour bar beneath and a dashed
    lambda = 1 reference.
    """
    sorted_tissues = data['sorted_tissues']
    lam_u, lam_pp, colors = data['lambda_u'], data['lambda_pp'], data['colors']

    fig, axes = plt.subplots(figsize=figsize, nrows=2,
                             gridspec_kw={'height_ratios': [15, 1], 'hspace': 0.01},
                             dpi=dpi)

    for i in range(len(lam_u)):
        axes[0].plot([i, i], [0, lam_u[i]], c='tab:gray', linestyle='--', zorder=0)

    axes[0].scatter(np.arange(49), lam_u, c='tab:red', label='u-sQTLs')
    axes[0].scatter(np.arange(49), lam_pp, c='tab:blue', label='p-sQTLs')

    axes[1].bar(range(49), [1] * 49, width=1, color=colors)

    axes[0].plot([0, 49], [1, 1], c='tab:gray', alpha=0.9, linestyle='--')

    axes[0].set_xticks([])
    axes[1].set_xticks([])
    axes[1].spines[['top', 'bottom', 'right', 'left']].set_visible(False)
    axes[1].set_yticks([])

    axes[1].set_xlim([-1, 49.5])
    axes[0].set_xlim([-1, 49.5])
    axes[0].set_ylabel(r'$\lambda$ inflation' + '\nin eQTL effect')
    axes[1].set_xlabel('Tissues (ranked by significance)')
    axes[0].legend(frameon=False)
    axes[0].set_ylim([0, 21])
    axes[0].text(45, 1.5, r'$\lambda$ = 1')

    axes[1].set_xticks(np.arange(49), sorted_tissues, rotation=45, ha='right',
                       rotation_mode='anchor', size=7)
    return axes


def plot_sup_fig_lambda(sqtl_qq, prepare_qq, figsize=(12, 12), dpi=300,
                        ncols=7, nrows=7):
    """sup_fig_lambda -- QQ plot of eQTL p-values per tissue, 7 x 7.

    Grey: the null from the eQTL nominal pass (see the bug note in
    SupFig_QTL_helpers.get_var_eqtls -- in the published figure this is Testis
    in every panel). Blue: p-sQTLs. Red: u-sQTLs. Red dashed line is y = x.

    Unlike the other 7 x 7 grids here, the source does NOT blank the unused
    49th cell, so it is left as drawn.
    """
    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi)
    i = j = 0

    for tissue in sqtl_qq.keys():
        eqq = prepare_qq(sqtl_qq[tissue]['eQTLs'])
        pqq = prepare_qq(sqtl_qq[tissue]['ppsQTLs'])
        uqq = prepare_qq(sqtl_qq[tissue]['usQTLs'])

        ax[i, j].plot([0, 4], [0, 4], 'r--')
        for qq, color in ((eqq, 'tab:gray'), (pqq, 'tab:blue'), (uqq, 'tab:red')):
            ax[i, j].scatter(qq[0], qq[1], edgecolor=color, s=10, alpha=0.75,
                             facecolor='none', rasterized=True)

        if j == 0:
            ax[i, j].set_ylabel('-log10(p-value)', size=5)
        if i == nrows - 1:
            ax[i, j].set_xlabel('Quantiles', size=5)

        ax[i, j].tick_params(labelsize=4, length=0.5, pad=1)
        ax[i, j].set_title(tissue, size=5, pad=2)

        j += 1
        if j >= ncols:
            j = 0
            i += 1

    return ax
