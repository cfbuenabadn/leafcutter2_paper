"""
Plotting functions for the Fig2.ipynb supplementary panels.

Created with Figure code cleaner.

Every function takes plot-ready data (as produced by SupFigures_helpers) and
draws one panel. No data loading, no computation.
"""

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


# =========================================================================== #
# Fig2.ipynb group
# =========================================================================== #

def plot_median_scatter(median_data, fit, xlabel, annotation=None, annotation_xy=None,
                        figsize=(3, 3), dpi=300):
    """sup_fig4A / sup_fig4B / sup_fig4UPF1 / sup_fig4UPF2 / sup_fig4UPF3B.

    Per-tissue median of `xvar` against the median unproductive-read
    percentage, one point per GTEx tissue in its GTEx colour, with the robust
    linear fit and its subsample band.

    All five panels are this same plot; only the x variable, the axis label and
    the hand-written rho/p annotation differ between them in Fig2.ipynb.
    """
    data, xvar, colors = median_data['data'], median_data['xvar'], median_data['colors']

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.scatter(data[xvar], data.pct, c=colors)

    mean_x = fit['mean_x']
    x = np.linspace(fit['x_min'] - mean_x, fit['x_max'] - mean_x, 100)

    ax.plot(x + mean_x, (fit['slope'] * (x + mean_x)) + fit['const'],
            color='black', linestyle='--')
    ax.fill_between(x + mean_x,
                    (fit['slope_90'] * np.array(x)) + fit['const'] + (fit['slope'] * mean_x),
                    (fit['slope_10'] * x) + fit['const'] + (fit['slope'] * mean_x),
                    color='tab:gray', alpha=0.2, label='Confidence Interval', linewidth=0)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('Median percentager\nunproductive splicing reads')
    ax.spines[['top', 'right']].set_visible(False)

    if annotation is not None and annotation_xy is not None:
        ax.text(annotation_xy[0], annotation_xy[1], annotation)

    return ax


def plot_per_tissue_grid(series, xlabel, ylabel='% NMD juncs', highlight='Artery-Tibial',
                         figsize=(9, 9), dpi=300, ncols=7, nrows=7):
    """sup_fig5A / sup_fig5B -- one small panel per GTEx tissue, 7 x 7.

    Same plot as plot_median_scatter but per tissue and per sample, with the
    Spearman rho and p annotated in each cell. `highlight` is drawn with a black
    fit line instead of red (Artery-Tibial in Fig2.ipynb); the unused 49th cell
    is blanked.
    """
    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi)
    i = j = 0

    for s in series:
        a = ax[i, j]
        a.scatter(s['x'], s['y'], edgecolor=s['color'], s=10, alpha=0.75,
                  facecolor='none', rasterized=True)

        if j == 0:
            a.set_ylabel(ylabel, size=5)
        if i == nrows - 1:
            a.set_xlabel(xlabel, size=5)
        a.tick_params(labelsize=4, length=0.5, pad=1)

        mean_x = s['mean_x']
        x = np.linspace(s['x_min'] - mean_x, s['x_max'] - mean_x, 100)
        line_color = 'black' if s['tissue'] == highlight else 'red'
        a.plot(x + mean_x, (s['slope'] * (x + mean_x)) + s['const'],
               color=line_color, linestyle='--')
        a.fill_between(x + mean_x,
                       (s['slope_90'] * np.array(x)) + s['const'] + (s['slope'] * mean_x),
                       (s['slope_10'] * x) + s['const'] + (s['slope'] * mean_x),
                       color='tab:gray', alpha=0.2, label='Confidence Interval', linewidth=0)

        a.text(np.quantile(a.get_xlim(), 0.05), np.quantile(a.get_ylim(), 0.05),
               r'$\rho=$' + f"{s['spearman_rho']:.2f}" + '\np-val ' + r'$=$'
               + f"{s['spearman_pval']:.2e}", fontsize=5)
        a.set_title(s['tissue'], size=5, pad=2)

        j += 1
        if j >= ncols:
            j = 0
            i += 1

    # 49 tissues in a 7x7 grid leaves the last cell empty
    ax[-1, -1].set_xticks([])
    ax[-1, -1].set_yticks([])
    ax[-1, -1].spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    ax[-1, -1].set_facecolor('none')

    return ax
