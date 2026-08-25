"""
Plotting functions for the Python panel of Figure 5.

Created with Figure code cleaner.
Source notebook: ../Figure5.ipynb

Every function takes fully plot-ready data (as produced by
`Figure5_Python_helpers.make_Fig5D_right_data`) and draws one panel. No data loading,
no computation.
"""

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def fmt_p(p):
    """P-value formatting used by Fig5D_annotated."""
    if p >= 1e-3:
        return f'{p:.3g}'
    m, e = f'{p:.2e}'.split('e')
    return f'{m}$\\times$10$^{{{int(e)}}}$'


def fmt_annotation(s):
    """'beta = x.xx / P = ...' block, from the nominal stats of the plotted donors."""
    return f'$\\beta$ = {s["beta"]:.2f}\nP = {fmt_p(s["pval"])}'


def plot_Fig5D_right(Fig5D_right_data, panels=(('eQTL', 'tab:blue'), ('u-sQTL', 'tab:red')),
                         headroom=0.20):
    """Fig5D_annotated -- TSPAN14 eQTL and u-sQTL boxplots by genotype.

    Per-genotype n goes under each x tick label; the nominal effect size and P
    -- computed from the plotted donors -- are centred at the top of each panel, with `headroom` of extra y range to
    make space. The two panels share a y range and only the left one keeps its
    axis, so they read as one plot.
    """
    counts = Fig5D_right_data['counts']
    geno_order = Fig5D_right_data['geno_order']
    stats_by_panel = Fig5D_right_data['stats_by_panel']

    fig, axes = plt.subplots(ncols=2, figsize=(4, 4), width_ratios=[1, 1], dpi=500)

    for ax, (variable, color) in zip(axes, panels):
        sub = Fig5D_right_data['panels'][variable]
        sns.boxplot(data=sub, x='geno', y='value', ax=ax, order=geno_order,
                    boxprops=dict(facecolor='white', edgecolor=color), showfliers=False,
                    whiskerprops=dict(color=color),
                    capprops=dict(color=color),
                    medianprops=dict(linewidth=2, color='black'), linewidth=2)
        sns.stripplot(data=sub, x='geno', y='value', alpha=0.9, order=geno_order,
                      edgecolor=color, facecolor='none', linewidth=1, ax=ax,
                      rasterized=True)

    axes[0].spines[['top', 'right']].set_visible(False)
    axes[1].spines[['left', 'top', 'right']].set_visible(False)
    axes[1].set_yticks([])
    axes[1].set_ylabel('')

    ax0y = axes[0].get_ylim()
    ax1y = axes[1].get_ylim()

    ybottom = np.min([ax0y[0], ax1y[0]])
    # np.max, not np.min: with np.min the top eQTL point (2.04) fell outside the axes
    ytop = np.max([ax0y[1], ax1y[1]])
    ylim = (ybottom, ytop + headroom * (ytop - ybottom))   # headroom for the annotation

    for ax, (variable, _) in zip(axes, panels):
        ax.set_ylim(ylim)
        ax.set_title(variable, pad=14)
        ax.set_xlabel('')
        ax.set_xticks(range(len(geno_order)))
        ax.set_xticklabels([f'{g}\nn={counts[g]}' for g in geno_order])
        ax.text(0.5, 0.99, fmt_annotation(stats_by_panel[variable]),
                transform=ax.transAxes, ha='center', va='top', fontsize=8)

    axes[0].set_ylabel('Normalized phenotype')

    return axes
