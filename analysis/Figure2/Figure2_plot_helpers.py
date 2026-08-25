"""
Plotting functions for Figure 2.

Created with Figure code cleaner.
Source notebook: ../Fig2.ipynb

Every function here takes fully plot-ready data (as produced by
`Figure2_helpers.run_all`) and draws one panel. No data loading, no
computation.

Panel letters follow the published Figure 2 caption: plot_fig2a (was fig2A),
plot_fig2b (the volcano, was fig2C), plot_gene_boxplots (Fig. 2e, was
fig2_boxplots).
"""

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def plot_fig2a(panels, tissue_names, savefig=False):
    """Fig. 2a -- per-tissue unproductive-read percentage and UPF3A TPM.

    Parameters
    ----------
    panels : list of dict
        One entry per tissue, already in left-to-right plotting order, with
        keys 'color', 'pct' (sorted) and 'upf3a' (ordered by pct).
    tissue_names : list of str
        Display names, same order as `panels`.
    """
    fig, ax = plt.subplots(nrows=3, figsize=(10, 4),
                           gridspec_kw={'height_ratios': [15, 10, 1.5], 'hspace': 0.005},
                           dpi=300)

    c_list = []
    for pos, panel in enumerate(panels):
        c = panel['color']
        c_list.append(c)

        y = panel['pct']
        x = np.linspace(pos - 0.4, pos + 0.4, len(y))
        ax[0].scatter(x, y, alpha=0.2, edgecolor='none', facecolor=c, s=5, rasterized=True)

        median_y = np.median(y)
        ax[0].plot([pos - 0.4, pos + 0.4], [median_y, median_y],
                   c='black', linewidth=0.75, alpha=0.95)

        y_tpm = panel['upf3a']
        ax[1].scatter(x, y_tpm, alpha=0.2, edgecolor='none', facecolor=c, s=5, rasterized=True)

    ax[0].set_xlim((-1, 49.25))
    ax[0].set_ylim([0, 2.6])
    ax[1].set_ylim([0, 240])
    ax[1].set_xlim((-1, 49.25))
    ax[2].set_xlim((-1, 49.25))

    ax[1].set_xticks([])
    ax[0].set_xticks([])
    ax[2].set_xticks(np.arange(49), tissue_names, rotation=45,
                     ha='right', rotation_mode='anchor', size=7)

    ax[2].bar(range(49), [1] * 49, width=1, color=c_list)

    ax[2].spines[['top', 'bottom', 'right', 'left']].set_visible(False)
    ax[2].set_yticks([])
    ax[0].spines[['top', 'right']].set_visible(False)
    ax[1].set_ylabel('UPF3A\n(TPM)', size=8)
    ax[0].set_ylabel('Percentage unproductive\njunction reads', size=8)

    ax[1].text(3, 210, 'Spearman correlation:', fontsize=8)
    ax[1].text(3, 180, r'$\rho=0.52$', fontsize=8)
    ax[1].text(3, 150, 'p-val < 1e-16', fontsize=8)

    if savefig:
        for ext in ['png', 'pdf', 'svg']:
            plt.savefig(f'{savefig}.{ext}', dpi=300, bbox_inches='tight')

    return fig, ax


# Fig. 2a draws no boxes: the black horizontal bar on each tissue is the median of
# that tissue's per-sample unproductive-read percentage, over the same points.
FIG2A_CENTRE_LEGEND = (
    "Each point is one sample; the black horizontal bar is the median of that "
    "tissue's distribution. No error bars are shown."
)


def plain_pvalue(p):
    """Plain-text p-value for legend prose (never rounds to zero)."""
    if p == 0:
        return 'P < 1 x 10-308'      # below double-precision underflow
    if p < 1e-3:
        mantissa, exponent = f'{p:.1e}'.split('e')
        return f'P = {mantissa} x 10{int(exponent)}'
    return f'P = {p:.3g}'


def boxplot_legend(stats, gene='GABBR1'):
    """Concise legend text for Fig. 2e, in published-legend register.

    Covers the checklist items in two added sentences: unit of study and n
    (biological, distinct donors, n below each box), the box-plot definition
    (centre, bounds, whiskers, percentile, minima/maxima via the overlaid
    points), and the test with its sidedness, statistic, exact P and
    multiple-comparison status. seaborn/matplotlib defaults in force in
    `plot_gene_boxplots`: whis=1.5, showfliers=False.
    """
    test = stats[stats['primary']]

    if stats['primary'] == 'wilcoxon':
        comparison = (
            "Brain and non-brain tissues were compared by a two-sided Wilcoxon "
            "signed-rank test on per-donor medians (n = {n} paired donors, each "
            "contributing one value per group; no technical replicates were "
            "used): W = {stat:.0f}, {p}."
        ).format(n=test['n_pairs'], stat=test['statistic'],
                 p=plain_pvalue(test['pvalue']))
    else:
        comparison = (
            "Brain (n = {n_a} donors) and non-brain (n = {n_b} donors) tissues "
            "were compared by a two-sided Mann-Whitney U test on per-donor "
            "medians (no technical replicates were used): U = {stat:.0f}, {p}."
        ).format(n_a=test['n_a'], n_b=test['n_b'], stat=test['statistic'],
                 p=plain_pvalue(test['pvalue']))

    return (
        "Boxplots of {gene} expression level across GTEx tissues. High expression "
        "of {gene} in brain tissue is associated with high exon inclusion. Each "
        "point is one sample from a distinct GTEx donor, with n given below each "
        "box; boxes depict interquartiles with the median (red line) and whiskers "
        "extending to the most extreme value no greater than 1.5x IQR from the "
        "hinge, with all points overlaid. {comparison} A single comparison was "
        "performed, so no adjustment for multiple comparisons was applied."
    ).format(gene=gene, comparison=comparison)


def format_pvalue(p):
    """Format a p-value for display, without ever rounding it to `p = 0`."""
    if p == 0:
        return r'$p < 10^{-308}$'   # below double-precision underflow
    if p < 1e-3:
        mantissa, exponent = f'{p:.1e}'.split('e')
        return rf'$p = {mantissa} \times 10^{{{int(exponent)}}}$'
    return rf'$p = {p:.3g}$'


def plot_gene_boxplots(boxplot_df, gene, tissue_names, palette, stats=None, savefig=False):
    """Fig. 2e -- per-tissue TPM strip + box plot for a single gene.

    Parameters
    ----------
    stats : dict, optional
        Output of `Figure2_helpers.make_boxplot_stats`. When given, the per-box
        sample size is appended to each tick label and the brain-vs-non-brain
        group test named by `stats['primary']` is annotated on the axes.
    """
    fig, ax = plt.subplots(figsize=(6, 3), dpi=300)

    sns.stripplot(data=boxplot_df, x='tissue', y=gene, hue='tissue',
                  palette=palette, alpha=0.2, rasterized=True, zorder=0, jitter=0.2)

    sns.boxplot(data=boxplot_df, x='tissue', y=gene,
                boxprops=dict(facecolor='none', edgecolor='black'),
                medianprops={'color': 'tab:red'},
                showfliers=False, width=.6, zorder=1)

    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel(rf'$\mathit{{{gene}}}$ (TPM)')

    labels = list(tissue_names)
    if stats is not None:
        labels = [f'{name}\nn = {n}' for name, n in zip(labels, stats['n_per_tissue'])]

    ax.set_xticks(np.arange(10), labels, rotation=45,
                  ha='right', rotation_mode='anchor', size=7)

    if stats is not None:
        # Boundary between the pooled brain tissues and the pooled non-brain ones.
        ax.axvline(len(stats['group_a']) - 0.5, color='black',
                   linestyle=':', linewidth=0.75, alpha=0.6, zorder=2)

        test = stats[stats['primary']]
        annotation = '{label}\n{test}: {stat_name} = {statistic:.0f}\n{p}'.format(
            label=test['label'], test=test['test'], stat_name=test['stat_name'],
            statistic=test['statistic'], p=format_pvalue(test['pvalue']))
        ax.text(0.98, 0.98, annotation, transform=ax.transAxes,
                ha='right', va='top', fontsize=7)

    if savefig:
        for ext in ['png', 'pdf', 'svg']:
            plt.savefig(f'{savefig}.{ext}', dpi=300, bbox_inches='tight')

    return fig, ax


def fig2b_legend(source_data, min_n=50, fdr=0.1, alpha=0.05):
    """Concise legend text for Fig. 2b, with every count filled in.

    Takes `fig2b_source_data` (from `Figure2_helpers.make_fig2b_source_data`).
    The first sentence is the existing published description; the rest adds the
    test, its sidedness, the statistic and degrees of freedom, why n varies, and
    the multiple-comparison adjustment.
    """
    up = source_data.loc[source_data.splicing_class == 'unproductive']
    pr = source_data.loc[source_data.splicing_class == 'productive']
    sig_col = f'significant_at_fdr_{fdr}'

    return (
        "Scatter plot showing the Spearman's correlation between delta PSI and "
        "differential gene expression Z-scores for tissue pairs (x-axis) and the "
        "significance of the correlation (y-axis). Correlations are two-sided "
        "Spearman's rho over the n genes with both a significant splicing change "
        "and a significant expression change in that pair, which varies because "
        "splicing is only observed where expression is sufficient in both tissues "
        "(df = n - 2; unproductive, left: n = {up_lo}-{up_hi}, median {up_med}; "
        "productive, right: n = {pr_lo}-{pr_hi}, median {pr_med}); pairs with "
        "n < {min_n} are not shown, leaving {n_pairs} tissue pairs. Exact n, rho, "
        "{conf}% confidence intervals and P values are in Source Data. The dashed "
        "line marks unadjusted P = 0.01 and coloured points have "
        "Benjamini-Hochberg-adjusted P <= {fdr} ({up_sig} of {n_pairs} pairs "
        "unproductive, {pr_sig} of {n_pairs} productive), adjusted separately "
        "within each panel."
    ).format(
        n_pairs=len(up), min_n=min_n, fdr=fdr, conf=int((1 - alpha) * 100),
        up_lo=int(up.n_genes.min()), up_hi=int(up.n_genes.max()),
        up_med=int(up.n_genes.median()),
        pr_lo=int(pr.n_genes.min()), pr_hi=int(pr.n_genes.max()),
        pr_med=int(pr.n_genes.median()),
        up_sig=int(up[sig_col].sum()), pr_sig=int(pr[sig_col].sum()),
    )


def plot_fig2b(series, savefig=False):
    """Fig. 2b -- volcano plots of splicing-vs-expression correlation across tissue pairs.

    `series` holds the four plot-ready (x, y) scatter sets:
    unproductive_ns / unproductive_sig / productive_ns / productive_sig.
    """
    fig, ax = plt.subplots(ncols=2, figsize=(6.2, 3),
                           gridspec_kw={'wspace': 0.1}, dpi=600)

    s = series['unproductive_ns']
    ax[0].scatter(s['x'], s['y'], edgecolor='tab:gray', facecolor='none',
                  s=10, alpha=0.5, rasterized=True, label='Not significant')

    s = series['unproductive_sig']
    ax[0].scatter(s['x'], s['y'], edgecolor='tab:red', facecolor='none',
                  s=10, alpha=0.5, rasterized=True, label='Unproductive splicing')

    ax[0].plot([-0.6, 0.2], [2, 2], linestyle='--', c='black')
    ax[0].text(-0.1, 3, r'$p=0.01$')
    ax[0].set_xlabel('Unproductive splicing v expression\n' + r'Spearman $\rho$')
    ax[0].set_ylabel(r'$-log10(p)$')
    ax[0].set_ylim([-0.5, 20])
    ax[0].spines[['top', 'right']].set_visible(False)

    s = series['productive_ns']
    ax[1].scatter(s['x'], s['y'], edgecolor='tab:gray', facecolor='none',
                  s=10, alpha=0.5, rasterized=True)

    # Off-axis proxy points, only there to build the legend (as in Fig2.ipynb).
    ax[1].scatter([0], [-50], edgecolor='tab:red', facecolor='none',
                  s=10, alpha=0.5, rasterized=True, label='Unproductive splicing')

    s = series['productive_sig']
    ax[1].scatter(s['x'], s['y'], edgecolor='tab:blue', facecolor='none',
                  s=10, alpha=0.5, rasterized=True, label='Productive splicing')

    ax[1].scatter([0], [-50], edgecolor='tab:gray', facecolor='none',
                  s=10, alpha=0.5, rasterized=True, label='Not significant')

    ax[1].plot([-0.2, 0.2], [2, 2], linestyle='--', c='black')
    ax[1].set_xlabel('Productive splicing v expression\n' + r'Spearman $\rho$')
    ax[1].set_ylim([-0.5, 20])
    ax[1].set_yticks([])
    ax[1].spines[['top', 'right', 'left']].set_visible(False)

    ax[1].legend(frameon=False)

    if savefig:
        for ext in ['png', 'pdf', 'svg']:
            plt.savefig(f'{savefig}.{ext}', dpi=600, bbox_inches='tight')

    return fig, ax
