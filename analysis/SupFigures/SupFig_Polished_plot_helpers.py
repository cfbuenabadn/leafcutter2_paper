"""Plotting for the three polished supplementary grids.

Every grid is 7 x 7 over 49 GTEx tissues. 7 x 7 = 49 exactly, so there is NO
spare panel -- the earlier versions of sup_fig5A and sup_fig5B ended with

    ax[-1,-1].set_xticks([]); ax[-1,-1].set_yticks([])
    ax[-1,-1].spines[[...]].set_visible(False)

which blanks the axes of the 49th tissue's own panel rather than tidying an
unused one. That block is not reproduced here.
"""
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

NROW = NCOL = 7


def fmt_p(p):
    """P value in full, never as an asterisk and never rounded to zero."""
    if p == 0:
        return 'P < 1e-308'
    if p < 1e-3:
        m, e = f'{p:.1e}'.split('e')
        return rf'P = {m}$\times 10^{{{int(e)}}}$'
    return f'P = {p:.3g}'


def _grid(figsize=(9, 9), dpi=300):
    fig, ax = plt.subplots(NROW, NCOL, figsize=figsize, dpi=dpi)
    return fig, ax


def _panel_axis(ax, i, j, xlabel, ylabel, size=5, tick=4):
    if j == 0:
        ax.set_ylabel(ylabel, size=size)
    if i == NROW - 1:
        ax.set_xlabel(xlabel, size=size)
    ax.tick_params(labelsize=tick, length=0.5, pad=1)


def stat_text(d, metric):
    """rho / P / n block for one panel, for the chosen correlation metric."""
    st = d[metric]
    label = 'Pearson' if metric == 'pearson' else 'Spearman'
    return (label + r' $\rho = $' + f"{st['rho']:.2f}" + '\n'
            + fmt_p(st['pvalue']) + '\n' + f"n = {d['n']}")


def plot_scatter_grid(data, xlabel, ylabel='% NMD juncs',
                      highlight='Artery-Tibial', metric='spearman'):
    """sup_fig5A / sup_fig5B -- one scatter + RLM line per tissue.

    Every panel keeps its axes, including the 49th.
    """
    fig, ax = _grid()
    for k, d in enumerate(data):
        i, j = divmod(k, NCOL)
        a = ax[i, j]
        a.scatter(d['x'], d['y'], edgecolor=d['color'], s=10, alpha=0.75,
                  facecolor='none', rasterized=True)
        a.plot(d['fit_x'], d['fit_y'], linestyle='--',
               color='black' if d['tissue'] == highlight else 'red')
        a.text(0.03, 0.03, stat_text(d, metric),
               transform=a.transAxes, ha='left', va='bottom', fontsize=4.5)
        a.set_title(d['tissue'], size=5, pad=2)
        _panel_axis(a, i, j, xlabel, ylabel)
    fig.tight_layout()
    return fig, ax


def plot_selfsorted_grid(data, xlabel='Unproductive splicing',
                         ylabel='Cumulative distribution', metric='pearson'):
    """Per tissue, cumulative %UP by that tissue's own expression quintile.

    Five curves per panel, shaded light-to-dark in the tissue's GTEx colour --
    the same palette the single-tissue version uses -- so the lowest-expression
    quintile is palest and the highest darkest.
    """
    # Layout follows UP_splicing_by_exppression.All: a large grid with ticks
    # only on the outer edge, so the 49 panels stay readable.
    fig, ax = _grid(figsize=(21, 21))
    for k, d in enumerate(data):
        i, j = divmod(k, NCOL)
        a = ax[i, j]
        palette = sns.light_palette(d['color'], n_colors=5)
        for (x, y), c in zip(d['curves'], palette):
            a.plot(x, y, c=c, linewidth=3)
        a.set_xlim([-2.1, 2.1])
        if i == NROW - 1:
            a.set_xticks([-1, 0, 1, 2], ['0.1%', '1%', '10%', '100%'])
        else:
            a.set_xticks([])
        if j != 0:
            a.set_yticks([])
        a.text(0.03, 0.97, stat_text(d, metric),
               transform=a.transAxes, ha='left', va='top', fontsize=9)
        a.set_title(d['tissue'], size=12, pad=3)
        _panel_axis(a, i, j, xlabel, ylabel, size=12, tick=9)
    fig.tight_layout()
    return fig, ax


def plot_corr_matrix(data, metric='pearson'):
    """corr_across_tissues -- 49 x 49 rho, rows and columns in Ward order.

    Row: the tissue whose unproductive splicing is measured. Column: the tissue
    whose expression ranks the genes. Tissue colour bars on both margins.
    """
    import pandas as pd
    m = data[metric]
    order = m['order']
    df = pd.DataFrame(m['matrix'], index=data['tissues'], columns=data['tissues'])
    df = df.iloc[order, order]
    colors = [data['colors'][k] for k in order]

    g = sns.clustermap(df, col_cluster=False, row_cluster=False,
                       col_colors=colors, row_colors=colors,
                       cmap='inferno_r', figsize=(20, 20))
    g.ax_heatmap.set_xticks([])
    g.ax_heatmap.tick_params(labelsize=16, which='both')
    g.ax_heatmap.set_xlabel('Gene expression order', fontsize=20)
    g.ax_heatmap.yaxis.set_label_position('left')
    g.ax_heatmap.set_ylabel('UP splicing', labelpad=50, fontsize=20)
    g.cax.yaxis.set_label_position('left')
    label = 'Pearson' if metric == 'pearson' else 'Spearman'
    g.cax.set_ylabel(label + ' ' + r'$\rho$', labelpad=15, fontsize=20)
    g.cax.tick_params(labelsize=16)
    return g


def save_panel(name, plots_dir, dpi=300):
    for ext in ('png', 'pdf', 'svg'):
        plt.savefig(f'{plots_dir}/{name}.{ext}', dpi=dpi, bbox_inches='tight')
