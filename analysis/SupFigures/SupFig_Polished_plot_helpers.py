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


def _panel_axis(ax, i, j, xlabel, ylabel):
    if j == 0:
        ax.set_ylabel(ylabel, size=5)
    if i == NROW - 1:
        ax.set_xlabel(xlabel, size=5)
    ax.tick_params(labelsize=4, length=0.5, pad=1)


def plot_scatter_grid(data, xlabel, ylabel='% NMD juncs',
                      highlight='Artery-Tibial'):
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
        a.text(0.03, 0.03,
               r'$\rho = $' + f"{d['rho']:.2f}" + '\n' + fmt_p(d['pvalue'])
               + '\n' + f"n = {d['n']}",
               transform=a.transAxes, ha='left', va='bottom', fontsize=4.5)
        a.set_title(d['tissue'], size=5, pad=2)
        _panel_axis(a, i, j, xlabel, ylabel)
    fig.tight_layout()
    return fig, ax


def plot_selfsorted_grid(data, xlabel='Unproductive splicing',
                         ylabel='Cumulative distribution'):
    """Per tissue, cumulative %UP by that tissue's own expression quintile.

    Five curves per panel, shaded light-to-dark in the tissue's GTEx colour --
    the same palette the single-tissue version uses -- so the lowest-expression
    quintile is palest and the highest darkest.
    """
    fig, ax = _grid()
    for k, d in enumerate(data):
        i, j = divmod(k, NCOL)
        a = ax[i, j]
        palette = sns.light_palette(d['color'], n_colors=5)
        for (x, y), c in zip(d['curves'], palette):
            a.plot(x, y, c=c, linewidth=1.5)
        a.set_xlim([-2.1, 2.1])
        a.set_xticks([-1, 0, 1, 2], ['0.1%', '1%', '10%', '100%'], size=4)
        a.text(0.03, 0.97,
               r'Pearson $\rho = $' + f"{d['rho']:.2f}" + '\n' + fmt_p(d['pvalue'])
               + '\n' + f"n = {d['n']}",
               transform=a.transAxes, ha='left', va='top', fontsize=4.5)
        a.set_title(d['tissue'], size=5, pad=2)
        _panel_axis(a, i, j, xlabel, ylabel)
    fig.tight_layout()
    return fig, ax


def save_panel(name, plots_dir, dpi=300):
    for ext in ('png', 'pdf', 'svg'):
        plt.savefig(f'{plots_dir}/{name}.{ext}', dpi=dpi, bbox_inches='tight')
