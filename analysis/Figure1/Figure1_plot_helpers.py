"""Plotting for the salvaged Figure 1 panels."""
import numpy as np
from matplotlib import pyplot as plt

from Figure1_helpers import (QUARTILES, CATEGORIES, CATEGORY_LABELS,
                             CATEGORY_COLORS)


def fmt_p(p):
    if not np.isfinite(p):
        return 'n/a'
    if p == 0:
        return 'P < 1e-308'
    if p < 1e-3:
        m, e = f'{p:.1e}'.split('e')
        return rf'P = {m}$\times10^{{{int(e)}}}$'
    return f'P = {p:.3g}'


def plot_fig1d(d, quartiles=None, ax=None):
    """Stacked bars: LeafCutter2 class composition by usage quartile."""
    quartiles = quartiles or QUARTILES
    if ax is None:
        _, ax = plt.subplots(figsize=(3.2, 3), dpi=300)
    bottom = np.zeros(len(quartiles))
    for cat in CATEGORIES:
        sub = d.set_index(['quartile', 'category']).loc[
            [(q, cat) for q in quartiles], 'pct_junctions'].to_numpy()
        ax.bar(range(len(quartiles)), sub, bottom=bottom, width=0.75,
               color=CATEGORY_COLORS[cat], label=CATEGORY_LABELS[cat],
               edgecolor='white', linewidth=0.5)
        bottom += sub
    n = d.set_index(['quartile', 'category']).loc[
        [(q, CATEGORIES[0]) for q in quartiles], 'n_in_quartile'].to_numpy()
    ax.set_xticks(range(len(quartiles)),
                  [f'{q}\nn = {v:,}' for q, v in zip(quartiles, n)], size=7)
    ax.set_ylim(0, 100)
    ax.set_ylabel('% of junctions')
    ax.set_xlabel('Junction usage quartile')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc='upper center',
              bbox_to_anchor=(0.5, 1.22), ncol=3)
    return ax


def plot_fig1g(series, comparisons=None, ncols=2, categories=('productive', 'unproductive')):
    """ECDF of log2 fold change per class, one axes per perturbation."""
    labels = comparisons or sorted({s['comparison_label'] for s in series},
                                   key=lambda l: [s['comparison_label'] for s in series].index(l))
    nrows = int(np.ceil(len(labels) / ncols))
    fig, ax = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.8 * nrows), dpi=300)
    ax = np.atleast_1d(ax).ravel()
    for k, label in enumerate(labels):
        a = ax[k]
        for s in series:
            if s['comparison_label'] != label or s['category'] not in categories:
                continue
            a.plot(s['ecdf_x'], s['ecdf_y'], color=s['color'], linewidth=1.8,
                   label=f"{s['category_label']} (n = {s['n']:,})")
        a.axvline(0, color='grey', linestyle=':', linewidth=0.8)
        a.set_title(label, size=8)
        a.set_xlabel(r'log$_2$ fold change')
        if k % ncols == 0:
            a.set_ylabel('Cumulative fraction')
        a.spines[['top', 'right']].set_visible(False)
        a.legend(frameon=False, fontsize=6, loc='upper left')
    for a in ax[len(labels):]:
        a.set_visible(False)
    fig.tight_layout()
    return fig, ax


def plot_fig1h(h, cmap='Blues'):
    """Heatmap of the NMD-efficiency rules: delta log2FD with its P value.

    Sequential ramp, cells labelled 'log2FD:' / 'P:', matching the published
    panel rather than plot_rules.py's stale diverging version.
    """
    rules = list(dict.fromkeys(h.rule))
    classes = ['Unproductive', 'Productive']
    M = np.full((len(rules), len(classes)), np.nan)
    for i, r in enumerate(rules):
        for j, c in enumerate(classes):
            row = h[(h.rule == r) & (h['class'] == c)]
            if len(row):
                M[i, j] = row.delta_log2fd.iloc[0]

    fig, ax = plt.subplots(figsize=(4.2, 2.6), dpi=300)
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=np.nanmax(M), aspect='auto')
    for i, r in enumerate(rules):
        for j, c in enumerate(classes):
            row = h[(h.rule == r) & (h['class'] == c)]
            if not len(row) or not np.isfinite(row.delta_log2fd.iloc[0]):
                ax.text(j, i, 'N/A', ha='center', va='center', size=7, color='0.35')
                continue
            v = row.delta_log2fd.iloc[0]
            shade = 'white' if v > 0.6 * np.nanmax(M) else 'black'
            ax.text(j, i, f'log2FD: {v:.2f}\n{fmt_p(row.p_value.iloc[0])}',
                    ha='center', va='center', size=6.5, color=shade)
    ax.set_xticks(range(len(classes)), classes, size=8)
    ax.set_yticks(range(len(rules)), [r.replace(' of the PTC', '') for r in rules], size=7)
    fig.colorbar(im, ax=ax, label=r'$\Delta$ log2FD, naRNA vs polyA', shrink=0.85)
    fig.tight_layout()
    return fig, ax


def save_panel(name, plots_dir='plots', dpi=300):
    import os
    os.makedirs(plots_dir, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        plt.savefig(f'{plots_dir}/{name}.{ext}', dpi=dpi, bbox_inches='tight')
