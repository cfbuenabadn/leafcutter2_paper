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


def plot_fig1d(d, quartiles=None, ax=None, label_min_pct=3.0):
    """Stacked bars: LeafCutter2 class composition, ALL plus each usage quartile.

    Segments above `label_min_pct` carry their percentage and count, as in the
    published panel.
    """
    quartiles = quartiles or ['ALL'] + QUARTILES
    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 3.4), dpi=300)
    t = d.set_index(['quartile', 'category'])
    bottom = np.zeros(len(quartiles))
    for cat in CATEGORIES:
        vals = t.loc[[(q, cat) for q in quartiles], 'pct_junctions'].to_numpy()
        cnts = t.loc[[(q, cat) for q in quartiles], 'n_junctions'].to_numpy()
        ax.bar(range(len(quartiles)), vals, bottom=bottom, width=0.72,
               color=CATEGORY_COLORS[cat], label=CATEGORY_LABELS[cat],
               edgecolor='white', linewidth=0.5)
        for i, (v, c) in enumerate(zip(vals, cnts)):
            if v >= label_min_pct:
                ax.text(i, bottom[i] + v / 2, f'{v:.1f}%\n(n={c:,})', ha='center',
                        va='center', fontsize=6, color='white', fontweight='bold')
        bottom += vals
    n = t.loc[[(q, CATEGORIES[0]) for q in quartiles], 'n_in_quartile'].to_numpy()
    ax.set_xticks(range(len(quartiles)),
                  [f'{q}\nn = {v:,}' for q, v in zip(quartiles, n)], size=7)
    ax.set_ylim(0, 100)
    ax.set_ylabel('% of junctions')
    ax.set_xlabel('Junction usage quartile')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc='upper center',
              bbox_to_anchor=(0.5, 1.18), ncol=3)
    return ax


def plot_fig1e(d, quartiles=None, label_min_pct=3.0, palette=None):
    """Stacked bars: GENCODE composition of each class, one panel per class.

    Mirrors the published layout -- x is the usage quartile, the stack is the
    GENCODE v46 transcript type, and one subplot per LeafCutter2 class.
    """
    quartiles = quartiles or ['ALL'] + QUARTILES
    classes = [c for c in CATEGORIES if c in set(d.leafcutter2_category)]
    gencode = list(dict.fromkeys(d.sort_values('stack_order').gencode_annotation))
    if palette is None:
        base = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
                '#e6ab02', '#a6761d', '#666666']
        palette = {g: base[i % len(base)] for i, g in enumerate(gencode)}

    fig, axes = plt.subplots(1, len(classes), figsize=(3.6 * len(classes), 3.6),
                             dpi=300, squeeze=False)
    for ax, cls in zip(axes[0], classes):
        sub = d[d.leafcutter2_category == cls].set_index(['quartile', 'gencode_annotation'])
        bottom = np.zeros(len(quartiles))
        for g in gencode:
            vals = np.array([sub['pct_of_class'].get((q, g), 0.0) for q in quartiles])
            cnts = np.array([sub['n_junctions'].get((q, g), 0) for q in quartiles])
            if not vals.any():
                continue
            ax.bar(range(len(quartiles)), vals, bottom=bottom, width=0.72,
                   color=palette[g], label=g, edgecolor='white', linewidth=0.4)
            for i, (v, c) in enumerate(zip(vals, cnts)):
                if v >= label_min_pct:
                    ax.text(i, bottom[i] + v / 2, f'{v:.1f}%\n(n={c:,})', ha='center',
                            va='center', fontsize=5.5, color='white', fontweight='bold')
            bottom += vals
        ax.set_xticks(range(len(quartiles)), quartiles, size=7)
        ax.set_ylim(0, 100)
        ax.set_title(CATEGORY_LABELS.get(cls, cls), size=9)
        ax.set_xlabel('Usage quartile', size=8)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0][0].set_ylabel('% of junctions in class')
    for ax in axes[0][1:]:
        ax.set_yticklabels([])
    axes[0][-1].legend(frameon=False, fontsize=6, loc='center left',
                       bbox_to_anchor=(1.02, 0.5), title='GENCODE v46',
                       title_fontsize=6)
    fig.tight_layout()
    return fig, axes[0]


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



def plot_fig1i(panels, ncols=2):
    """Cumulative distribution of unproductive splicing, split two ways.

    Left: clusters binned by the length of their most-used productive intron.
    Right: genes binned by expression quintile. Both use a log x axis labelled
    in percent, and a light-to-dark ramp within each panel so the ordering of
    the bins is readable without the legend.
    """
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.4),
                             dpi=300, squeeze=False)
    for ax, p in zip(axes[0], panels):
        shades = plt.cm.Blues(np.linspace(0.35, 0.95, len(p['series'])))
        for s, c in zip(p['series'], shades):
            ax.step(s['x'], s['y'], where='post', color=c, linewidth=1.8,
                    label=f"{s['group']} (n = {s['n']:,})")
        ax.set_xscale('log')
        ax.set_xticks([0.001, 0.01, 0.1, 1], ['0%', '1%', '10%', '100%'])
        ax.set_xlabel('Unproductive splicing')
        ax.set_ylim(0, 1)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=6, loc='upper left',
                  title=p['legend_title'], title_fontsize=6.5)
        ax.text(0.98, 0.03,
                r'Spearman $\rho = $' + f"{p['rho']:.2f}" + '\n' + fmt_p(p['pvalue'])
                + '\n' + f"n = {p['n']:,}",
                transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5)
    axes[0][0].set_ylabel('Cumulative distribution')
    for ax in axes[0][1:]:
        ax.set_yticklabels([])
    fig.tight_layout()
    return fig, axes[0]


def save_panel(name, plots_dir='plots', dpi=300):
    import os
    os.makedirs(plots_dir, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        plt.savefig(f'{plots_dir}/{name}.{ext}', dpi=dpi, bbox_inches='tight')
