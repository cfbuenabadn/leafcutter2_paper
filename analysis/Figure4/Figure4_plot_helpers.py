"""
Plotting functions for Figure 4.

Created with Figure code cleaner.
Source notebooks: ../QTL_analysis.ipynb, ../hyprcoloc_results.ipynb,
                  ../coloc_plots.ipynb, ../Fig4_example.ipynb

Every function here takes fully plot-ready data (as produced by
`Figure4_helpers.run_all`) and draws one panel. No data loading, no
computation. Function names follow the manuscript caption lettering.
"""

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def plot_fig4A(fig4A_data):
    """fig4A (panel a) -- p-sQTL and u-sQTL counts per tissue, with sample sizes.

    Parameters
    ----------
    fig4A_data : dict
        Keys 'pr_counts', 'up_counts', 'n_samples', 'colors', 'tissue_names',
        all in the same tissue order (total sQTLs, descending).
    """
    pr_counts = fig4A_data['pr_counts']
    up_counts = fig4A_data['up_counts']
    n_samples = fig4A_data['n_samples']
    c_list = fig4A_data['colors']
    tissue_names = fig4A_data['tissue_names']

    fig, axes = plt.subplots(figsize=(10, 4.5), nrows=3,
                             gridspec_kw={'height_ratios': [15, 5, 1], 'hspace': 0.005},
                             dpi=300)
    axes[0].scatter(np.arange(49) + 0.5, pr_counts, c='tab:blue', label='p-sQTLs', zorder=1)

    axes[0].scatter(np.arange(49) + 0.5, up_counts, c='tab:red', label='u-sQTLs', zorder=2)

    for i in np.arange(49):
        axes[0].plot([i + 0.5, i + 0.5], [0, pr_counts[i]],
                     linestyle='--', c='tab:gray', zorder=0, alpha=0.5)
        axes[1].text(i + 0.2, 80, str(n_samples[i]), rotation=90, fontsize=8)

    axes[1].bar(np.arange(49) + 0.5, n_samples, width=0.8, color='lightgray')

    axes[2].bar(np.arange(49) + 0.5, [1] * 49, width=1, color=c_list)

    axes[0].set_xticks([])
    axes[1].set_xticks([])
    axes[1].set_ylabel('Samples')
    axes[2].set_xticks([])
    axes[2].spines[['top', 'bottom', 'right', 'left']].set_visible(False)
    axes[2].set_yticks([])

    axes[2].set_xlim([-1, 49.5])
    axes[1].set_xlim([-1, 49.5])
    axes[0].set_ylim([0, 7300])
    axes[0].set_xlim([-1, 49.5])
    axes[0].set_ylabel('Total sQTLs')
    axes[0].legend(frameon=False)
    axes[0].spines[['top', 'right']].set_visible(False)

    axes[2].set_xticks(np.arange(49) + 0.5, tissue_names, rotation=45, ha='right',
                       rotation_mode='anchor', size=7)

    return axes


def plot_fig4C(fig4C_data):
    """fig4C (panel c) -- per-tissue Spearman rho between sQTL and eQTL betas.

    Marker size encodes the u-sQTL p-value bin; crosses mark non-significant
    tissues. `fig4C_data['bins']` is (n1, n2, n3, n4), counting tissues with
    p <= 1e-10, <= 1e-4, <= 0.05 and > 0.05 respectively.
    """
    rho_u = fig4C_data['rho_u']
    rho_pp = fig4C_data['rho_pp']
    c_list = fig4C_data['colors']
    tissue_names = fig4C_data['tissue_names']
    n1, n2, n3, n4 = fig4C_data['bins']

    fig, axes = plt.subplots(figsize=(3.5, 6), ncols=2,
                             gridspec_kw={'width_ratios': [1, 15], 'wspace': 0}, dpi=300)

    axes[1].scatter(rho_u[:n4], np.arange(n4) + 0.5, c='tab:red', marker='x')
    axes[1].scatter(rho_pp[:n4], np.arange(n4) + 0.5, c='tab:blue', marker='x')

    axes[1].scatter(rho_u[n4:(n4 + n3)], np.arange(n3) + 0.5 + n4, c='tab:red', s=20)
    axes[1].scatter(rho_pp[n4:(n4 + n3)], np.arange(n3) + 0.5 + n4, c='tab:blue', marker='x')

    axes[1].scatter(rho_u[(n4 + n3):(n4 + n3 + n2)], np.arange(n2) + 0.5 + n4 + n3,
                    c='tab:red', s=60)
    axes[1].scatter(rho_pp[(n4 + n3):(n4 + n3 + n2)], np.arange(n2) + 0.5 + n4 + n3,
                    c='tab:blue', marker='x')

    axes[1].scatter(rho_u[(n4 + n3 + n2):], np.arange(n1) + 0.5 + n4 + n3 + n2,
                    c='tab:red', s=120)
    axes[1].scatter(rho_pp[(n4 + n3 + n2):], np.arange(n1) + 0.5 + n4 + n3 + n2,
                    c='tab:blue', marker='x')

    axes[0].barh(np.arange(49) + 0.5, [1] * 49, height=1, color=c_list)

    axes[1].plot([0, 0], [0, 49], c='tab:gray', alpha=0.9, linestyle='--')

    axes[1].scatter([0], [100], c='tab:gray', s=120, label=r'$\leq$ 1e-10')
    axes[1].scatter([0], [100], c='tab:gray', s=60, label=r'$\leq$ 1e-4')
    axes[1].scatter([0], [100], c='tab:gray', s=20, label=r'$\leq$ 0.05')
    axes[1].scatter([0], [100], c='tab:gray', marker='x', label=r'> 0.05')

    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].spines[['top', 'bottom', 'right', 'left']].set_visible(False)
    axes[1].set_yticks([])

    axes[0].set_ylim([-1, 49.5])
    axes[1].set_ylim([-1, 49.5])
    axes[1].set_xlabel(r'Spearman $\rho$' + '\nsQTL vs eQTL' + r' $\beta$')

    axes[1].legend(frameon=False, title='p-value')
    axes[1].spines[['top', 'right']].set_visible(False)

    axes[0].set_yticks(np.arange(49) + 0.5, tissue_names, rotation=0, ha='right',
                       rotation_mode='anchor', size=7)

    return axes


def format_rho_annotation(stats, p_digits=2, include_n=True):
    """On-panel annotation carrying the exact rho and P value.

    Replaces the rounded, hand-typed r'$\\rho = -0.33$' + 'p-val $< 1e^{-36}$'
    strings in the original notebooks, which the editor asked to be given as
    exact values. `include_n` adds a third line with the exact n -- the number of
    points entering the test, not a replicate count.
    """
    rho = stats['rho']
    p = stats['pvalue']
    if p == 0:
        p_str = r'$P < 10^{-300}$'
    else:
        exponent = int(np.floor(np.log10(p)))
        mantissa = p / (10 ** exponent)
        if -3 < exponent < 3:
            # same 3-significant-figure rule as fmt_p / Fig5D_annotated
            p_str = rf'$P = {p:.3g}$'
        else:
            p_str = rf'$P = {mantissa:.{p_digits}f} \times 10^{{{exponent}}}$'
    out = rf'$\rho = {rho:.3f}$' + '\n' + p_str
    if include_n:
        # n outside math mode: mathtext inserts a thin space after a comma,
        # which turns "1,391" into "1, 391".
        out += '\n' + f'$n$ = {stats["n"]:,}'
    return out


def plot_beta_scatter(fit, color, xlim, ylim, text, text_xy, xlabel, ylabel, title,
                      alpha=0.1, figsize=(3, 3), x_fit=(-3.2, 3.2), band='subsample',
                      stats=None, line=None, annotate_n=True):
    """sQTL-beta vs eQTL-beta scatter: fig4B_usQTL, fig4B_psQTL, fig4_colocs_sc2.

    `fit` is the dict returned by `Figure4_helpers.make_rho_fit`. `ylim=None`
    leaves the y range on autoscale, as fig4B_psQTL does.

    `line` selects the dashed line:

      * `'ols'` -- conventional least-squares fit of y on x. Used by
        fig4B_usQTL and fig4B_psQTL.
      * `'rho'` -- slope = Spearman's rho of the full sample, intercept =
        mean(y). Not a regression fit; this is the construction the original
        notebooks drew, kept so the hyprcoloc panel still reproduces.

    `band` selects the shaded region. Each option carries its own legend wording
    in `fit['band']`:

      * `None` (or `'none'`) -- no shaded region. Used by fig4B_usQTL and
        fig4B_psQTL: with the line reported as a least-squares fit and the
        association tested by Spearman's rho, a band adds nothing that the
        quoted rho, P and CI do not already state.
      * `'subsample'` -- 10th-90th percentile of rho over 100 random subsamples
        of 100 points. A dispersion band at an arbitrary n, not an uncertainty
        interval; reproduces the original panels.
      * `'bootstrap95'` -- genuine 95% CI on rho (percentile bootstrap).
      * `'ols_ci95'` / `'ols_pi95'` -- least-squares fit with the 95% CI of the
        fitted mean, or the 95% prediction interval for a single observation.

    When `line` is left as `None` it is inferred from `band`, so existing calls
    keep their previous behaviour.

    Pass `stats` (from `Figure4_helpers.make_rho_stats`) to have the on-panel
    annotation written from the exact rho, P and n instead of the `text` string;
    `annotate_n=False` drops the n line. A stats-driven annotation is anchored by
    its bottom edge, so it grows upward from `text_xy` and stays inside the axes
    whatever the line count; a literal `text` keeps matplotlib's default baseline
    anchoring, so panels passing `text` are unchanged.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    plt.scatter(fit['x'], fit['y'], facecolor=color, edgecolor='none', alpha=alpha,
                rasterized=True)

    if band == 'none':
        band = None
    if band not in (None, 'subsample', 'bootstrap95', 'ols_ci95', 'ols_pi95'):
        raise ValueError("band must be None, 'subsample', 'bootstrap95', "
                         "'ols_ci95' or 'ols_pi95'")

    if line is None:
        line = 'ols' if band in ('ols_ci95', 'ols_pi95') else 'rho'
    if line not in ('ols', 'rho'):
        raise ValueError("line must be 'ols' or 'rho'")

    if line == 'ols':
        ax.plot(fit['x_grid'], fit['ols']['fitted'], color='black', linestyle='--')
    else:
        x = np.linspace(x_fit[0], x_fit[1], 100)
        ax.plot(x, fit['slope'] * x + fit['intercept'], color='black', linestyle='--')

    if band in ('subsample', 'bootstrap95'):
        hi_slope, lo_slope = ((fit['slope_90'], fit['slope_10']) if band == 'subsample'
                              else (fit['boot_97p5'], fit['boot_2p5']))
        x_b = np.linspace(x_fit[0], x_fit[1], 100)
        ax.fill_between(x_b, (hi_slope * np.array(x_b)) + fit['intercept'],
                        (lo_slope * x_b) + fit['intercept'],
                        color='tab:gray', alpha=0.2, label='Confidence Interval',
                        linewidth=0)
    elif band in ('ols_ci95', 'ols_pi95'):
        ols = fit['ols']
        lo, hi = ((ols['ci_lo'], ols['ci_hi']) if band == 'ols_ci95'
                  else (ols['pi_lo'], ols['pi_hi']))
        ax.fill_between(fit['x_grid'], hi, lo, color='tab:gray', alpha=0.2, linewidth=0,
                        label='95% CI' if band == 'ols_ci95' else '95% PI')

    ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.spines[['top', 'right']].set_visible(False)

    if stats is not None:
        text = format_rho_annotation(stats, include_n=annotate_n)
        text_va = 'bottom'
    else:
        text_va = 'baseline'
    ax.text(text_xy[0], text_xy[1], text, va=text_va)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax


def plot_fig4D(fig4D_data):
    """fig4D (panel d) -- percentage of GWAS loci colocalizing with a u-sQTL.

    Boxes are the per-tissue percentages for each trait; the red dots are the
    percentage of loci colocalizing in any tissue.
    """
    boxes = fig4D_data['boxes']
    overall_pct = fig4D_data['overall_pct']
    trait_names = fig4D_data['trait_names']

    fig, ax = plt.subplots(figsize=(4, 2), dpi=300)
    sns.boxplot(data=boxes, ax=ax, boxprops=dict(facecolor='white', edgecolor='black'))
    ax.scatter(range(19), overall_pct, c='tab:red', s=20)

    ax.set_xticks(range(19), trait_names, rotation=45, ha='right',
                  rotation_mode='anchor', size=7)
    ax.set_ylabel('Percentage of GWAS loci\ncolocalized w/ non p-sQTLs', size=8)
    ax.spines[['top', 'right']].set_visible(False)

    return ax


def fmt_p(p):
    """P-value formatting used by Fig5D_annotated in Figure5.ipynb."""
    if p >= 1e-3:
        return f'{p:.3g}'
    m, e = f'{p:.2e}'.split('e')
    return f'{m}$\\times$10$^{{{int(e)}}}$'


def fmt_qtl_annotation(beta, pval, p_label='P'):
    """'beta = x.xx / P = ...' block, in the Fig5D_annotated format.

    `p_label` is the only thing that changes between reporting the nominal P at
    the plotted variant ('P', the default here) and the phenotype-level
    permutation-adjusted P ('adj. P', which is what Fig5D_annotated shows).
    """
    if beta is None or pval is None:
        return f'$\\beta$ = ??\n{p_label} = ??'
    return f'$\\beta$ = {beta:.2f}\n{p_label} = {fmt_p(pval)}'


def plot_genotype_boxplot(df, color='black', beta=None, pval=None, p_label='P',
                          headroom=0.20):
    """fig4E_sQTL / fig4E_eQTL -- genotype boxplot + strip plot of a phenotype.

    `df` is the table built by `Figure4_helpers.make_genotype_boxplot_df`
    (columns 'genotype' and 'qqnorm', samples ordered hom-ref/het/hom-alt).

    Passing `beta` and `pval` annotates the panel in the same format as
    `Fig5D_annotated` in Figure5.ipynb: the per-genotype n goes under each x tick
    label, and the effect size with its P value is centred at the top of the
    axes, with `headroom` of extra y range added to make space.

    Both should describe the **plotted variant**, so that the annotation matches
    the genotypes on the x axis; pass `p_label='adj. P'` with a phenotype-level
    permutation-adjusted P only if that is what you intend to show.
    """
    fig, ax = plt.subplots(figsize=(2.5, 2.5), dpi=300)

    sns.boxplot(x='genotype', y='qqnorm', data=df, width=.8, showfliers=False,
                linewidth=2, boxprops=dict(facecolor='white', edgecolor=color),
                whiskerprops=dict(color=color),
                capprops=dict(color=color),
                medianprops=dict(linewidth=2, color='black'), ax=ax)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.stripplot(x='genotype', y='qqnorm', data=df, alpha=0.9, edgecolor=color,
                  facecolor='none', linewidth=1, ax=ax, rasterized=True)

    if beta is not None or pval is not None:
        counts = df.genotype.value_counts()
        order = [t.get_text() for t in ax.get_xticklabels()]
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f'{g}\nn={counts[g]}' for g in order])

        ybottom, ytop = ax.get_ylim()
        ax.set_ylim(ybottom, ytop + headroom * (ytop - ybottom))
        ax.text(0.5, 0.99, fmt_qtl_annotation(beta, pval, p_label),
                transform=ax.transAxes, ha='center', va='top', fontsize=8)

    return ax


def plot_isoform_annotations(annotation_exons, gene, colores=None, start=None, end=None,
                             figsize=None, lwidth=5, iso_order=None, axes=None, xlim=None):
    """fig4E_ASB16 -- one row per transcript: exon blocks joined by a line."""
    gene_exons = annotation_exons.loc[annotation_exons.gene_id == gene]

    if iso_order is None:
        isoforms = sorted(gene_exons.transcript_id.unique())
    else:
        isoforms = iso_order

    isoform_dict = {}
    for i, iso in enumerate(isoforms):
        isoform_name = f'isoform_{str(i+1)}'
        df = gene_exons.loc[gene_exons.transcript_id == iso].copy()
        df['transcript_id'] = f'{gene}.{isoform_name}'
        isoform_dict.update({isoform_name: {'df': df}})

    try:
        chrom = list(annotation_exons.chrom)[0]
    except:
        chrom = list(annotation_exons['#chrom'])[0]
    if start is None:
        start = str(np.min([int(list(gene_exons.start)[0]), int(list(gene_exons.start)[0])]) - 1000)
    if end is None:
        end = str(np.max([int(list(gene_exons.end)[-1]), int(list(gene_exons.end)[-1])]) + 1000)

    coords = [f'{chrom}:{start}', f'{chrom}:{end}']

    plot_gene_isoforms(isoform_dict, coords, color_list=colores, figsize=figsize,
                       lwidth=lwidth, axes=axes, xlim=xlim)


def plot_gene_isoforms(isoforms_dict, coordinates, color_list=None, axes=None,
                       figsize=None, lwidth=5, xlim=None):
    if xlim is None:
        xlim1 = int(coordinates[0].split(':')[1])
        xlim2 = int(coordinates[-1].split(':')[1])
    else:
        xlim1 = xlim[0]
        xlim2 = xlim[1]

    if color_list is None:
        color_list = sns.color_palette("tab10")

    K = len(isoforms_dict)

    if figsize is None:
        figsize = (20, 3)

    if axes is None:
        fig, axes = plt.subplots(K, 1, figsize=figsize)

    for i in range(K):
        isoform_df = isoforms_dict[f'isoform_{str(i+1)}']['df']
        ax = axes[i]
        color = color_list[i]
        plot_isoform(isoform_df, ax, color, xlim1, xlim2, lwidth=lwidth)


def plot_isoform(isoform_df, ax, color, xlim1, xlim2, lwidth):
    is_first = True
    for idx, row in isoform_df.iterrows():
        start = int(row.start)
        end = int(row.end)
        if is_first:
            first = end
            is_first = False

        ax.fill_between([start, end], [0, 0], [1, 1], color=color, zorder=2)

    ax.plot([first, start], [0.5, 0.5], c=color, linewidth=lwidth)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[['bottom', 'top', 'right', 'left']].set_visible(False)
    ax.set_xlim([xlim1, xlim2])


def plot_fig4F_locuscompare(fig4F_data):
    """fig4F (panel f) -- LocusCompare panels at the ASB16 / bipolar-disorder locus.

    Three scatters of -log10 P: u-sQTL vs GWAS, eQTL vs GWAS, and u-sQTL vs
    eQTL, each coloured by r^2 with the lead variant rs7212573, which is ringed
    in black. Only the first panel draws in r^2 order (low r^2 underneath); the
    other two keep the SNP order, as in Fig4_example.ipynb.
    """
    pts = fig4F_data['points']
    pts_sorted = fig4F_data['points_sorted']
    lead = fig4F_data['lead']

    fig, ax = plt.subplots(ncols=3, figsize=(12, 3), width_ratios=[1, 1, 1.2], dpi=300)

    ax[0].scatter(pts_sorted.sQTL, pts_sorted.GWAS, c=pts_sorted.r2, s=15,
                  cmap='inferno', vmin=0, rasterized=True)
    ax[0].set_xlabel('ASB16 u-sQTL, -log10(P)\nBrain - Cerebellar Hemisphere')
    ax[0].set_ylabel('GWAS, -log10(P)\nBipolar disorder')
    ax[0].scatter([lead['sQTL']], [lead['GWAS']], edgecolor='black', facecolor='none', s=50)
    ax[0].text(7, 7.5, lead['label'])

    ax[1].scatter(pts.eQTL, pts.GWAS, c=pts.r2, s=15, cmap='inferno', vmin=0,
                  rasterized=True)
    ax[1].set_xlabel('ASB16 eQTL, -log10(P)\nBrain - Cerebellar Hemisphere')
    ax[1].set_ylabel('GWAS, -log10(P)\nBipolar disorder')
    ax[1].scatter([lead['eQTL']], [lead['GWAS']], edgecolor='black', facecolor='none', s=50)
    ax[1].text(4.5, 7.5, lead['label'])

    c = ax[2].scatter(pts.sQTL, pts.eQTL, c=pts.r2, s=15, cmap='inferno', vmin=0,
                      rasterized=True)
    plt.colorbar(c, label=r'$R^2$')
    ax[2].set_xlabel('ASB16 u-sQTL, -log10(P)\nBrain - Cerebellar Hemisphere')
    ax[2].set_ylabel('ASB16 eQTL, -log10(P)\nBrain - Cerebellar Hemisphere')
    ax[2].scatter([lead['sQTL']], [lead['eQTL']], edgecolor='black', facecolor='none', s=50)
    ax[2].text(7, 6.5, lead['label'])

    for a in ax:
        a.spines[['top', 'right']].set_visible(False)

    fig.subplots_adjust(wspace=0.35, hspace=5)
    return ax
