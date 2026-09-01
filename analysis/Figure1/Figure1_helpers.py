"""Data for the Figure 1 panels that can be rebuilt from the LeafCutter2 run.

Figure 1 was made in a different project directory by Yang Li and Quinn Hauck.
The code was handed over as a bundle (yang_lc2/leafcutter2_paper_figures/), and
the pipeline outputs it consumes still live in Yang's directory. This module
reads those outputs directly and rebuilds the panels that are recoverable:

    1c   the BED12 junction tracks behind the SRSF4 browser view.  The panel
         itself is an IGV screenshot, so only the track data is recovered.
    1d   junction classification by usage quartile (stacked bars)
    1g   log2 fold-change ECDFs, unproductive vs productive, four perturbations
    1h   the three NMD efficiency rules, as tested on naRNA vs polyA

Not recoverable here:
    1a, 1b   schematics
    1e, 1f   1e is the GENCODE composition panel (source present, see
             fig1e_data); 1f is the simulation benchmark, which lives in
             bfairkun/20260825_leaf2simulation_paper
    1i       nothing in either directory produces it

Everything is read-only from SOURCE; nothing is written back there.
"""
import os
import pickle

import numpy as np
import pandas as pd

# The upstream project. Read-only.
SOURCE = '/project2/yangili1/yangili/chRNA/nmd_splicing_rules/add_on_script'
PARSED = f'{SOURCE}/parse_classifications_output'
CLUSTERING = f'{SOURCE}/leaf2_2026/clustering'

JUNCTION_CLASSIFICATIONS = f'{CLUSTERING}/leafcutter2_junction_classifications.txt'
NUC_RULE_DISTANCES = f'{CLUSTERING}/leafcutter2_nuc_rule_distances.txt'
EXON_STATS = f'{CLUSTERING}/leafcutter2_exon_stats.txt'
LOG2FC_PICKLE = f'{SOURCE}/intron_log2fc.pickle'

# Panel 1h's productive column is a matched coding-junction control set, built by
# a fork of the classifier (pc_PosControlClassifier.py). It lives in Quinn
# Hauck's directory; 260128_pos_control is the set the published panel used.
POS_CONTROL = ('/project2/yangili1/qhauck/nmd_splicing_rules/add_on_script/'
               'rules_test/260128_pos_control')

QUARTILES = ['Q1', 'Q2', 'Q3', 'Q4']
CATEGORIES = ['productive', 'unproductive', 'utr']
CATEGORY_LABELS = {'productive': 'Productive', 'unproductive': 'Unproductive',
                   'utr': 'Near-UTR'}
# from the handover bundle's log2fc_common.CLASS_COLORS
CATEGORY_COLORS = {'productive': '#377eb8', 'unproductive': '#e41a1c',
                   'utr': '#999999'}

# get_default_comparisons() in the bundle, in the order the panel shows them
COMPARISONS = [
    (('Colombo.SM67.and.SMG7', 'Colombo.Control'), 'SMG6/SMG7 dKD vs control (HeLa)'),
    (('naRNA', 'SteadyState'), 'naRNA-seq vs polyA RNA-seq (LCL)'),
    (('Darman.Darman.Cycloxhexamide', 'Darman.DMSO.Control'), 'Cycloheximide vs DMSO (NALM6)'),
    (('LCL.Monosome', 'LCL.HeavyPolysome'), 'Monosome vs heavy polysome (LCL)'),
]


def junctions_table(quartile='ALL'):
    """One row per junction: coordinates, gene, LeafCutter2 class, usage."""
    return pd.read_csv(f'{PARSED}/junctions_{quartile}.tsv', sep='\t')


# --------------------------------------------------------------------------- #
# 1c -- BED12 tracks behind the SRSF4 browser view
# --------------------------------------------------------------------------- #

BED12_COLUMNS = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                 'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                 'blockSizes', 'blockStarts']


def make_fig1c_tracks(gene='SRSF4'):
    """The BED12 junction tracks, restricted to the panel's gene.

    Panel 1c is an IGV screenshot, so this recovers the underlying track rather
    than the image. One table per usage quartile plus ALL, each carrying the
    LeafCutter2 class of every junction drawn, since the track is coloured by it.
    """
    out = {}
    for q in ['ALL'] + QUARTILES:
        bed = pd.read_csv(f'{PARSED}/junctions_{q}.bed12', sep='\t',
                          header=None, names=BED12_COLUMNS)
        cls = junctions_table(q)[['coordinates', 'gene', 'leafcutter2_category',
                                  'gencode_annotation', 'usage']]
        merged = bed.merge(cls, left_on='name', right_on='coordinates', how='left')
        out[q] = merged.loc[merged.gene == gene].reset_index(drop=True) if gene else merged
    return out


# --------------------------------------------------------------------------- #
# 1d / 1e -- classification by usage quartile, and GENCODE composition
# --------------------------------------------------------------------------- #

def make_fig1d_data():
    """Counts and percentages of each LeafCutter2 class, per usage quartile.

    ALL is included as its own bar, as in the published panel, and comes first.
    """
    rows = []
    for q in ['ALL'] + QUARTILES:
        t = junctions_table(q)
        n = len(t)
        for cat in CATEGORIES:
            k = int((t.leafcutter2_category == cat).sum())
            rows.append({'quartile': q, 'category': cat,
                         'category_label': CATEGORY_LABELS[cat],
                         'color': CATEGORY_COLORS[cat],
                         'n_junctions': k, 'pct_junctions': 100 * k / n,
                         'n_in_quartile': n})
    return pd.DataFrame(rows)


def make_fig1e_data():
    """GENCODE v46 transcript-type composition of each class, per usage quartile.

    Restricted to junctions GENCODE v46 annotates: `not_in_gencode` is 94% of
    the unproductive class and would swamp the comparison, and the panel asks
    what GENCODE calls the junctions it does annotate. One row per (quartile,
    LeafCutter2 class, GENCODE type); `pct_of_class` is within that quartile and
    class, so the stacked bars sum to 100.
    """
    rows = []
    for q in ['ALL'] + QUARTILES:
        t = junctions_table(q)
        t = t.loc[t.gencode_annotation != 'not_in_gencode']
        g = (t.groupby(['leafcutter2_category', 'gencode_annotation']).size()
               .rename('n_junctions').reset_index())
        g['pct_of_class'] = 100 * g.n_junctions / g.groupby(
            'leafcutter2_category').n_junctions.transform('sum')
        g.insert(0, 'quartile', q)
        rows.append(g)
    out = pd.concat(rows, ignore_index=True)
    out['category_label'] = out.leafcutter2_category.map(CATEGORY_LABELS)
    # stack order: most abundant GENCODE type at the bottom, pooled over ALL
    order = (out.loc[out.quartile == 'ALL'].groupby('gencode_annotation')
                .n_junctions.sum().sort_values(ascending=False).index.tolist())
    out['stack_order'] = out.gencode_annotation.map({g: i for i, g in enumerate(order)})
    return out.sort_values(['quartile', 'leafcutter2_category', 'stack_order'])


# --------------------------------------------------------------------------- #
# 1g -- log2 fold-change ECDFs
# --------------------------------------------------------------------------- #

def _ecdf(x, steps=200):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    grid = np.linspace(np.nanpercentile(x, 0.5), np.nanpercentile(x, 99.5), steps)
    return grid, np.array([(x <= g).mean() for g in grid])


def make_fig1g_data(comparisons=None):
    """Per comparison and LeafCutter2 class: the log2FCs and their ECDF.

    log2FCs come from intron_log2fc.pickle, the cache written by
    cl_plot_log2fc_dist.py --save-pickle. Classes come from the junction table,
    joined on (chrom, start, end).
    """
    if comparisons is None:
        comparisons = COMPARISONS
    with open(LOG2FC_PICKLE, 'rb') as fh:
        log2fc = pickle.load(fh)

    t = junctions_table('ALL')
    cls = dict(zip(zip(t.chrom, t.start, t.end), t.leafcutter2_category))

    out = []
    for comp, label in comparisons:
        if comp not in log2fc:
            continue
        vals = log2fc[comp]
        by_class = {c: [] for c in CATEGORIES}
        for key, v in vals.items():
            c = cls.get(key)
            if c is not None and np.isfinite(v):
                by_class[c].append(v)
        for c in CATEGORIES:
            arr = np.asarray(by_class[c], dtype=float)
            if arr.size == 0:
                continue
            gx, gy = _ecdf(arr)
            out.append({'comparison': comp, 'comparison_label': label,
                        'category': c, 'category_label': CATEGORY_LABELS[c],
                        'color': CATEGORY_COLORS[c],
                        'log2fc': arr, 'ecdf_x': gx, 'ecdf_y': gy,
                        'n': int(arr.size), 'median_log2fc': float(np.median(arr))})
    return out


# --------------------------------------------------------------------------- #
# 1h -- the NMD efficiency rules
# --------------------------------------------------------------------------- #

RULES = [
    ('ejc_distance', 'PTC further than 50 nt from the last EJC', 50),
    ('Exons_before', 'Number of introns upstream of the PTC', None),
    ('Exons_after', 'Number of introns downstream of the PTC', None),
]


def make_fig1h_data(comparison=('naRNA', 'SteadyState')):
    """One-sided Wilcoxon test per rule, on naRNA-vs-polyA log2FC.

    Reproduces rules_test/260128_updated_intron_fd.Rmd: for each rule the
    unproductive junctions are split into a low and a high group, and the two
    groups' log2FCs are compared with a one-sided Mann-Whitney U test. The
    reported effect is mean(high) - mean(low), the panel's "Delta log2FD".

    The productive column is the negative control: the same test on junctions
    LeafCutter2 calls coding. The EJC rule has no productive counterpart --
    productive junctions have no PTC -- which is why that cell is N/A.
    """
    from scipy.stats import mannwhitneyu

    with open(LOG2FC_PICKLE, 'rb') as fh:
        log2fc = pickle.load(fh)[comparison]

    def _key(coord):
        chrom, rest = coord.split(':')
        a, b = rest.split('-')
        return (chrom, int(a), int(b))

    def _features(jc_path, exon_path, nuc_path=None):
        jc = pd.read_csv(jc_path, sep='\t')
        jc = jc.loc[~jc.UTR].drop_duplicates('Intron_coord')
        f = jc.merge(pd.read_csv(exon_path, sep='\t'),
                     on=['Gene_name', 'Intron_coord'], how='left')
        if nuc_path is not None:
            f = f.merge(pd.read_csv(nuc_path, sep='\t'),
                        on=['Gene_name', 'Intron_coord'], how='left')
        else:
            f['ejc_distance'] = np.nan
        f['log2fc'] = [log2fc.get(_key(c), np.nan) for c in f.Intron_coord]
        return f.loc[np.isfinite(f.log2fc)]

    # The main run's exon-stats file only covers unproductive junctions, so the
    # productive control comes from the matched set built by the classifier fork.
    unproductive = _features(JUNCTION_CLASSIFICATIONS, EXON_STATS, NUC_RULE_DISTANCES)
    productive = _features(f'{POS_CONTROL}/final_pos_control_junction_classifications.txt',
                           f'{POS_CONTROL}/final_pos_control_exon_stats.txt')

    rows = []
    for cls_name, feats in (('Unproductive', unproductive), ('Productive', productive)):
        sub = feats
        for col, label, thresh in RULES:
            if cls_name == 'Productive' and col == 'ejc_distance':
                rows.append({'rule': label, 'class': cls_name, 'n_low': np.nan,
                             'n_high': np.nan, 'mean_low': np.nan, 'mean_high': np.nan,
                             'delta_log2fd': np.nan, 'p_value': np.nan,
                             'note': 'not defined: productive junctions have no PTC'})
                continue
            v = sub.loc[np.isfinite(sub[col])]
            # `high` is always the group the rule predicts is degraded MORE
            # efficiently, so a positive delta means the rule behaves as
            # expected. For the 50-nt rule that is the group whose PTC is
            # FURTHER than 50 nt from the last EJC: a PTC within 50 nt escapes
            # NMD, so proximity predicts less degradation, not more.
            if thresh is not None:
                high, low = v.loc[v[col] > thresh, 'log2fc'], v.loc[v[col] <= thresh, 'log2fc']
            else:
                q1, q3 = v[col].quantile([0.25, 0.75])
                low, high = v.loc[v[col] <= q1, 'log2fc'], v.loc[v[col] >= q3, 'log2fc']
            if len(low) < 3 or len(high) < 3:
                continue
            U, p = mannwhitneyu(high, low, alternative='greater')
            rows.append({'rule': label, 'class': cls_name,
                         'n_low': int(len(low)), 'n_high': int(len(high)),
                         'mean_low': float(low.mean()), 'mean_high': float(high.mean()),
                         'delta_log2fd': float(high.mean() - low.mean()),
                         'U': float(U), 'p_value': float(p), 'note': ''})
    return pd.DataFrame(rows)



# --------------------------------------------------------------------------- #
# 1i -- unproductive splicing vs intron length and vs gene expression
# --------------------------------------------------------------------------- #

# Bin order as drawn, palest to darkest
LENGTH_BINS = ['<1kb', '1-5kb', '5-20kb', '20-50kb', '>50kb']
RPKM_BINS = ['Q1 - lowly expressed', 'Q2', 'Q3', 'Q4', 'Q5 - highly expressed']

# Zeros cannot go on a log axis; the source notebook floors them here.
ZERO_FLOOR = 0.0005


def make_fig1i_data(data_dir='figure_data'):
    """ECDF of unproductive splicing, split by intron length and by expression.

    Reads the two tables written by Figure1_fig1i.R, which ports the
    computation out of the Geuvadis EUR analysis: 373 lymphoblastoid cell
    lines, one row per cluster (length panel) or per gene (expression panel),
    not per sample.
    """
    from scipy.stats import spearmanr

    panels = []
    for name, fname, value_col, group_col, bins, xlab in (
        ('intron length', 'fig1i_intron_length.tsv', 'meanUnprod',
         'intron_length_bin', LENGTH_BINS, 'intron length'),
        ('gene expression', 'fig1i_expression.tsv', 'medUPratio',
         'rpkm_bin', RPKM_BINS, 'gene RPKM quintile'),
    ):
        t = pd.read_csv(os.path.join(data_dir, fname), sep='\t')
        covar = 'intron_length' if 'intron_length' in t.columns else 'meanrpkm'
        rho, p = spearmanr(t[value_col], t[covar])
        series = []
        for b in bins:
            v = t.loc[t[group_col] == b, value_col].to_numpy(dtype=float)
            if v.size == 0:
                continue
            v = np.where(v == 0, ZERO_FLOOR, v)
            x = np.sort(v)
            series.append({'group': b, 'n': int(x.size), 'x': x,
                           'y': np.arange(1, x.size + 1) / x.size})
        panels.append({'panel': name, 'legend_title': xlab, 'series': series,
                       'rho': float(rho), 'pvalue': float(p), 'n': int(len(t))})
    return panels


# --------------------------------------------------------------------------- #

PLOT_READY_VARS = ['fig1c_tracks', 'fig1d_data', 'fig1e_data', 'fig1g_data',
                   'fig1h_data', 'fig1i_data']


def run_all(data_dir='figure_data', gene='SRSF4'):
    os.makedirs(data_dir, exist_ok=True)
    data = {
        'fig1c_tracks': make_fig1c_tracks(gene),
        'fig1d_data': make_fig1d_data(),
        'fig1e_data': make_fig1e_data(),
        'fig1g_data': make_fig1g_data(),
        'fig1h_data': make_fig1h_data(),
        'fig1i_data': make_fig1i_data(data_dir),
    }
    for k, v in data.items():
        with open(os.path.join(data_dir, f'{k}.pickle'), 'wb') as fh:
            pickle.dump(v, fh)
    return data


def load_plot_data(data_dir='figure_data', vars=PLOT_READY_VARS):
    out = {}
    for k in vars:
        with open(os.path.join(data_dir, f'{k}.pickle'), 'rb') as fh:
            out[k] = pickle.load(fh)
    return out
