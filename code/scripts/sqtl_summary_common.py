"""
Shared helpers for the two sQTL summary tables in analysis_files/.

Lifted verbatim from ../analysis/QTL_analysis.ipynb, where these tables were
built by cells whose `to_csv` calls were commented out and whose output was
moved into analysis_files/ by hand. See code_report.txt.

Used by MakeTotalsQTLTable.py (Fig 4a) and MakesQTLStatsTable.py (Fig 4c).
"""

import gzip

import numpy as np
import pandas as pd

# Unproductive-splicing phenotype tables live in the separate SpliFi project.
NOISY_PHENO_DIR = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx'


def get_intron_list(tissue, noisy_pheno_dir=NOISY_PHENO_DIR):
    """Median read count per intron, with its cluster, from the leafcutter2 noise table."""
    intron_list = []
    cluster_list = []
    median_count_list = []
    with gzip.open(
        f'{noisy_pheno_dir}/{tissue}/leafcutter_perind_numers.counts.noise_by_intron.gz'
    ) as fh:
        fh.readline()
        for line in fh:
            line = line.decode().rstrip().split(' ')
            intron_list.append(line[0])
            cluster_list.append(line[0].split(':')[0] + ':' + line[0].split(':')[-2])
            median_count_list.append(int(np.median([int(y) for y in line[1:]])))

    intron_median_counts = pd.DataFrame()
    intron_median_counts['phenotype_id'] = intron_list
    intron_median_counts['cluster'] = cluster_list
    intron_median_counts['median_counts'] = median_count_list
    intron_median_counts['intron'] = intron_median_counts.phenotype_id.apply(
        lambda x: ':'.join(x.split(':')[:-1]))
    return intron_median_counts


def get_perm_counts(tissue, perm_file, noisy_pheno_dir=NOISY_PHENO_DIR):
    """sQTL permutation pass joined to intron counts, with intron/cluster PSI columns."""
    perm = pd.read_csv(perm_file, sep='\t')
    intron_median_counts = get_intron_list(tissue, noisy_pheno_dir)
    perm_counts = pd.merge(perm, intron_median_counts,
                           right_on=['phenotype_id', 'intron', 'cluster'],
                           left_on=['phe_id', 'intron', 'cluster'])

    counts_per_cluster = perm_counts.groupby('cluster').median_counts.sum().reset_index()
    counts_per_cluster.columns = ['cluster', 'cluster_counts']
    counts_per_gene = perm_counts.groupby('gene_id').median_counts.sum().reset_index()
    counts_per_gene.columns = ['gene_id', 'gene_counts']

    perm_counts = pd.merge(perm_counts, counts_per_cluster, on='cluster')
    perm_counts = pd.merge(perm_counts, counts_per_gene, on='gene_id')
    perm_counts['intron_psi'] = perm_counts.median_counts / (perm_counts.cluster_counts + 1e-10)
    perm_counts['cluster_psi'] = perm_counts.cluster_counts / (perm_counts.gene_counts + 1e-10)

    UP_per_cluster = perm_counts.groupby('cluster').itype.apply(
        lambda x: (x == 'UP').sum()).reset_index()
    UP_per_cluster.columns = ['cluster', 'UP_juncs']
    return pd.merge(perm_counts, UP_per_cluster, on='cluster')


def run_tabix_on_nom(nom_file, coords):
    """QTLtools nominal-pass records overlapping `coords` ('chrom:start:end')."""
    import tabix

    chrom, start, end = coords.split(':')
    tb = tabix.open(nom_file)
    nom = tb.query(chrom, int(start), int(end))

    columns = ['#phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis',
               'dist_phe_var', 'var_id', 'var_chr', 'var_from', 'var_to', 'nom_pval',
               'r_squared', 'slope', 'slope_se', 'best_hit']

    nom_bed, nom_bed_cols = [], []
    for idx, record in enumerate(nom):
        nom_bed.append(pd.Series(record))
        nom_bed_cols.append(f'record{str(idx)}')
    nom_bed = pd.concat(nom_bed, axis=1)
    nom_bed.columns = nom_bed_cols
    nom_bed = nom_bed.T
    nom_bed.columns = columns
    nom_bed['#phe_id'] = nom_bed['#phe_id'].apply(lambda x: x.split('.')[0])
    return nom_bed


def get_expression_effect(perm_qtl, eqtl_nom_dir):
    """eQTL nominal record for the gene / variant of one sQTL permutation hit."""
    var_id = perm_qtl.var_id
    chrom = perm_qtl.var_chr
    coords = f'{chrom}:{int(perm_qtl.var_from) - 1}:{int(perm_qtl.var_to) + 1}'
    gene = perm_qtl.gene_id.split('.')[0]

    nom_df = run_tabix_on_nom(f'{eqtl_nom_dir}/{chrom}.txt.gz', coords)
    return nom_df.loc[(nom_df['#phe_id'] == gene) & (nom_df['var_id'] == var_id)]


def get_sqtl_expression_effects(perm_counts, eqtl_nom_dir, ctype='PR,UP', itype='UP',
                                qmax=1e-2, clu_psi_min=0.1):
    """Matched (sQTL beta, eQTL beta, eQTL nominal p) for every significant sQTL."""
    sqtl_slopes, eqtl_slopes, eqtl_pvals = [], [], []

    perm_select = perm_counts.loc[(perm_counts.ctype == ctype) & (perm_counts.itype == itype)
                                  & (perm_counts.q <= qmax)
                                  & (perm_counts.cluster_psi >= clu_psi_min)]

    for idx, row in perm_select.iterrows():
        try:
            eqtl_df = get_expression_effect(row, eqtl_nom_dir)
            if eqtl_df.shape[0] > 0:
                sqtl_slopes.append(float(row.slope))
                eqtl_slopes.append(float(eqtl_df.iloc[0].slope))
                eqtl_pvals.append(float(eqtl_df.iloc[0].nom_pval))
        except:
            continue

    return sqtl_slopes, eqtl_slopes, eqtl_pvals


def p_lambda(p_values, df=1, q=0.5):
    """Genomic-control inflation factor from a distribution of p-values."""
    from scipy.stats import chi2

    if len(p_values) == 0:
        raise ValueError("The p-values list cannot be empty.")
    p_values = np.asarray(p_values)
    if np.any((p_values < 0) | (p_values > 1)):
        raise ValueError("All p-values must be between 0 and 1.")

    chi_squared_stats = chi2.isf(p_values, df)
    median_observed = np.quantile(chi_squared_stats, q)
    median_expected = chi2.ppf(q, df)
    return median_observed / median_expected
