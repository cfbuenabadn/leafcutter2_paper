"""
Build analysis_files/sQTL_stats.tsv.gz -- the per-tissue sQTL-vs-eQTL effect-size
statistics behind Figure 4c (Supplementary Table 9).

For each GTEx tissue and each of three sQTL classes, correlates the sQTL effect
size with the eQTL effect size of the same variant on the host gene:

    u_sqtl   unproductive introns in mixed clusters   (ctype PR,UP / itype UP)
    p_sqtl   productive introns in mixed clusters     (ctype PR,UP / itype PR)
    pp_sqtl  productive introns in productive clusters(ctype PR    / itype PR)

Reports Pearson and Spearman rho with p-values, the genomic-control lambda at
the 50th and 90th percentiles, and n.

Reproduces the commented-out `sqtl_stats.to_csv('sQTL_stats.tsv.gz', ...)` in
../analysis/QTL_analysis.ipynb.

Slow: every significant sQTL costs one tabix query against the eQTL nominal pass.

Usage:
    python scripts/MakesQTLStatsTable.py <out.tsv.gz> <tissue> [<tissue> ...]
"""

import sys

import pandas as pd

from sqtl_summary_common import get_perm_counts, get_sqtl_expression_effects, p_lambda

PERM_TEMPLATE = 'results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz'
EQTL_NOM_TEMPLATE = 'results/eqtl/GTEx/{tissue}/cis_100000/nom'

QMAX = 1e-1
CLU_PSI_MIN = 0.1

# sQTL class -> (ctype, itype), as in QTL_analysis.ipynb
CLASSES = {
    'u_sqtl': ('PR,UP', 'UP'),
    'p_sqtl': ('PR,UP', 'PR'),
    'pp_sqtl': ('PR', 'PR'),
}

COLUMNS = ['tissue', 'sqtl_type', 'lambda_inflation', 'lambda_90', 'n',
           'pearson_pval', 'pearson_rho', 'spearman_pval', 'spearman_rho']


def tissue_stats(tissue):
    from scipy.stats import pearsonr, spearmanr

    perm_counts = get_perm_counts(tissue, PERM_TEMPLATE.format(tissue=tissue))
    eqtl_nom_dir = EQTL_NOM_TEMPLATE.format(tissue=tissue)

    rows = []
    for name, (ctype, itype) in CLASSES.items():
        sqtl_beta, eqtl_beta, eqtl_p = get_sqtl_expression_effects(
            perm_counts, eqtl_nom_dir, ctype=ctype, itype=itype,
            qmax=QMAX, clu_psi_min=CLU_PSI_MIN)

        pear = pearsonr(sqtl_beta, eqtl_beta)
        spear = spearmanr(sqtl_beta, eqtl_beta)
        rows.append({
            'tissue': tissue,
            'sqtl_type': name,
            'lambda_inflation': p_lambda(eqtl_p),
            'lambda_90': p_lambda(eqtl_p, q=0.9),
            'n': float(len(sqtl_beta)),
            'pearson_pval': pear[1],
            'pearson_rho': pear[0],
            'spearman_pval': spear[1],
            'spearman_rho': spear[0],
        })
    return rows


def main():
    out_path, tissues = sys.argv[1], sorted(sys.argv[2:])

    rows = []
    for tissue in tissues:
        rows.extend(tissue_stats(tissue))
        print(f'  {tissue} done ({len(rows)} rows so far)', flush=True)

    # Sorted by tissue then sqtl_type, matching the notebook's pivot output.
    sqtl_stats = pd.DataFrame(rows, columns=COLUMNS)
    sqtl_stats = sqtl_stats.sort_values(['tissue', 'sqtl_type']).reset_index(drop=True)
    sqtl_stats.to_csv(out_path, sep='\t', index=False, header=True)
    print(f'wrote {out_path}: {sqtl_stats.shape[0]} rows '
          f'({sqtl_stats.tissue.nunique()} tissues x {sqtl_stats.sqtl_type.nunique()} classes)')


if __name__ == '__main__':
    main()
