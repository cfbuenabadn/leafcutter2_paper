"""
Build analysis_files/Total_sQTLs.tsv.gz -- the per-tissue sQTL counts behind
Figure 4a.

One row per GTEx tissue, columns PR / UP / Other / total: the number of
significant sQTLs (permutation-pass Storey q <= 0.1) whose cluster is entirely
productive (PR), mixed productive/unproductive (PR,UP -> 'UP'), or anything else.

Reproduces the commented-out `qtl_dataframe.to_csv('Total_sQTLs.tsv.gz', ...)`
in ../analysis/QTL_analysis.ipynb.

Usage:
    python scripts/MakeTotalsQTLTable.py <out.tsv.gz> <tissue> [<tissue> ...]
"""

import sys

import pandas as pd

from sqtl_summary_common import get_perm_counts

PERM_TEMPLATE = 'results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz'
QMAX = 0.1


def tissue_counts(tissue):
    perm_counts = get_perm_counts(tissue, PERM_TEMPLATE.format(tissue=tissue))
    sig = perm_counts.q <= QMAX
    return {
        'PR': int((sig & (perm_counts.ctype == 'PR')).sum()),
        'UP': int((sig & (perm_counts.ctype == 'PR,UP')).sum()),
        'Other': int((sig & (perm_counts.ctype != 'PR') & (perm_counts.ctype != 'PR,UP')).sum()),
    }


def main():
    out_path, tissues = sys.argv[1], sorted(sys.argv[2:])

    sqtl_dict = {}
    for tissue in tissues:
        sqtl_dict[tissue] = tissue_counts(tissue)
        print(f'  {tissue}: {sqtl_dict[tissue]}', flush=True)

    qtl_dataframe = pd.DataFrame(sqtl_dict).T
    qtl_dataframe['total'] = qtl_dataframe.sum(axis=1)
    qtl_dataframe.to_csv(out_path, sep='\t', header=True, index=True)
    print(f'wrote {out_path}: {qtl_dataframe.shape[0]} tissues')


if __name__ == '__main__':
    main()
