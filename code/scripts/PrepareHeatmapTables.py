import sys
import os
import gzip
import argparse
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import os
from rpy2.robjects import pandas2ri

def process_comparison(rds_file):
    a_v_b = rds_file.split('.rds')[0]
    readRDS = ro.r['readRDS']
    rds = readRDS(f'/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge/{rds_file}')
    
    ds = rds[0]
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(ds)
        output = pd_from_r_df
        
    dge = rds[1]
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(dge)
        dge = pd_from_r_df
        
    intron_list = []
    cluster_list = []
    median_count_list = []
    
    with gzip.open(
        f'/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds/GTEx/{a_v_b}/ds_perind_numers.counts.noise_by_intron.lf1.gz'
    ) as fh:
        fh.readline()
        for line in fh:
            line = line.decode().rstrip().split(' ')
            intron_list.append(line[0])
            cluster_list.append(line[0].split(':')[0] + ':' + line[0].split(':')[-1])
            median_counts = int(np.median([int(y) for y in line[1:]]))
            median_count_list.append(median_counts)
    
    intron_median_counts = pd.DataFrame()
    intron_median_counts['intron'] = intron_list
    intron_median_counts['cluster'] = cluster_list
    intron_median_counts['median_counts'] = median_count_list
    
    ds = pd.merge(output, intron_median_counts, right_on=['intron', 'cluster'], left_on=['intron', 'cluster'])
    
    cluster_counts = pd.DataFrame(ds.groupby('cluster').median_counts.sum()).reset_index()
    cluster_counts.columns = ['cluster', 'cluster_counts']

    ds = pd.merge(ds, cluster_counts, right_on='cluster', left_on='cluster')

    dfCluster_max = pd.DataFrame(ds.groupby('gene_id').cluster_counts.max()).reset_index()
    dfCluster_max.columns = ['gene_id', 'max_cluster_counts']

    ds = pd.merge(ds, dfCluster_max, left_on='gene_id', right_on='gene_id')
    ds['cluster_fraction'] = ds.cluster_counts/ds.max_cluster_counts
    
    ds_dge_table = pd.merge(ds, dge, left_on='gene_id', right_on = 'gene_id').drop_duplicates()
    
    tissue1, tissue2 = a_v_b.split('_v_')
    junc_cols = list(ds_dge_table.columns[:6])
    
    psi_cols = junc_cols + [tissue1, tissue2]
    
    PSI = ds_dge_table[psi_cols]
    
    a_v_b_columns = junc_cols + [a_v_b]
    
    delta_psi = ds_dge_table[junc_cols + ['deltapsi']]
    logFC_psi = ds_dge_table[junc_cols + ['logef']]
    logFC_exp = ds_dge_table[junc_cols + ['log2FoldChange']]
    cluster_counts = ds_dge_table[junc_cols + ['max_cluster_counts']]
    cluster_fraction = ds_dge_table[junc_cols + ['cluster_fraction']]
    psi_p = ds_dge_table[junc_cols + ['p']]
    exp_p = ds_dge_table[junc_cols + ['pvalue']]
        

    delta_psi.columns = a_v_b_columns
    logFC_psi.columns = a_v_b_columns
    logFC_exp.columns = a_v_b_columns
    cluster_counts.columns = a_v_b_columns
    cluster_fraction.columns = a_v_b_columns
    psi_p.columns = a_v_b_columns
    exp_p.columns = a_v_b_columns
    
    
    PSI.to_csv(f'tmp/ds_v_dge/{a_v_b}.psi.tsv.gz', sep='\t', index=False, header=True)
    delta_psi.to_csv(f'tmp/ds_v_dge/{a_v_b}.delta_psi.tsv.gz', sep='\t', index=False, header=True)
    logFC_psi.to_csv(f'tmp/ds_v_dge/{a_v_b}.logFC_psi.tsv.gz', sep='\t', index=False, header=True)
    logFC_exp.to_csv(f'tmp/ds_v_dge/{a_v_b}.logFC_exp.tsv.gz', sep='\t', index=False, header=True)
    cluster_counts.to_csv(f'tmp/ds_v_dge/{a_v_b}.cluster_counts.tsv.gz', sep='\t', index=False, header=True)
    cluster_fraction.to_csv(f'tmp/ds_v_dge/{a_v_b}.cluster_fraction.tsv.gz', sep='\t', index=False, header=True)
    psi_p.to_csv(f'tmp/ds_v_dge/{a_v_b}.psi_p.tsv.gz', sep='\t', index=False, header=True)
    exp_p.to_csv(f'tmp/ds_v_dge/{a_v_b}.exp_p.tsv.gz', sep='\t', index=False, header=True)
    
    print('success!')



def main(args):
    infile = args.inputFile
    rds_file = infile.split('/')[-1]
    process_comparison(rds_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect DS v DGE cluster information.')
    parser.add_argument('-i', dest='inputFile', help='Rds file containing DS and DGE results.')
    args = parser.parse_args()
    main(args)