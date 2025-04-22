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
        
#     dge = rds[1]
#     with (ro.default_converter + pandas2ri.converter).context():
#         pd_from_r_df = ro.conversion.get_conversion().rpy2py(dge)
#         dge = pd_from_r_df
        
    intron_list = []
    cluster_list = []
    tissue1_median_count_list = []
    tissue2_median_count_list = []
    
    tissue1, tissue2 = a_v_b.split('_v_')
    
    with gzip.open(
        f'/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds/GTEx/{a_v_b}/ds_perind_numers.counts.noise_by_intron.lf1.gz'
    ) as fh:
        header_lf1 = fh.readline().decode().rstrip().split(' ')#[1:]
        header_tissue1 = np.array([x.startswith(tissue1) for x in header_lf1])
        header_tissue2 = np.array([x.startswith(tissue2) for x in header_lf1])
        for line in fh:
            line = line.decode().rstrip().split(' ')
            intron_list.append(line[0])
            cluster_list.append(line[0].split(':')[0] + ':' + line[0].split(':')[-1])
            
            counts = np.array([int(y) for y in line[1:]])
#             print(header_tissue2)
#             print(header_tissue1)
#             print(counts)
#             print(header_lf1)
#             print(len(counts))
#             print(len(header_tissue1))
#             print(len(header_tissue2))
            tissue1_counts = int(np.median(counts[header_tissue1]))
            tissue2_counts = int(np.median(counts[header_tissue2]))
            
            tissue1_median_count_list.append(tissue1_counts)
            tissue2_median_count_list.append(tissue2_counts)
    
    intron_median_counts = pd.DataFrame()
    intron_median_counts['intron'] = intron_list
    intron_median_counts['cluster'] = cluster_list
    intron_median_counts[tissue1] = tissue1_median_count_list
    intron_median_counts[tissue2] = tissue2_median_count_list
    
    bed = output[['intron', 'cluster', 'itype', 'ctype', 'gene_name', 'gene_id']]
    
    median_counts_output = pd.merge(bed, intron_median_counts, right_on=['intron', 'cluster'], left_on=['intron', 'cluster'])
    median_counts_output.to_csv(f'tmp/ds_v_dge/{a_v_b}.cluster_counts_per_tissue.tsv.gz', sep='\t', index=False, header=True)
    
def main(args):
    infile = args.inputFile
    rds_file = infile.split('/')[-1]
    process_comparison(rds_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect DS v DGE cluster information.')
    parser.add_argument('-i', dest='inputFile', help='Rds file containing DS and DGE results.')
    args = parser.parse_args()
    main(args)