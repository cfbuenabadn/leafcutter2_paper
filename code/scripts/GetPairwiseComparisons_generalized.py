import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
import os
from rpy2.robjects import pandas2ri
import gzip
from scipy.stats import spearmanr, pearsonr
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
#from tqdm import tqdm

pairwise_comparisons = os.listdir('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge_confounder/rds_files/')

import numpy as np
from scipy.stats import rankdata


def process_comparison(rds_file):
    a_v_b = rds_file.split('.rds')[0]
    readRDS = ro.r['readRDS']
    rds = readRDS(f'/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge_confounder/rds_files/{rds_file}')
    
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
    
    
    ds_sig_p = ds.loc[(ds.itype=='PR') & (ds.ctype=='PR') & (ds['p.adjust'] <= 1e-2) & (ds.deltapsi.abs() >= 0.1) & (ds.cluster_counts >= 10) & (ds.cluster_fraction >= 0.1)]

    ds_sig_p = pd.merge(ds_sig_p, dge.loc[(dge.padj <= 1e-2) & (dge.log2FoldChange.abs() >= 1)], 
                      left_on='gene_id', right_on = 'gene_id').drop_duplicates()

    ds_sig_p = ds_sig_p.groupby('gene_id')[['logef', 'log2FoldChange']].median()
    
    
    
    Y = ds_sig_p['log2FoldChange']
    X = ds_sig_p['logef']
    X = sm.add_constant(X)

    model = sm.RLM(Y,X)
    results = model.fit()

    slope_p = results.params.loc['logef']
    const_p = results.params.loc['const']
    pval_p = results.pvalues.loc['logef']
    
    rlm_out_p = (slope_p, pval_p)
    spearman_rho_p = spearmanr(ds_sig_p.logef, ds_sig_p.log2FoldChange)
    pearson_rho_p = pearsonr(ds_sig_p.logef, ds_sig_p.log2FoldChange)



    ds_sig = ds.loc[(ds.itype=='UP') & (ds.ctype=='PR,UP') & (ds['p.adjust'] <= 1e-2) & (ds.deltapsi.abs() >= 0.1) & (ds.cluster_counts >= 10) & (ds.cluster_fraction >= 0.1)]

    ds_sig = pd.merge(ds_sig, dge.loc[(dge.padj <= 1e-2) & (dge.log2FoldChange.abs() >= 1)], 
                      left_on='gene_id', right_on = 'gene_id').drop_duplicates()

    ds_sig = ds_sig.groupby('gene_id')[['logef', 'log2FoldChange']].median()


    Y = ds_sig['log2FoldChange']
    X = ds_sig['logef']
    X = sm.add_constant(X)

    model = sm.RLM(Y,X)
    results = model.fit()

    slope = results.params.loc['logef']
    const = results.params.loc['const']
    pval = results.pvalues.loc['logef']
    
    rlm_out = (slope, pval)
    spearman_rho = spearmanr(ds_sig.logef, ds_sig.log2FoldChange)
    pearson_rho = pearsonr(ds_sig.logef, ds_sig.log2FoldChange)


#     fig, ax = plt.subplots(figsize=(3, 3))
#     plt.scatter(ds_sig.log2FoldChange, ds_sig.logef, alpha=0.2)
#     rho = spearmanr(ds_sig.logef, ds_sig.log2FoldChange)[0]
#     x = np.array([-5, 5])
#     y = const + (x*(rho))
#     plt.plot(x, y, c='red', linewidth=3)
#     plt.xlabel('logFC PSI')
#     plt.ylabel('logFC expression')
#     plt.title(a_v_b)

#     plt.show()
    
#     print(pearson_rho)
    
    return a_v_b, rlm_out, spearman_rho, pearson_rho, ds_sig.shape[0], rlm_out_p, spearman_rho_p, pearson_rho_p, ds_sig_p.shape[0]

a_v_b_list = []
slope_list = []
spearman_list = []
pearson_list = []
n_list = []

slope_list_p = []
spearman_list_p = []
pearson_list_p = []
n_list_p = []
for rds_file in pairwise_comparisons:
    print(rds_file)
    try:
        a_v_b, slope, spear, pear, n, slope_p, spear_p, pear_p, n_p = process_comparison(rds_file)
        slope_list.append(slope)
        a_v_b_list.append(a_v_b)
        spearman_list.append(spear)
        pearson_list.append(pear)
        n_list.append(n)
        
        slope_list_p.append(slope_p)
        spearman_list_p.append(spear_p)
        pearson_list_p.append(pear_p)
        n_list_p.append(n_p)
    except:
        continue
        
summary_df = pd.DataFrame()
summary_df['tissue_comparison'] = a_v_b_list
summary_df['slope'] = [x[0] for x in slope_list]
summary_df['slope_pval'] = [x[1] for x in slope_list]
summary_df['spearman'] = [x[0] for x in spearman_list]
summary_df['spearman_pval'] = [x[1] for x in spearman_list]
summary_df['pearson'] = [x[0] for x in pearson_list]
summary_df['pearson_pval'] = [x[1] for x in pearson_list]
summary_df['n'] = n_list


summary_df['slope_p'] = [x[0] for x in slope_list_p]
summary_df['slope_pval_p'] = [x[1] for x in slope_list_p]
summary_df['spearman_p'] = [x[0] for x in spearman_list_p]
summary_df['spearman_pval_p'] = [x[1] for x in spearman_list_p]
summary_df['pearson_p'] = [x[0] for x in pearson_list_p]
summary_df['pearson_pval_p'] = [x[1] for x in pearson_list_p]
summary_df['n_p'] = n_list_p

summary_df.to_csv('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge_confounder/pairwise_comparisons.tab.gz', sep='\t', header=True, index=False)

