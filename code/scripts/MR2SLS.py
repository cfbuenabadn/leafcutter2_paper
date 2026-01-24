import sys
import pandas as pd
import numpy as np
import gzip
import pysam
import linearmodels
from matplotlib import pyplot as plt
import seaborn as sns
import os
from linearmodels.iv import IV2SLS
import statsmodels.api as sm
from tqdm import tqdm


arguments = sys.argv
tissue = arguments[1]
def get_intron_list(tissue):
    intron_list = []
    cluster_list = []
    median_count_list = []
    with gzip.open(
        f'/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/{tissue}/leafcutter_perind_numers.counts.noise_by_intron.gz'
    ) as fh:
        fh.readline()
        for line in fh:
            line = line.decode().rstrip().split(' ')
            intron_list.append(line[0])
            cluster_list.append(line[0].split(':')[0] + ':' + line[0].split(':')[-2])
            median_counts = int(np.median([int(y) for y in line[1:]]))
            median_count_list.append(median_counts)

    intron_median_counts = pd.DataFrame()
    intron_median_counts['phenotype_id'] = intron_list
    intron_median_counts['cluster'] = cluster_list
    intron_median_counts['median_counts'] = median_count_list
    intron_median_counts['intron'] = intron_median_counts.phenotype_id.apply(lambda x: ':'.join(x.split(':')[:-1]))
    
    return intron_median_counts
    
def get_perm_counts(tissue):
#     results/sqtl/GTEx/Ovary/cis_100000/perm/PermutationPass.Qval.txt.gz
    perm = pd.read_csv(f'../code/results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz', sep='\t')
    intron_median_counts = get_intron_list(tissue)
    perm_counts = pd.merge(perm, intron_median_counts, right_on=['phenotype_id', 'intron', 'cluster'], 
         left_on=['phe_id', 'intron', 'cluster'])
    counts_per_cluster = perm_counts.groupby('cluster').median_counts.sum().reset_index()
    counts_per_cluster.columns = ['cluster', 'cluster_counts']

    counts_per_gene = perm_counts.groupby('gene_id').median_counts.sum().reset_index()
    counts_per_gene.columns = ['gene_id', 'gene_counts']
    
    perm_counts = pd.merge(perm_counts, counts_per_cluster, left_on='cluster', right_on='cluster')
    perm_counts = pd.merge(perm_counts, counts_per_gene, left_on='gene_id', right_on='gene_id')
    perm_counts['intron_psi'] = perm_counts.median_counts/(perm_counts.cluster_counts+1e-10)
    perm_counts['cluster_psi'] = perm_counts.cluster_counts/(perm_counts.gene_counts+1e-10)
    UP_per_cluster = perm_counts.groupby('cluster').itype.apply(lambda x: (x == 'UP').sum()).reset_index()
    UP_per_cluster.columns = ['cluster', 'UP_juncs']
    perm_counts = pd.merge(perm_counts, UP_per_cluster, left_on='cluster', right_on='cluster')
    
    return perm_counts


def get_expression_effect(tissue, perm_qtl):
    var_id = perm_qtl.var_id
    chrom = perm_qtl.var_chr
    start = str(int(perm_qtl.var_from) - 1)
    
    end = str(int(perm_qtl.var_to) + 1)
    
    coords = f'{chrom}:{start}:{end}'
    
    gene = perm_qtl.gene_id.split('.')[0]
            
    nom_file = f'/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/eqtl/GTEx/{tissue}/cis_100000/nom/{chrom}.txt.gz'
    nom_df = run_tabix_on_nom(nom_file, coords, gene=gene)

    # nom_df.set_index('#phe_id').T
    # return nom_df
    nom_df = nom_df.loc[nom_df['var_id'] == var_id]

    return nom_df



def run_tabix_on_nom(nom_file, coords, gene):
    chrom, start, end = coords.split(':')
    start = int(start)
    end = int(end)

    tb = pysam.TabixFile(nom_file)
    
    columns = tb.header[0].rstrip().split('\t')

    fields_dict = {'#phe_id':columns[1:]}
    for line in tb.fetch(chrom, start, end):
        fields = line.strip().split('\t')
        gene_name = fields[0].split('.')[0]
        # print(fields)
        if gene_name == gene:
            fields_dict.update({gene_name:fields[1:]})

    fields_df = pd.DataFrame(fields_dict).set_index('#phe_id').T
    return fields_df

def fetch_phenotypes(stbx, etbx, intron, gene):
    gene = gene.split('.')[0]
    chrom, start, end, clu, annot = intron.split(':')
    start = int(start)
    end = int(end)

    intron_found = False
    gene_found = False
    
    for line in stbx.fetch(chrom, start, end):
        fields = line.strip().split('\t')
        field_intron = fields[3]
        if field_intron == intron:
            pheno_intron = [float(x) for x in fields[6:]]
            intron_found = True

    for line in etbx.fetch(chrom, start, end):
        fields = line.strip().split('\t')
        field_gene = fields[3].split('.')[0]
        
        if field_gene == gene:
            pheno_gene = [float(x) for x in fields[6:]]
            gene_found = True

    if gene_found and intron_found:
        return pheno_intron, pheno_gene
    else:
        return None, None


def fetch_variant(vcf_in, var_id):
    chrom, pos, ref, alt, tag = var_id.split('_') #"chr5_52897294_A_G_b38"
    pos = int(pos)
    for record in vcf_in.fetch(chrom, pos-1, pos+1):
        if (record.chrom == chrom) and (pos == record.pos) and (ref == record.ref) and (alt in record.alts):
            return record
    return None

def fetch_genotypes(vcf_in, var_id, samples):
    record = fetch_variant(vcf_in, var_id)
    if record is None:
        return None

    genotypes_list = []
    record_dict = dict(record.samples.items())
    for sample in samples:
        if sample in record_dict.keys():
            sample_gt = np.sum(record_dict[sample]['GT'])
        else:
            sample_gt = np.nan
        genotypes_list.append(sample_gt)

    return genotypes_list

def get_iv_input(stbx, etbx, intron, gene, vcf_in, var_id):
    assert stbx.header[0] == etbx.header[0]
    samples_list = etbx.header[0].rstrip().split('\t')[6:]

    pheno_intron, pheno_gene = fetch_phenotypes(stbx, etbx, intron, gene)
    genotype = fetch_genotypes(vcf_in, var_id, samples_list)

    if (pheno_intron is not None) and (pheno_gene is not None) and (genotype is not None):
        out_df = pd.DataFrame({'samples':samples_list, 'genotype':genotype, 'pheno_intron':pheno_intron, 'pheno_gene':pheno_gene})
        out_df = out_df.set_index('samples')
        return out_df

    else:
        return None


# perm_select

def process_perm_row(row, tissue, vcf_in):
    chrom = row.phe_chr
    intron = row.phe_id
    gene = row.gene_id.split('.')[0]
    var_id = row.var_id

    res_dir = '/project/yangili1/cfbuenabadn/SpliFi/code/results'

    sPheno = f'{res_dir}/pheno/noisy/GTEx/{tissue}/separateNoise/leafcutter.qqnorm_{chrom}.gz'
    ePheno = f'{res_dir}/eqtl/GTEx/{tissue}/qqnorm.sorted.{chrom}.bed.gz'
    etbx = pysam.TabixFile(ePheno)
    stbx = pysam.TabixFile(sPheno)

    pca = pd.read_csv(f'/project/yangili1/cfbuenabadn/SpliFi/code/results/eqtl/GTEx/{tissue}/qqnorm.sorted.bed.pca', sep=' ', index_col = 0).T

    iv_df = get_iv_input(stbx, etbx, intron, gene, vcf_in, var_id)

    merged_iv = iv_df.merge(pca, left_index = True, right_index=True)

    first_stage_F, intron_gene_beta, intron_gene_pval = get_2SLS_results(merged_iv)
    pval_1, pval_2, pval_3, beta_step1, beta_step3 = get_baron_kenny_results(merged_iv)

    abs_beta_diff = np.abs(beta_step3) - np.abs(beta_step1)
    beta_diff = beta_step3 - beta_step1

    out_row = [intron, gene, var_id, first_stage_F, intron_gene_beta, intron_gene_pval, pval_1, pval_2, pval_3, beta_step1, beta_step3, beta_diff, abs_beta_diff]
    return out_row


def get_tissue_mediation_stats(tissue, vcf_in):
    print(f'loading sQTLs for {tissue}')
    perm_counts = get_perm_counts(tissue)
    ctype='PR,UP'
    itype='UP'

    qmax = 1e-1
    clu_psi_min=0.1
    perm_select = perm_counts.loc[(perm_counts.ctype == ctype) & (perm_counts.itype == itype) & (perm_counts.q <= qmax) & (perm_counts.cluster_psi >= clu_psi_min)].copy()

    print('Working on UP introns...')
    
    up_mediation_stats = []
    i = 0
    for idx, row in perm_select.iterrows():
        try:
            row_stats = process_perm_row(row, tissue, vcf_in)
            up_mediation_stats.append(row_stats)
        except:
            continue

    up_mediation_stats = pd.DataFrame(up_mediation_stats, 
             columns = ['intron', 'gene', 'var', 'F', 'intron_gene_beta', 'itron_gene_pval', 'pval1', 
                        'pval2', 'pval3', 'beta_1', 
                        'beta_3', 'beta_diff', 'abs_beta_diff'])


    perm_select['gene'] = perm_select.gene_id.apply(lambda x: x.split('.')[0])

    UP_out = perm_select[['gene', 'phe_id', 'q', 'cluster_psi']].merge(up_mediation_stats,
                            left_on=['gene', 'phe_id'],  right_on = ['gene', 'intron'])

    UP_out.to_csv(f'results/mr_2sls/{tissue}_UP.tab.gz', sep='\t', header=True, index = False)
    
    #########
    print('Working on PR introns...')
    
    ctype='PR'
    itype='PR'

    qmax = 1e-1
    clu_psi_min=0.1
    perm_select = perm_counts.loc[(perm_counts.ctype == ctype) & (perm_counts.itype == itype) & (perm_counts.q <= qmax) & (perm_counts.cluster_psi >= clu_psi_min)]

    pr_mediation_stats = []
    i = 0
    for idx, row in tqdm(perm_select.iterrows(), leave=True, position=0):
        try:
            row_stats = process_perm_row(row, tissue, vcf_in)
            pr_mediation_stats.append(row_stats)
        except:
            continue

    pr_mediation_stats = pd.DataFrame(pr_mediation_stats, 
             columns = ['intron', 'gene', 'var', 'F', 'intron_gene_beta', 'intron_gene_pval', 'pval1', 
                        'pval2', 'pval3', 'beta_1', 
                        'beta_3', 'beta_diff', 'abs_beta_diff'])

    perm_select['gene'] = perm_select.gene_id.apply(lambda x: x.split('.')[0])

    PR_out = perm_select[['gene', 'phe_id', 'q', 'cluster_psi']].merge(pr_mediation_stats,
                            left_on=['gene', 'phe_id'],  right_on = ['gene', 'intron'])

    PR_out.to_csv(f'results/mr_2sls/{tissue}_PR.tab.gz', sep='\t', header=True, index = False)

    

    # out_dict = {'UP':up_mediation_stats, 'PR':pr_mediation_stats}
    # return out_dict


def baron_kenny_mediation(merged_iv):
    covariates = merged_iv.columns[3:]

    XG_cov   = sm.add_constant(pd.concat([merged_iv['genotype'], merged_iv[covariates]], axis=1))
    XG_I_cov  = sm.add_constant(pd.concat([merged_iv[['genotype', 'pheno_intron']], merged_iv[covariates]], axis=1))

    # Step 1: total genetic effect on expression
    step1 = sm.OLS(merged_iv['pheno_gene'], XG_cov).fit()

    # Step 2: UP splicing path
    step2 = sm.OLS(merged_iv['pheno_intron'], XG_cov).fit()

    # Step 3: mediated analysis
    step3 = sm.OLS(merged_iv['pheno_gene'], XG_I_cov).fit()

    return step1, step2, step3
    
    # print(model_c.summary())
    # print(model_a.summary())
    # print(model_bc.summary())



def MR2SLS(merged_iv):
    covariates = merged_iv.columns[3:]
    covs_formula = ' + '.join(covariates)
    formula = 'pheno_gene ~ 1 + ' + covs_formula +' + [pheno_intron ~ genotype]'


    # 2SLS with covariates
    model = IV2SLS.from_formula(
        formula,
        data=merged_iv
    )

    results = model.fit(cov_type='robust')
    return (results)


def get_2SLS_results(merged_iv):
    results = MR2SLS(merged_iv)
    first_stage_F = results.first_stage.diagnostics.loc['pheno_intron', 'f.stat']
    intron_gene_beta = results.params['pheno_intron']
    intron_gene_pval = results.pvalues['pheno_intron']

    return first_stage_F, intron_gene_beta, intron_gene_pval


def get_baron_kenny_results(merged_iv):
    step1, step2, step3 = baron_kenny_mediation(merged_iv)
    pval_1 = step1.pvalues.loc['genotype']
    pval_2 = step2.pvalues.loc['genotype']
    pval_3 = step3.pvalues.loc['genotype']

    beta_step1 = step1.params.loc['genotype']
    beta_step3 = step3.params.loc['genotype']

    return pval_1, pval_2, pval_3, beta_step1, beta_step3


if __name__ == '__main__':
    vcf_in = pysam.VariantFile('/project2/yangili1/cfbuenabadn/torino_paper/code/resources/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz')  # automatically detects format
    
    get_tissue_mediation_stats(tissue, vcf_in)
