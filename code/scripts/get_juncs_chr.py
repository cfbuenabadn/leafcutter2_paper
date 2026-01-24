#!/bin/env python


def process_row(row):
    chrom = row.chromosome
    start = str(row.start)
    end = str(row.end)
    strand = row.strand
    name = '.'

    out_dict = {}

    samples = row.samples.lstrip(',').split(',')
    for sample in samples:
        sample_id, counts = sample.split(':')

        # print([chrom, start, end, name, counts, strand])
        out_row = '\t'.join([chrom, start, end, name, counts, strand]) + '\n'

        out_dict.update({sample_id:[out_row]})

    return out_dict

def process_gene(gene_df):
    gene_dict = {}
    for idx, row in gene_df.iterrows():
        row_dict = process_row(row)
        for sample_id in row_dict.keys():
            if sample_id in gene_dict.keys():
                gene_dict[sample_id].extend(row_dict[sample_id])
            else:
                gene_dict.update({sample_id:row_dict[sample_id]})
    return gene_dict

if __name__ == "__main__":
    
    import pandas as pd
    import numpy as np
    import os
    import sys

    arguments = sys.argv
    chrom = arguments[1]
    tcga_samples = pd.read_csv('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/samples.tsv', sep='\t')

    bc_primary_tumor = list(tcga_samples.loc[
    (tcga_samples.cgc_sample_sample_type == 'Primary Tumor') & (tcga_samples['gdc_cases.project.name'] == 'Breast Invasive Carcinoma'
                                                               ) & (tcga_samples.study == 'BRCA')].rail_id)

    bc_normal_tissue = list(tcga_samples.loc[
        (tcga_samples.cgc_sample_sample_type == 'Solid Tissue Normal') & (
            tcga_samples['gdc_cases.project.name'] == 'Breast Invasive Carcinoma'
                                                                   ) & (tcga_samples.study == 'BRCA')].rail_id)

    chrom_dir = f'/project2/yangili1/qhauck/splicing_PCA/pca_testing/tcga_sampling/all_genes/{chrom}/snaptron_output/'
    genes = [x.split('_')[0] for x in os.listdir(chrom_dir)]

    for gene in genes:
        try:
            print(chrom_dir + f'{gene}_snaptron.tsv')
            gene_snap = pd.read_csv(chrom_dir + f'{gene}_snaptron.tsv', sep='\t')
            gene_dict = process_gene(gene_snap)
            bc_normal_samples = [x for x in gene_dict.keys() if int(x) in bc_normal_tissue]
            bc_tumor_samples = [x for x in gene_dict.keys() if int(x) in bc_primary_tumor]
    
            all_samples = bc_normal_samples + bc_tumor_samples
    
            for sample_id in all_samples:
                out_file = f'resources/TCGA/{chrom}/Breast_Cancer.{sample_id}.junc'
                with open(out_file, 'a') as fh:
                    fh.writelines(gene_dict[sample_id])
        except:
            continue


    with open(f'resources/TCGA/{chrom}/Breast_Cancer.done', 'w') as fh:
        fh.write('Done!')
    


