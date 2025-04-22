pairwise_ds_dge = []
with open('config/pairwise_comparisons.txt', 'r') as fh:
    for c in fh:
        pairwise_ds_dge.append(c.rstrip())
        
tissues = []
with open('config/tissues.txt', 'r') as fh:
    for tissue in fh:
        tissues.append(tissue.rstrip())

def much_more_mem_after_first_attempt(wildcards, attempt):
    if int(attempt) == 1:
        return 4000
    else:
        return 62000
    
import pandas as pd
gwas_df = pd.read_csv('config/gwas_traits.txt', sep='\t', names = ['trait', 'accession', 'source', 'assembly'])
gwas_traits = list(gwas_df.trait)
gwas_new_traits = list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait)

leadSNPs = pd.read_csv('resources/gwas/LeadSnpWindows.bed', sep='\t', names = ['chrom', 'start', 'end', 'gwas_loci'])
gwas_loci = list(leadSNPs.gwas_loci)


bigwig_tissues = ['Brain_Anterior_cingulate_cortex_BA24', 'Brain_Frontal_Cortex_BA9', 'Heart_Atrial_Appendage', 'Lung', 
                  'Skin_Not_Sun_Exposed_Suprapubic', 'Brain_Cortex', 'Brain_Putamen_basal_ganglia', 'Liver', 'Muscle_Skeletal',
                  'Whole_Blood']

ba24_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Brain_Anterior_cingulate_cortex_BA24/') if x.endswith('.gz')]

ba9_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Brain_Frontal_Cortex_BA9/') if x.endswith('.gz')]

bc_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Brain_Cortex/') if x.endswith('.gz')]

bputamen_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Brain_Putamen_basal_ganglia/') if x.endswith('.gz')]

heart_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Heart_Atrial_Appendage/') if x.endswith('.gz')]

lung_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Lung/') if x.endswith('.gz')]

skin_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Skin_Not_Sun_Exposed_Suprapubic/') if x.endswith('.gz')]

liver_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Liver/') if x.endswith('.gz')]

ms_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Muscle_Skeletal/') if x.endswith('.gz')]

wb_samples = [x.split('.')[0] for x in os.listdir('/project2/yangili1/GTEx_v8/bedGraph/Whole_Blood/') if x.endswith('.gz')]


hela_samples = pd.read_csv('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/config/HeLa_samples.tsv', sep='\t')