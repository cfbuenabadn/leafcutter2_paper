import gzip
import pandas as pd
import numpy as np
import re

gtex_table = '/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex_tables/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'

with gzip.open(gtex_table, 'rb') as fh:
    fh.readline()
    fh.readline()
    sample_list = fh.readline().decode().rstrip().split('\t')[2:]
    for line in fh:
        row = line.decode().rstrip().split('\t')
        if row[1] == 'UPF3A':
            upf3a = [float(x) for x in row[2:]]
        elif row[1] == 'UPF3B':
            upf3b = [float(x) for x in row[2:]]
        

gtex_samples = pd.read_csv('/project2/mstephens/cfbuenabadn/gtex-stm/data/sample.tsv', sep='\t')

gtex_sampleID_list = []
gtex_newID_list = []
for idx, sample_ in gtex_samples.iterrows():
    tissue = re.sub(' ', '', sample_.tissue_site_detail)
    sampleID = sample_['entity:sample_id']
    indID = '-'.join(sampleID.split('-')[:2])
    new_sampleID = tissue + '.' + indID + '.tsv.gz'
    gtex_sampleID_list.append(sampleID)
    gtex_newID_list.append(new_sampleID)
    
    
gtex_names = pd.DataFrame()
gtex_names['sampleID'] = gtex_sampleID_list
gtex_names['sample'] = gtex_newID_list

tpm_df = pd.DataFrame()
tpm_df['sampleID'] = sample_list
tpm_df['logUPF3A'] = np.log(np.array(upf3a))


gtex_upf3a = gtex_names.merge(tpm_df, left_on='sampleID', right_on='sampleID')

gtex_upf3a.to_csv('resources/GTEx/upf3a_logTPM.tsv.gz', sep='\t', header=True, index=False)
