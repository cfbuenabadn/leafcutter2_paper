import pandas as pd
import sys

gene_names = pd.read_csv('annotations/gencode.v34.gene_id_name.tab.gz', sep='\t', names=['gene_id', 'gene_name'])

arguments = sys.argv
condition = arguments[1]

counts_output = f'results/dge/HeLa/{condition}_v_controls_counts.tsv'
coldata = f'results/dge/HeLa/{condition}_v_controls_coldata.tsv'

cnts_file_1 = 'resources/HeLa/counts/HeLa_controls/Counts.txt'
cnts_file_2 = f'resources/HeLa/counts/{condition}/Counts.txt'



with open(cnts_file_1, 'r') as fh:
    with open(cnts_file_2, 'r') as fh2:
        fh.readline()
        fh2.readline()
        
        fh_cols = fh.readline().rstrip().split()
        fh2_cols = fh2.readline().rstrip().split()
        
        samples_1 = [x.split('/')[-3] for x in fh_cols[6:]]
        samples_2 = [x.split('/')[-3] for x in fh2_cols[6:]]
        
        with open(coldata, 'w') as cd_fh:
            line = 'sample_id\ttissue\n'
            cd_fh.write(line)
            for sample in samples_1:
                line = sample + '\tHeLa_controls\n'
                cd_fh.write(line)
            
            for sample in samples_2:
                line = sample + f'\t{condition}\n'
                cd_fh.write(line)
        
        with open(counts_output, 'w') as fh_out:
            first_line = 'Name\tDescription\t' + '\t'.join(samples_1 + samples_2) + '\n'
            fh_out.write(first_line)
            for line_1 in fh:
                line_2 = fh2.readline()
                row_1 = line_1.rstrip().split('\t')
                row_2 = line_2.rstrip().split('\t')

                gene = row_1[0]

                gene_name = gene_names.loc[gene_names.gene_id == gene].gene_name.iloc[0]

                cnts_1 = list(row_1[6:])
                cnts_2 = list(row_2[6:])
                
                line_out = [gene, gene_name] + cnts_1 + cnts_2
                line_out = '\t'.join(line_out) + '\n'
                
                fh_out.write(line_out)