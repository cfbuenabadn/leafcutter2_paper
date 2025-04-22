import pandas as pd
import numpy as np
import sys

arguments = sys.argv
input_groups = arguments[1]
gtex_upf3a = arguments[2]
output = arguments[3]

groups = pd.read_csv(input_groups, sep=' ', 
                  names = ['sample', 'tissue'])

gtex_upf3a = pd.read_csv(gtex_upf3a, sep='\t')

out_groups = groups.merge(gtex_upf3a[['sample', 'logUPF3A']], left_on='sample', right_on='sample')


out_data.to_csv(output, sep=' ', header=False, index=False)