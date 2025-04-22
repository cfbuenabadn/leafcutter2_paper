import pandas as pd
import numpy as np
import sys

arguments = sys.argv
input_groups = arguments[1]
input_counts = arguments[2]
output = arguments[3]

groups = pd.read_csv(input_groups, sep=' ', 
                  names = ['sample', 'tissue'])

data = pd.read_csv(input_counts, sep=' ')

total_up = data.loc[data.chrom.apply(lambda x: x.endswith(':UP')), data.columns[1:]].sum(axis=0)
total_juncs = data[data.columns[1:]].sum(axis=0)

datalog = pd.DataFrame(np.log(total_up/total_juncs))

# datalog = pd.DataFrame(np.log(data.loc[data.chrom.apply(lambda x: x.endswith(':UP')), data.columns[1:]].sum(axis=0)))
datalog.columns = ['logUP']

out_data = groups.merge(datalog, left_on='sample', right_index=True)

out_data.to_csv(output, sep=' ', header=False, index=False)