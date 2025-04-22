import sys

arguments = sys.argv
input_bed = arguments[1]
output = arguments[2]

import gzip
with gzip.open(input_bed, 'rb') as fh:
    with gzip.open(output, 'wb') as fh2:
        x = fh.readline()#.decode()
        fh2.write(x)
        for line in fh:
            row = line.decode().rstrip().split('\t')
            bed = row[:6]
            y = [str(float(z)*100) for z in row[6:]]
            out_line = '\t'.join(bed + y) + '\n'
            fh2.write(out_line.encode())