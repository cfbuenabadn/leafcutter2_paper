import gzip

import sys

arguments = sys.argv
input_vcf = arguments[1]
output_vcf = arguments[2]


with open(output_vcf, 'w') as vcf:

    with gzip.open(input_vcf, 'rb') as vcf2:
        for line in vcf2:
            linea = line.decode()
            if linea.startswith('#CHR'):
                linea = linea.replace('0_', '')
            if not linea.startswith('#'):
                linea = 'chr' + linea.replace("/", "|")
                
            vcf.write(linea)