N_PermutationChunks = 2
ChunkNumbers = range(0, 1+N_PermutationChunks) 

# top sQTL's nominal pass stats, used for beta-beta plot
rule NominaleQTL:
    message: '### Map QTL using nominal pass'
    input: 
        cov = '/project/yangili1/cfbuenabadn/SpliFi/code/results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.pca',
        bed = '/project/yangili1/cfbuenabadn/SpliFi/code/results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.bed.gz',
        vcf = '/project/yangili1/cfbuenabadn/SpliFi/code/results/geno/{datasource}/{group}/{chrom}.vcf.gz',
    output: temp('results/eqtl/{datasource}/{group}/cis_{window}/nom/chunks/{chrom}.{QTLTools_chunk_n}.txt')
    log: 'logs/NominaleQTL_{datasource}_{group}_{window}_{chrom}.{QTLTools_chunk_n}.log'
    params:
        NChunks = N_PermutationChunks
    resources: mem_mb = 54000
    wildcard_constraints:
        datasource = 'GTEx',
        group = '|'.join(tissues),
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)]),
        QTLTools_chunk_n = '|'.join([str(x) for x in ChunkNumbers])
    shell:
        '''
        module unload gsl && module load gsl/2.5
        {config[QTLtools]} cis \
            --std-err \
            --seed 123 \
            --nominal 1 \
            --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} \
            --vcf {input.vcf} --bed {input.bed} --cov {input.cov}  --out {output} \
            --window {wildcards.window} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        '''
        
        
rule Gather_QTLtools_cis_eqtls_pass:
    input:
        expand('results/eqtl/{{datasource}}/{{group}}/cis_{{window}}/nom/chunks/{{chrom}}.{QTLTools_chunk_n}.txt', QTLTools_chunk_n=ChunkNumbers)
    output:
       temp('results/eqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.temp.txt.gz'),
    log:
        "logs/Gather_QTLtools_cis_pass/eqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.log"
    wildcard_constraints:
        datasource = 'GTEx',
        group = '|'.join(tissues),
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)])
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """

rule tabixNominalPass_eQTLResults:
    """
    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
    for easy access with tabix
    """
    input:
        'results/eqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.temp.txt.gz',
    params:
        sort_temp = '-T ' + config['scratch'][:-1]
    wildcard_constraints:
        datasource = 'GTEx',
        group = '|'.join(tissues),
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)])
    output:
        txt = 'results/eqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.txt.gz',
        tbi = 'results/eqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.txt.gz.tbi',
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "logs/tabixNominalPass_eQTLResults/{datasource}/{group}/cis_{window}/nom/{chrom}.log"
    shadow: "shallow"
    shell:
        """
        (cat <(zcat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(zcat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k9,9 -k10,10n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
        tabix -b 10 -e10 -s9 {output.txt} &>> {log}
        """

rule collect_eQTLs:
    input:
        expand('results/eqtl/GTEx/{group}/cis_100000/nom/{chrom}.txt.gz', 
        group = tissues,
        chrom = ['chr' + str(x) for x in range(1, 23)])






