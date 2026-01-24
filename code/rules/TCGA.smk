rule GetTCGAJuncFiles:
    output:
        'resources/TCGA/{chrom}/Breast_Cancer.done'
    resources:
        mem_mb = 24000
    log:
        'logs/tcga.{chrom}.log'
    wildcard_constraints:
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)])
    shell:
        """
        (python scripts/get_juncs_chr.py {wildcards.chrom}) &> {log}
        """

rule CollectTCGAJuncs:
    input:
        expand('resources/TCGA/{chrom}/Breast_Cancer.done', chrom = ['chr' + str(x) for x in range(1, 23)])

rule CollectTCGAJuncsForLeafCutter:
    input:
        expand('resources/TCGA/{chrom}/Breast_Cancer.{{sampleID}}.junc', chrom = ['chr' + str(x) for x in range(1, 23)])
    output:
        'resources/TCGA/JuncFiles/Breast_Cancer.{sampleID}.junc'
    wildcard_constraints:
        sampleID = '|'.join(tcga_samples)
    resources:
        mem_mb = 24000
    log:
        'logs/tcga_collect.{sampleID}.log'
    shell:
        """
        cat {input} > {output}
        """

rule CollectTCGAJuncsAll:
    input:
        expand('resources/TCGA/JuncFiles/Breast_Cancer.{sampleID}.junc', sampleID = tcga_samples)








rule TCGALeafCutter2:
    output:
        ds_numers = 'results/tcga/ds/leafcutter2.junction_counts.gz',
    params:
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        other_params = '-k --keepleafcutter1 ', # not keeping constitutive introns
    log:
        "logs/tcga_lc2.log"
    resources:
        mem_mb = 52000
    shell:
        """
        (python /project/yangili1/cfbuenabadn/leafcutter2/scripts/leafcutter2.py \
            -j config/tcga_junc_files.txt \
            -r results/tcga/ds/ \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} {params.other_params}) &> {log}
        """


rule RunLeafcutterDiffSplicingTCGA:
    message: """Run differential splicing analysis using leafcutter1's leafcutter_ds.R script"""
    input:
        ds_numers_lf1 = 'results/tcga/ds/leafcutter1_files/leafcutter2_perind_numers.counts.filtered.gz',
        ds_sample_group = 'results/tcga/ds/ds_sample_group.filtered.txt'
    output:
        "results/tcga/ds/leafcutter_ds/ds_cluster_significance.txt"
    params:
        Rscript = 'submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/tcga/ds/leafcutter_ds/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 5,
        MIN_SAMPLES_PER_GROUP = 3,
        MIN_COVERAGE = 5
    resources: cpu = 4, mem_mb = 58000, time = 2100
    threads: 4
    log: 'logs/RunLeafcutterDiffSplicingGtex/tcga.log'
    shell:
        '''
        /software/R-4.1.0-el7-x86_64/bin/Rscript {params.Rscript} --num_threads {threads} \
            --output_prefix {params.outprefix} \
            --min_samples_per_intron={params.MIN_SAMPLES_PER_INTRON} \
            --min_samples_per_group={params.MIN_SAMPLES_PER_GROUP} \
            --min_coverage={params.MIN_COVERAGE} \
            {input.ds_numers_lf1} {input.ds_sample_group} &> {log}
        '''
#'results/tcga/ds/ds_sample_group.txt'

