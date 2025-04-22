'''
Differential splicing analysis &
Differential expression analysis &
related analysis & data preps & plots
'''

N_DIFFER = 80 # number of samples to choose for differnetial splicing and expression analysis


# ----------------------------------------------------------------------------------------
#           Cluster introns using all GTEx tissues at once
#         Then use the clusters to run leafcutter2 on specific tissues
#         This enables differential splicing analysis across tissues
# ----------------------------------------------------------------------------------------

# NOTE: Copied files from Chao, but changed names
rule LeafcutterForDSGtex:
    message: '### Run leafcutter2 on GTEx samples for differential splicing analysis'
    input:
        pre_clusters = '/project/yangili1/cfbuenabadn/SpliFi/code/results/ds/GTEx/all49tissues_refined_noisy',
        tissue1_flag = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/juncs/all49tissues/{ds_tissue_1}.done',
        tissue2_flag = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/juncs/all49tissues/{ds_tissue_2}.done',
    output:
        ds_numers = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        # ds_counts_lf1 is necessary for leafcutter_ds.R script
        ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    params:
        run_dir = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}',
        out_prefix = 'ds',
        tissue1_juncs = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/juncs/all49tissues/{ds_tissue_1}',
        tissue2_juncs = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/juncs/all49tissues/{ds_tissue_2}',
        NSamples = N_DIFFER, # select this number of samples per tissue type
        pre_clustered = '-c /project/yangili1/cfbuenabadn/SpliFi/code/results/ds/GTEx/all49tissues_refined_noisy',
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        other_params = '-k ', # not keeping constitutive introns
        py_script  = 'submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        py_script2 = 'scripts/makeSampleGroupFileForDifferentialSplicing.py'
    log: 'logs/LeafcutterForDSGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    resources: cpu = 1, time = 2100, mem_mb = 25000
    shell:
        '''
        # run leafcutter2 for differential analysis
        python {params.py_script} \
            -j <(cat <(ls {params.tissue1_juncs}*tsv.gz | head -{params.NSamples}) <(ls {params.tissue2_juncs}*tsv.gz | head -{params.NSamples})) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} \
            {params.pre_clustered} {params.other_params} &> {log}

        # make sample group file for differential splicing analysis
        python {params.py_script2} -i {output.ds_numers} -o {output.ds_numers_lf1} -s {output.ds_sample_group} &>> {log}

        '''


rule RunLeafcutterDiffSplicingGtex:
    message: """Run differential splicing analysis using leafcutter1's leafcutter_ds.R script"""
    input:
        ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    output:
        "results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_cluster_significance.txt"
        #flag = touch('results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/done')
        # produces two files:
        # 1. {outprefix}_effect_sizes.txt
        # 2. {outprefix}_manual_ds_cluster_significance.txt
    params:
        Rscript = 'submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 5,
        MIN_SAMPLES_PER_GROUP = 3,
        MIN_COVERAGE = 5
    resources: cpu = 4, mem_mb = 58000, time = 2100
    threads: 4
    log: 'logs/RunLeafcutterDiffSplicingGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        /software/R-4.1.0-el7-x86_64/bin/Rscript {params.Rscript} --num_threads {threads} \
            --output_prefix {params.outprefix} \
            --min_samples_per_intron={params.MIN_SAMPLES_PER_INTRON} \
            --min_samples_per_group={params.MIN_SAMPLES_PER_GROUP} \
            --min_coverage={params.MIN_COVERAGE} \
            {input.ds_numers_lf1} {input.ds_sample_group} &> {log}


        '''
        
def GetHeLaInput(wildcards):
    return hela_samples.loc[(hela_samples.sampleID == wildcards.sampleID) & (hela_samples.condition==wildcards.condition)].path_to_file
        
rule PrepareHeLaKDForLeafcutter:
    input:
        GetHeLaInput
    output:
        'resources/HeLa/juncs/{condition}.{sampleID}.tsv.gz'
    log:
        'logs/HeLa/transform_STAR_to_Leafcutter_input.{condition}.{sampleID}.log'
    shell:
        """
        python scripts/STAR2Junc.py {input} {output} &> {log}
        """
        
def GetHeLaConditionInput(wildcards):
    samples = list(hela_samples.loc[hela_samples.condition==wildcards.condition].sampleID)
    return expand('resources/HeLa/juncs/{condition}.{sampleID}.tsv.gz', condition = [wildcards.condition], sampleID=samples)
        
rule HeLaDoneFlag:
    input:
        GetHeLaConditionInput
    output:
        'resources/HeLa/juncs/{condition}.done.txt'
    shell:
        """
        touch {output}
        """
        
rule LeafcutterForDSHeLa:
    message: '### Run leafcutter2 on GTEx samples for differential splicing analysis'
    input:
        pre_clusters = '/project/yangili1/cfbuenabadn/SpliFi/code/results/ds/GTEx/all49tissues_refined_noisy',
        cond_flag = 'resources/HeLa/juncs/{condition}.done.txt',
        control_flag = 'resources/HeLa/juncs/HeLa_controls.done.txt'
    output:
        ds_numers = 'results/ds/HeLa/{condition}_v_controls/ds_perind_numers.counts.noise_by_intron.gz',
        ds_numers_lf1 = 'results/ds/HeLa/{condition}_v_controls/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/HeLa/{condition}_v_controls/ds_sample_group.txt'
    wildcard_constraints:
        condition = 'HeLa_SMG6|HeLa_SMG7|HeLa_UPF1|HeLa_dKD'
    params:
        run_dir = 'results/ds/HeLa/{condition}_v_controls',
        out_prefix = 'ds',
        tissue1_juncs = 'resources/HeLa/juncs/{condition}.',
        tissue2_juncs = 'resources/HeLa/juncs/HeLa_controls.',
        pre_clustered = '-c /project/yangili1/cfbuenabadn/SpliFi/code/results/ds/GTEx/all49tissues_refined_noisy',
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        other_params = '-k ', # not keeping constitutive introns
        py_script  = 'submodules/leafcutter2/scripts/leafcutter2_regtools.py',
        py_script2 = 'scripts/makeSampleGroupFileForDifferentialSplicing.py'
    log: 'logs/LeafcutterForDSGtex/{condition}_v_controls.log'
    resources: cpu = 1, time = 2100, mem_mb = 25000
    shell:
        '''
        # run leafcutter2 for differential analysis
        python {params.py_script} \
            -j <(cat <(ls {params.tissue1_juncs}*tsv.gz) <(ls {params.tissue2_juncs}*tsv.gz)) \
            -r {params.run_dir} \
            -o {params.out_prefix} \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} \
            {params.pre_clustered} {params.other_params} &> {log}

        # make sample group file for differential splicing analysis
        python {params.py_script2} -i {output.ds_numers} -o {output.ds_numers_lf1} -s {output.ds_sample_group} &>> {log}

        '''

rule collect_HeLa_input:
    input:
        expand('results/ds/HeLa/{condition}_v_controls/ds_sample_group.txt', 
            condition = ['HeLa_SMG6', 'HeLa_SMG7', 'HeLa_UPF1', 'HeLa_dKD'])
            
            
rule RunLeafcutterDiffSplicingHeLa:
    message: """Run differential splicing analysis using leafcutter1's leafcutter_ds.R script"""
    input:
        ds_numers_lf1 = 'results/ds/HeLa/{condition}_v_controls/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/HeLa/{condition}_v_controls/ds_sample_group.txt'
    output:
        "results/ds/HeLa/{condition}_v_controls/ds_cluster_significance.txt"
    params:
        Rscript = 'submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/HeLa/{condition}_v_controls/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 2,
        MIN_SAMPLES_PER_GROUP = 2,
        MIN_COVERAGE = 5
    resources: cpu = 4, mem_mb = 58000, time = 2100
    threads: 4
    log: 'logs/RunLeafcutterDiffSplicingGtex/{condition}_v_controls.log'
    wildcard_constraints:
        condition = 'HeLa_SMG6|HeLa_SMG7|HeLa_UPF1|HeLa_dKD'
    shell:
        '''
        /software/R-4.1.0-el7-x86_64/bin/Rscript {params.Rscript} --num_threads {threads} \
            --output_prefix {params.outprefix} \
            --min_samples_per_intron={params.MIN_SAMPLES_PER_INTRON} \
            --min_samples_per_group={params.MIN_SAMPLES_PER_GROUP} \
            --min_coverage={params.MIN_COVERAGE} \
            {input.ds_numers_lf1} {input.ds_sample_group} &> {log}


        '''

# NOTE to make leafcutter's leafcutter_ds.R script work, 
# must modify the _perind_numers.counts.noise_by_intron.gz file
# First row should not have the 'chrom' first column
# first columns can only be like `chr1:827775:829002:clu_1_+` because 
# the R function only expect to split by ":" into 4 columns.



rule MakeUPF3A_counfounder:
    output:
        'resources/GTEx/upf3a_logTPM.tsv.gz'
    log:
        "logs/leafcutter_DS_UPF3A.log"
    resources: mem_mb = 8000
    shell:
        """
        python scripts/prepare_upf3a_table.py &> {log}
        """
    


rule make_UPF3Acounfounder_groups:
    input:
        groups = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt',
        conf = 'resources/GTEx/upf3a_logTPM.tsv.gz'
    output:
        'results/ds/GTEx_UPF3A/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    log:
        "logs/leafcutter_DS_UPF3A/prepare_groups/{ds_tissue_2}_v_{ds_tissue_1}.log"
    resources: mem_mb = 12000
    shell:
        """
        python scripts/prepare_groups_with_upf3a_counfounder.py {input.groups} {input.counts} {output} &> {log}
        """
        
        


rule make_counfounder_groups:
    input:
        groups = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt',
        counts = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz'
    output:
        'results/ds/GTEx_confounder/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    log:
        "logs/leafcutter_DS_confounder/prepare_groups/{ds_tissue_2}_v_{ds_tissue_1}.log"
    resources: mem_mb = 12000
    shell:
        """
        python scripts/prepare_groups_with_counfounder.py {input.groups} {input.counts} {output} &> {log}
        """



rule RunLeafcutterDiffSplicingGtex_confounder:
    message: """Run differential splicing analysis using leafcutter1's leafcutter_ds.R script"""
    input:
        ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/GTEx_confounder/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    output:
        "results/ds/GTEx_confounder/{ds_tissue_2}_v_{ds_tissue_1}/ds_cluster_significance.txt"
    params:
        Rscript = 'submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/GTEx_confounder/{ds_tissue_2}_v_{ds_tissue_1}/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 5,
        MIN_SAMPLES_PER_GROUP = 3,
        MIN_COVERAGE = 5
    resources: cpu = 4, mem_mb = 58000, time = 2100
    threads: 4
    log: 'logs/RunLeafcutterDiffSplicingGtex_conf/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        /software/R-4.1.0-el7-x86_64/bin/Rscript {params.Rscript} --num_threads {threads} \
            --output_prefix {params.outprefix} \
            --min_samples_per_intron={params.MIN_SAMPLES_PER_INTRON} \
            --min_samples_per_group={params.MIN_SAMPLES_PER_GROUP} \
            --min_coverage={params.MIN_COVERAGE} \
            {input.ds_numers_lf1} {input.ds_sample_group} &> {log}


        '''


use rule RunLeafcutterDiffSplicingGtex_confounder as RunLeafcutterDiffSplicingGtex_UPF3A with:
    input:
        ds_numers_lf1 = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.lf1.gz',
        ds_sample_group = 'results/ds/GTEx_UPF3A/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt'
    output:
        "results/ds/GTEx_UPF3A/{ds_tissue_2}_v_{ds_tissue_1}/ds_cluster_significance.txt"
    params:
        Rscript = 'submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/GTEx_UPF3A/{ds_tissue_2}_v_{ds_tissue_1}/ds', # note you need to include path!



## -----------------------------------------------------------------------------
##   GTEx expression data
## -----------------------------------------------------------------------------

rule ExtractGTExGeneExpression:
    input: 
      tpm = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
      cnt = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
    output: 
      tpm = 'resources/GTEx/expression/{tissue}_gene_tpm.tsv.gz',
      cnt = 'resources/GTEx/expression/{tissue}_gene_reads.tsv.gz',
    params: 
        py_script = 'scripts/extract_gtex_gene_expression.py',
        junc_meta = config['Dataset']['GTEx']['Junc_meta']
    log: 'logs/ExtractGTExGeneExpression/{tissue}.log'
    resources: mem_mb = 24000
    shell:
        '''
        echo extracting tpm from {input.tpm} ...&> {log}
        python {params.py_script} \
            -I {input.tpm} \
            -M {params.junc_meta} \
            -O {output.tpm} \
            -T {wildcards.tissue} &>> {log}

        echo extracting raw counts from {input.cnt} ...
        python {params.py_script} \
            -I {input.cnt} \
            -M {params.junc_meta} \
            -O {output.cnt} \
            -T {wildcards.tissue} &>> {log}

        echo "Number of lines in {output.tpm} : $(zcat {output.tpm} | wc -l)" &>> {log}
        echo "Number of lines in {output.cnt} : $(zcat {output.cnt} | wc -l)" &>> {log}
        '''

# prepare GTEX dge data
# NOTE: Copied files from Chao, but changed names
# are not matched with ds samples! 
# Should use the ds_sample_group.txt file to get matching samples for each tissue
rule PrepareGTExDGE:
    input: 
        cnt1 = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/expression/{dge_tissue1}_gene_reads.tsv.gz',
        cnt2 = '/project/yangili1/cfbuenabadn/SpliFi/code/resources/GTEx/expression/{dge_tissue2}_gene_reads.tsv.gz',
        ds_samples = 'results/ds/GTEx/{dge_tissue2}_v_{dge_tissue1}/ds_sample_group.txt'
    output: 
        #NOTE: tissue2 is intended to be numerator, and tissue 1 denominator, in subsequent dge step
        cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
        coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv',
    params:
        R_script = 'scripts/prepare_GTEx_dge.R',
        outdir = 'results/dge/GTEx'
    log: 'logs/PrepareGTExDGE/{dge_tissue2}_v_{dge_tissue1}.log'
    resources: mem_mb = 24000
    shell:
        '''
        Rscript {params.R_script} {input.cnt1} {input.cnt2} {input.ds_samples} {params.outdir} &> {log}
        ls {output.cnt} {output.coldata} &>> {log}
        '''

# run dge
rule DgeGtex:
  input: 
    cnt = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_counts.tsv',
    coldata = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_coldata.tsv'
  output:
    dge = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}_dge_genes.tsv',
  params:
    R_script = 'scripts/dge.R',
    outprefix = 'results/dge/GTEx/{dge_tissue2}_v_{dge_tissue1}',
    min_reads = 10,
    min_samples = 10,
  log: 'logs/DgeGtex/{dge_tissue2}_v_{dge_tissue1}.log'
  resources: mem_mb = 24000
  shell:
    '''
    /software/R-4.1.0-el7-x86_64/bin/Rscript {params.R_script} \
            {input.cnt} {input.coldata} {params.outprefix} \
            {params.min_reads} {params.min_samples}  &> {log}
    ls {output.dge} &>> {log}

    '''








def GetBamForPhenotype(wildcards):
    return list(hela_samples.loc[hela_samples.condition==wildcards.condition].bam_file)


rule featureCounts:
    input:
        bam = GetBamForPhenotype,
        annotations = 'annotations/gencode.v34.primary_assembly.annotation.gtf'
    output:
        "resources/HeLa/counts/{condition}/Counts.txt"
    params:
        extraParams = "-s 2",
        paired = ""
    threads:
        8
    wildcard_constraints:
        Phenotype = "|".join(list(hela_samples.condition.unique()))
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{condition}.log"
    shell:
        """
        featureCounts {params.paired} {params.extraParams} -T {threads} --ignoreDup --primary -a {input.annotations} -o {output} {input.bam} &> {log}
        """


rule collectfeatureCounts:
    input:
        expand("resources/HeLa/counts/{condition}/Counts.txt", condition = list(hela_samples.condition.unique()))


rule PrepareHeLaCountsForDGE:
    input:
        "resources/HeLa/counts/HeLa_controls/Counts.txt",
        "resources/HeLa/counts/{condition}/Counts.txt"
    output:
        'results/dge/HeLa/{condition}_v_controls_counts.tsv',
        'results/dge/HeLa/{condition}_v_controls_coldata.tsv'
    log:
        "logs/prepare_dge_hela.{condition}.log"
    resources:
        mem = 8000
    wildcard_constraints:
        condition = 'HeLa_SMG6|HeLa_SMG7|HeLa_UPF1|HeLa_dKD'
    shell:
        """
        python scripts/prepare_HeLa_dge.py {wildcards.condition} &> {log}
        """
        
        
rule DgeHeLa:
  input: 
    cnt = 'results/dge/HeLa/{condition}_v_controls_counts.tsv',
    coldata = 'results/dge/HeLa/{condition}_v_controls_coldata.tsv'
  output:
    dge = 'results/dge/HeLa/{condition}_v_controls_dge_genes.tsv',
  params:
    R_script = 'scripts/dge.R',
    min_reads = 10,
    min_samples = 2,
  log: 'logs/DgeHeLa/{condition}_v_controls.log'
  resources: mem_mb = 24000
  shell:
    '''
    /software/R-4.1.0-el7-x86_64/bin/Rscript {params.R_script} \
            {input.cnt} {input.coldata} results/dge/HeLa/{wildcards.condition}_v_controls \
            {params.min_reads} {params.min_samples}  &> {log}
    ls {output.dge} &>> {log}

    '''

# -----------------------------------------------------------------------------
#   collect ds and dge results
# -----------------------------------------------------------------------------

rule Ds_Dge_Results:
  input: 
      'results/dge/GTEx/{a_v_b}_dge_genes.tsv',
      "results/ds/GTEx/{a_v_b}/ds_cluster_significance.txt"
  output: 'results/ds_v_dge/{a_v_b}.rds'
  params:
    R_script = 'scripts/prepDGE_DS_AnalysesData.R',
    contrast = lambda w: w.a_v_b,
    # must be v26 because it matches gene names from GTEx V8
    gtf = config['annotation']['csv']['v26_genes'],
    dsPrefix = 'results/ds/GTEx/',
    dgePrefix = 'results/dge/GTEx/',
    outPrefix = 'results/ds_v_dge/'
  resources: mem_mb = 24000
  log:
      'logs/collect_ds_v_dge.{a_v_b}.log'
  shell:
    '''
    (/software/R-4.1.0-el7-x86_64/bin/Rscript {params.R_script} {params.contrast} {params.gtf} {params.dsPrefix} {params.dgePrefix} {params.outPrefix}) &> {log}
    ls {output}
    '''

rule collect_ds_v_dge:
    input:
        expand('results/ds_v_dge/{a_v_b}.rds', a_v_b = pairwise_ds_dge)
        
rule collect_ds_v_dge_muscle_brain:
    input:
        "results/ds_v_dge/Brain-Cortex_v_Muscle-Skeletal.rds"
        
        
        
        
        
        
rule Ds_Dge_Results_HeLa:
  input: 
      'results/dge/HeLa/{condition}_v_controls_dge_genes.tsv',
      "results/ds/HeLa/{condition}_v_controls/ds_cluster_significance.txt"
  output: 'results/HeLa_ds_v_dge/{condition}_v_controls.rds'
  params:
    R_script = 'scripts/prepDGE_DS_AnalysesData.R',
    # must be v26 because it matches gene names from GTEx V8
    gtf = config['annotation']['csv']['v26_genes'],
    dsPrefix = 'results/ds/HeLa/',
    dgePrefix = 'results/dge/HeLa/',
    outPrefix = 'results/HeLa_ds_v_dge/'
  wildcard_constraints:
        condition = 'HeLa_SMG6|HeLa_SMG7|HeLa_UPF1|HeLa_dKD'
  resources: mem_mb = 24000
  log:
      'logs/collect_ds_v_dge.{condition}_v_controls.log'
  shell:
    '''
    (/software/R-4.1.0-el7-x86_64/bin/Rscript {params.R_script} {wildcards.condition}_v_controls {params.gtf} {params.dsPrefix} {params.dgePrefix} {params.outPrefix}) &> {log}
    ls {output}
    '''

rule collect_ds_v_dge_HeLa:
    input:
        expand('results/HeLa_ds_v_dge/{condition}_v_controls.rds', 
        condition = ['HeLa_SMG6', 'HeLa_SMG7', 'HeLa_UPF1', 'HeLa_dKD'])        
        
        
rule collect_ds_v_dge_confounder:
    input:
        expand('results/ds_v_dge_confounder/rds_files/{a_v_b}.rds', a_v_b = pairwise_ds_dge),
        expand('results/ds/GTEx/{a_v_b}/ds_perind_numers.counts.noise_by_intron.lf1.gz', a_v_b = pairwise_ds_dge),
        
rule collect_ds_v_dge_UPF3A:
    input:
        expand('results/ds_v_dge_UPF3A/{a_v_b}.rds', a_v_b = pairwise_ds_dge)
        

rule Ds_Dge_Results_confounder:
  input: 
      'results/dge/GTEx/{a_v_b}_dge_genes.tsv',
      "results/ds/GTEx_confounder/{a_v_b}/ds_cluster_significance.txt"
  output: 'results/ds_v_dge_confounder/rds_files/{a_v_b}.rds'
  params:
    R_script = 'scripts/prepDGE_DS_AnalysesData.R',
    contrast = lambda w: w.a_v_b,
    # must be v26 because it matches gene names from GTEx V8
    gtf = config['annotation']['csv']['v26_genes'],
    dsPrefix = 'results/ds/GTEx_confounder/',
    dgePrefix = 'results/dge/GTEx/',
    outPrefix = 'results/ds_v_dge_confounder/rds_files/'
  resources: mem_mb = 24000
  log:
      'logs/collect_ds_v_dge.{a_v_b}.log'
  shell:
    '''
    (/software/R-4.1.0-el7-x86_64/bin/Rscript {params.R_script} {params.contrast} {params.gtf} {params.dsPrefix} {params.dgePrefix} {params.outPrefix}) &> {log}
    ls {output}
    '''
    
rule Ds_Dge_Results_UPF3A:
  input: 
      'results/dge/GTEx/{a_v_b}_dge_genes.tsv',
      "results/ds/GTEx_UPF3A/{a_v_b}/ds_cluster_significance.txt"
  output: 'results/ds_v_dge_UPF3A/{a_v_b}.rds'
  params:
    R_script = 'scripts/prepDGE_DS_AnalysesData.R',
    contrast = lambda w: w.a_v_b,
    # must be v26 because it matches gene names from GTEx V8
    gtf = config['annotation']['csv']['v26_genes'],
    dsPrefix = 'results/ds/GTEx_UPF3A/',
    dgePrefix = 'results/dge/GTEx/',
    outPrefix = 'results/ds_v_dge_UPF3A/'
  resources: mem_mb = 24000
  log:
      'logs/collect_ds_v_dge.{a_v_b}.log'
  shell:
    '''
    (/software/R-4.1.0-el7-x86_64/bin/Rscript {params.R_script} {params.contrast} {params.gtf} {params.dsPrefix} {params.dgePrefix} {params.outPrefix}) &> {log}
    ls {output}
    '''
         


rule PrepareTablesForHeatmap:
    input:
        rds = 'results/ds_v_dge/{a_v_b}.rds',
        lf1 = 'results/ds/GTEx/{a_v_b}/ds_perind_numers.counts.noise_by_intron.lf1.gz'
    output:
        'tmp/ds_v_dge/{a_v_b}.psi.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.delta_psi.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.logFC_psi.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.logFC_exp.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.cluster_counts.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.cluster_fraction.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.psi_p.tsv.gz',
        'tmp/ds_v_dge/{a_v_b}.exp_p.tsv.gz'
    resources: mem_mb = 24000
    log:
        'logs/process_ds_v_dge.{a_v_b}.log'
    shell:
        '''
        python scripts/PrepareHeatmapTables.py -i {input.rds} &> {log}
        '''
        
rule PrepareCountsTablesForHeatmap:
    input:
        rds = 'results/ds_v_dge/{a_v_b}.rds',
        lf1 = 'results/ds/GTEx/{a_v_b}/ds_perind_numers.counts.noise_by_intron.lf1.gz'
    output:
        'tmp/ds_v_dge/{a_v_b}.cluster_counts_per_tissue.tsv.gz',
    resources: mem_mb = 24000
    log:
        'logs/process_ds_v_dge_counts.{a_v_b}.log'
    shell:
        '''
        python scripts/get_counts_per_cluster_per_tissue.py -i {input.rds} &> {log}
        '''


rule collect_tables_for_heatmap:
    input:
        expand('tmp/ds_v_dge/{a_v_b}.psi.tsv.gz', a_v_b = pairwise_ds_dge),
        expand('tmp/ds_v_dge/{a_v_b}.cluster_counts_per_tissue.tsv.gz', a_v_b = pairwise_ds_dge)
        



rule PrepareTablesForHeatmap_generalized:
    input:
        rds = 'results/ds_v_dge_confounder/rds_files/{a_v_b}.rds',
        lf1 = 'results/ds/GTEx/{a_v_b}/ds_perind_numers.counts.noise_by_intron.lf1.gz'
    output:
        'results/ds_v_dge_confounder/tables/{a_v_b}.psi.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.delta_psi.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.logFC_psi.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.logFC_exp.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.cluster_counts.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.cluster_fraction.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.psi_p.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.exp_p.tsv.gz',
        'results/ds_v_dge_confounder/tables/{a_v_b}.cluster_counts_per_tissue.tsv.gz',
    params:
        'results/ds_v_dge_confounder/tables/{a_v_b}'
    priority: -2
    resources: mem_mb = 24000
    log:
        'logs/process_ds_v_dge_confounder.{a_v_b}.log'
    shell:
        '''
        python scripts/PrepareHeatmapTables_generalized.py {input.rds} {input.lf1} {params} &> {log}
        '''


rule GetPairwiseComparisonsTable_confounder:
    output:
        '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge_confounder/pairwise_comparisons.tab.gz'
    threads: 4
    resources: cpu=1, mem_mb=48000
    log:
        "logs/ds_v_dge_confounder.log"
    priority: -1
    shell:
        """
        (python scripts/GetPairwiseComparisons_generalized.py) &> {log}
        """
        
rule collect_tables_for_fig2:
    input:
        '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge_confounder/pairwise_comparisons.tab.gz',
        expand('results/ds_v_dge_confounder/tables/{a_v_b}.psi.tsv.gz', a_v_b = pairwise_ds_dge),

# -----------------------------------------------------------------------------
#   plot sashimi
# -----------------------------------------------------------------------------

def get_bedgraph_input(wildcards):
    import glob

    bedgraph_dir = '/project2/yangili1/GTEx_v8/bedGraph'
    #tissueTrans = GTEX_BEDGRAPH_TISSUES.get(wildcards.tissue)
    bedFiles = glob.glob(f'{bedgraph_dir}/{wildcards.tissue}/*.bed.gz') # abs paths

    return bedFiles

rule BedgraphToBW:
    '''
    Note the wc for tissue here is the original tissue name, not transformed
    '''
    input: get_bedgraph_input
    output: touch('resources/GTEx/BigWig/{tissue}/done')
    params:
        outPrefix = 'resources/GTEx/BigWig/{tissue}',
        chromSizes = 'resources/hg38_w_chrEBV.chrom.sizes'
    threads: 8
    resources: cpu=8, mem_mb=58000, time=1200
    log:
        'logs/bed2bw.{tissue}.log'
    shell:
        '''
        bedFiles="{input}"

        numJobs=0
        maxJobs={threads}

        for f in ${{bedFiles[@]}}; do
            tmpf={params.outPrefix}/$(basename "$f" .bed.gz).tmp.bed &>> {log}
            tmp_sortedf={params.outPrefix}/$(basename "$f" .bed.gz).tmp.sorted.bed &>> {log}
            outf={params.outPrefix}/$(basename $f .bed.gz).bw &>> {log}
            
            if [[ $numJobs -le $maxJobs ]]; then
                echo scripts/bedGraphToBigWig $f ...
                (bgzip -f -d -c $f > $tmpf && bedtools sort -i $tmpf > $tmp_sortedf && scripts/bedGraphToBigWig $tmp_sortedf {params.chromSizes} $outf && rm $tmpf && rm $tmp_sortedf &
                numJobs=$((numJobs + 1))) &>> {log}
            else
                wait
                numJobs=$((numJobs - 1))
                
                echo bedGraphToBigWig $f ...
                (bgzip -f -d -c $f > $tmpf && bedtools sort -i $tmpf > $tmp_sortedf && scripts/bedGraphToBigWig $tmp_sortedf {params.chromSizes} $outf && rm $tmpf && rm $tmp_sortedf &
                numJobs=$((numJobs + 1))) &>> {log}
            fi
        done

        wait
        echo "All done!"

        '''
      

#-----------------------------------------------------------------------------------------
#   plot sashimi for differential splicing
#-----------------------------------------------------------------------------------------

rule getIntronsForSashimi:
    input:
        rds = '../data/ds_v_dge/{ds_tissue_2}_v_{ds_tissue_1}_data.rds' # separate smk rules to make rds files 
    output: 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt'
    params:
        py_script = 'workflow/scripts/getIntrons.R',
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        minDeltaPsi = 0.2,
        FDR = 0.001,
        minL2FC = 0.58,
    log: 'logs/getIntronsForSashimi/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        Rscript {params.py_script} {input.rds} {params.contrast} {output} {params.minDeltaPsi} {params.FDR} {params.minL2FC} &> {log}
        '''


def getPlotIntrons(wildcards):
    with open(f'plots/sashimi/ds/{wildcards.ds_tissue_2}_v_{wildcards.ds_tissue_1}/introns.txt') as f:
        introns = [x.strip() for x in f.readlines()]
    return introns

rule PrepSashimiDsGtex:
    message: '### Prepare sashimi plots for differential splicing'
    input: 
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt',
        ds_effect_sizes = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_effect_sizes.txt',
        ds_intron_count = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        introns = 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt',
    output:
        flag = touch('plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/prep.done')
        # flag = touch('plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/{plotIntron}.prep.done')
    params:
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        outDir = 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}',
        bwPrefix = 'resources/GTEx/BigWig',
        plotIntrons = getPlotIntrons,
        iniTemplate = 'config/template-sashimi-diffsplice.ini',
        shellTemplate = 'config/template-plot-sashimi-cmd.sh',
        pyscript = 'workflow/scripts/prepSashimi.py',
    threads: 4
    resources: cpu=1, mem_mb=25000, time=1200
    group: 'PrepSashimiDsGtex'
    log: 'logs/PrepSashimiDsGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        #module load parallel

        introns="{params.plotIntrons}"

        cmd="python {params.pyscript}" 
        cmd+=" --contrast {params.contrast} "
        cmd+=" --outDir {params.outDir} "
        cmd+=" --bwPrefix {params.bwPrefix} "
        cmd+=" --dsSampleGroupFile {input.ds_sample_group} "
        cmd+=" --dsEffectFile {input.ds_effect_sizes} "
        cmd+=" --intronCountsFile {input.ds_intron_count} "
        cmd+=" --plotIntron {{}} "
        cmd+=" --iniTemplate {params.iniTemplate} "
        cmd+=" --plotShellTemplate {params.shellTemplate} "
                                 
        parallel -j {threads} "$cmd" ::: $introns &> {log}

        '''

def getPLotSashimiDsGtexParams(wildcards):
    w = wildcards
    folder = f'plots/sashimi/ds/{w.ds_tissue_2}_v_{w.ds_tissue_1}'
    clu = w.plotIntron.split(':')[-1]

    try:
        shellFiles = [os.path.basename(x) for x in glob.glob(f'{folder}/*.sh')]
        shell = [x for x in shellFiles if clu in x][0]
    except:
        shell = None

    return f'{shell}'

rule PlotSashimiDsGtex:
    input: 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/prep.done'
    output: touch('plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}/plot.done')
    params:
        folder = 'plots/sashimi/ds/{ds_tissue_2}_v_{ds_tissue_1}',
        # shell = getPLotSashimiDsGtexParams
    conda: 'pygenometracks'
    group: 'PrepSashimiDsGtex'
    threads: 4
    log: 'logs/PlotSashimiDsGtex/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell: 
        '''
        # module load parallel

        log=$(realpath {log})
        echo "plot sashimi plots" > $log
        cd {params.folder}
        shellFiles=$(ls *.sh)
        parallel -j {threads} "sh {{}}" ::: $shellFiles &>> $log

        #run in login node after: pdfunite *.pdf all_plots.pdf
        '''



#-----------------------------------------------------------------------------------------
#   plot sashimi for selected SRSF genes
#-----------------------------------------------------------------------------------------

rule getIntronsForSashimi_SRSF:
    input: 
        rds = '../data/ds_v_dge/{ds_tissue_2}_v_{ds_tissue_1}_data.rds' # separate smk rules to make rds files 
    output: 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt'
    params:
        R_script = 'workflow/scripts/getIntrons_withGene.R',
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        minDeltaPsi = 0,
        FDR = 1,
        minL2FC = 0,
        gene_list = 'plots/sashimi/SRSF/SRSF.genes'
    log: 'logs/getIntronsForSashimi_SRSF/{ds_tissue_2}_v_{ds_tissue_1}.log'
    shell:
        '''
        Rscript {params.R_script} {input.rds} {params.contrast} {output} {params.minDeltaPsi} {params.FDR} {params.minL2FC} {params.gene_list} &> {log}
        '''

def getPlotIntrons_SRSF(wildcards):
    with open(f'plots/sashimi/SRSF/{wildcards.ds_tissue_2}_v_{wildcards.ds_tissue_1}/introns.txt') as f:
        introns = [x.strip() for x in f.readlines()]
    return introns


use rule PrepSashimiDsGtex as PrepSashimiDsGtex_SRSF with:
    input: 
        ds_sample_group = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_sample_group.txt',
        ds_effect_sizes = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_effect_sizes.txt',
        ds_intron_count = 'results/ds/GTEx/{ds_tissue_2}_v_{ds_tissue_1}/ds_perind_numers.counts.noise_by_intron.gz',
        introns = 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/introns.txt',
    output:
        flag = touch('plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/prep.done')
    params:
        contrast = '{ds_tissue_2}_v_{ds_tissue_1}',
        outDir = 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}',
        bwPrefix = 'resources/GTEx/BigWig',
        plotIntrons = getPlotIntrons_SRSF,
        iniTemplate = 'config/template-sashimi-diffsplice.ini',
        shellTemplate = 'config/template-plot-sashimi-cmd.sh',
        pyscript = 'workflow/scripts/prepSashimi.py',
    log: 'logs/PrepSashimiDsGtex/SRSF/{ds_tissue_2}_v_{ds_tissue_1}.log'


#use rule PlotSashimiDsGtex as PlotSashimiDsGtex_SRSF with:
#    input: 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/prep.done'
#    output: touch('plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}/plot.done')
#    params:
#        folder = 'plots/sashimi/SRSF/{ds_tissue_2}_v_{ds_tissue_1}',
#    conda: 'pygenometracks'
#    group: 'PrepSashimiDsGtex'
#    threads: 4
#    log: 'logs/PlotSashimiDsGtex/SRSF/{ds_tissue_2}_v_{ds_tissue_1}.log'




#-----------------------------------------------------------------------------------------
# ad hoc
#-----------------------------------------------------------------------------------------

rule adhoc_test_ds_step1: # testing differential splicing with neg control
    input: 
        ds_numers_lf1 = 'results/ds/GTEx/Brain-Cerebellum_v_Liver/ds_perind_numers.counts.noise_by_intron.lf1.gz',
    output:
        neg_control_numers = 'results/ds/GTEx/ds_test/BC_v_Liver/perind.numers.gz',
        neg_control_groups = 'results/ds/GTEx/ds_test/BC_v_Liver/sample_group.txt'
    run:
        import gzip
        outf1 = gzip.open(output.neg_control_numers, 'wt')
        outf2 = open(output.neg_control_groups, 'w')
        with gzip.open(input.ds_numers_lf1, 'rt') as f:
            i = 0
            for ln in f:
                if i == 0:
                    header = ln.split()
                    cols = ([x.replace('Brain-Cerebellum', 'BC-group1') for x in header[0:100]] +
                            [x.replace('Brain-Cerebellum', 'BC-group2') for x in header[100:200]])
                    groups = ['BC-group1' for x in header[0:100]] + ['BC-group2' for x in header[100:200]]
                    for c,g in zip(cols, groups):
                        outf2.write(f'{c} {g}\n')
                    outf1.write(' '.join(cols) + '\n')
                if i > 0:
                    outln = ln.split()[:201]
                    outf1.write(' '.join(outln) + '\n')
                i += 1

        outf1.close()
        outf2.close()


use rule RunLeafcutterDiffSplicingGtex as adhoc_test_ds_step2 with:
    input:
        ds_numers_lf1 = 'results/ds/GTEx/ds_test/BC_v_Liver/perind.numers.gz',
        ds_sample_group = 'results/ds/GTEx/ds_test/BC_v_Liver/sample_group.txt'
    output:
        flag = touch('results/ds/GTEx/ds_test/BC_v_Liver/ds.done')
        # produces two files:
        # 1. {outprefix}_effect_sizes.txt
        # 2. {outprefix}_manual_ds_cluster_significance.txt
    params:
        Rscript = 'workflow/submodules/leafcutter/scripts/leafcutter_ds.R', 
        outprefix = 'results/ds/GTEx/ds_test/BC_v_Liver/ds', # note you need to include path!
        MIN_SAMPLES_PER_INTRON = 5,
        MIN_SAMPLES_PER_GROUP = 3,
        MIN_COVERAGE = 5
    log: 'results/ds/GTEx/ds_test/BC_v_Liver/log'



#-----------------------------------------------------------------------------------------
# data for plotting heatmap of unprodcutive splicing
#-----------------------------------------------------------------------------------------
# only plot several tissues from brain
#"../../../analysis/2024-04-15-GTEx-psi-heatmap3.ipynb"

rule PrepGTExHeatmapData:
    '''
    Intro - Plot heatmap of unproductive splicing introns across GTEx tissues
    heatmap data:

    rows: top unproductive introns
    columns: GTEx tissues
    values: PSI values
    howto:

    pick top unproductive introns:
    get top deltaPSI unprod introns from a set of contrasts
    union these introns
    get PSI values for these introns across GTEx tissues

    '''
    message: '### Prepare data for plotting heatmap of unproductive splicing'
    output: 'plotdata/gtex-psi-heatmap/heatmap_data.rds'
    params:
        R_script = 'workflow/scripts/prepGTExHeatmapData.R',
        input_rds_dir = 'plotdata/ds_v_dge',
        out_dir = 'plotdata/GTEx-PSI-Heatmap',
        FDR_ds = 1e-3,
        FDR_dge = 0.05,
        minDeltaPSI = 0.1,
    shell:
        '''
        Rscript {params.R_script} {params.input_rds_dir} {output} {params.FDR_ds} {params.FDR_dge} {params.minDeltaPSI}
        ls {output}
        '''

rule GetPairwiseComparisonsTable:
    output:
        '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/ds_v_dge_pairwise_comparisons.tab.gz'
    threads: 4
    resources: cpu=1, mem_mb=48000
    log:
        "logs/ds_v_dge.log"
    shell:
        """
        (python scripts/GetPairwiseComparisons.py) &> {log}
        """