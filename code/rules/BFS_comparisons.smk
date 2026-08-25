rule RINLeafCutter2:
    output:
        ds_numers = 'results/RIN/{RIN}/leafcutter2.junction_counts.const.gz',
    params:
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        junc_files = '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/config/GTEx_{RIN}_files.txt',
        res_dir = 'results/RIN/{RIN}/',
        other_params = '--includeconst'
    log:
        "logs/RIN/{RIN}/tcga_lc2.log"
    resources:
        mem_mb = 52000
    wildcard_constraints:
        RIN = 'RIN5-6|RIN6-7|RIN7-8|RIN8-9|RIN9-10'
    shell:
        """
        (python /project/yangili1/cfbuenabadn/leafcutter2/scripts/leafcutter2.py \
            -j {params.junc_files} \
            -r {params.res_dir} \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} {params.other_params}) &> {log}
        """

rule CollectRINLeafCutter2:
    input:
        expand('results/RIN/{RIN}/leafcutter2.junction_counts.const.gz', 
                RIN = ['RIN5-6', 'RIN6-7', 'RIN7-8', 'RIN8-9', 'RIN9-10'])


rule LeafCutter2LCLs:
    output:
        ds_numers = 'results/LCLs/{RNAtype}/leafcutter2.junction_counts.const.gz',
    params:
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        junc_files = '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/config/{RNAtype}.LCL_juncs.txt',
        res_dir = 'results/LCLs/{RNAtype}/',
        other_params = '--includeconst'
    log:
        "logs/LCLs/{RNAtype}/lc2.log"
    resources:
        mem_mb = 52000
    wildcard_constraints:
        RNAtype = 'steady_state|nascent'
    shell:
        """
        (python /project/yangili1/cfbuenabadn/leafcutter2/scripts/leafcutter2.py \
            -j {params.junc_files} \
            -r {params.res_dir} \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} {params.other_params}) &> {log}
        """


rule LeafCutter2LCLs_recursive:
    output:
        ds_numers = 'results/LCLs_recursive/{RNAtype}/leafcutter2.junction_counts.const.gz',
    params:
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        junc_files = '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/config/{RNAtype}.LCL_juncs_updated.txt',
        res_dir = 'results/LCLs_recursive/{RNAtype}/',
        other_params = '--includeconst'
    log:
        "logs/LCLs/{RNAtype}/lc2.log"
    resources:
        mem_mb = 52000
    wildcard_constraints:
        RNAtype = 'steady_state|nascent'
    shell:
        """
        (leafcutter2 \
            -j {params.junc_files} \
            -r {params.res_dir} \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} {params.other_params}) &> {log}
        """


rule LeafCutter2HeLa_recursive:
    output:
        ds_numers = 'results/HeLa_recursive/{condition}/leafcutter2.junction_counts.const.gz',
    params:
        gtf = config['annotation']['gtf']['v43'],
        genome = config['genome38'],
        max_juncs = 1000, # maximum number of introns per gene
        res_dir = 'results/HeLa_recursive/{condition}/',
        other_params = '--includeconst',
        condition_juncs = 'resources/HeLa/juncs/{condition}.',
    log:
        "logs/HeLa_recursive/{condition}/lc2.log"
    resources:
        mem_mb = 52000
    wildcard_constraints:
        condition = 'controls|HeLa_dKD'
    shell:
        """
        (leafcutter2 \
            -j <(cat <(ls {params.condition_juncs}*tsv.gz)) \
            -r {params.res_dir} \
            -A {params.gtf} \
            -G {params.genome} \
            --max_juncs {params.max_juncs} {params.other_params}) &> {log}
        """



rule CollectLeafCutter2LCLs:
    input:
        expand('results/LCLs/{RNAtype}/leafcutter2.junction_counts.const.gz', 
                RNAtype = ['steady_state', 'nascent'])


rule LC2Dylan:
    output:
        ds_numers = 'results/Dylan/leafcutter2.junction_counts.gz',
    params:
        gtf = "/project/yangili1/dylan_stermer/ReferenceGenomes/GRCh38_GencodeRelease44Comprehensive/Reference.gtf",
        genome = "/project/yangili1/dylan_stermer/ReferenceGenomes/GRCh38_GencodeRelease44Comprehensive/Reference.fa",
        max_juncs = 1000, # maximum number of introns per gene
        junc_files = 'config/dylan_juncs.txt',
        res_dir = 'results/Dylan/',
        other_params = '-t transcript_type'
    log:
        "logs/Dylan/lc2.log"
    resources:
        mem_mb = 52000
    shell:
        """
        (python /project/yangili1/cfbuenabadn/leafcutter2/scripts/leafcutter2.py \
            -j {params.junc_files} \
            -r {params.res_dir} \
            -o leafcutter2 \
            -A {params.gtf} \
            -G {params.genome} \
            {params.other_params}) &> {log}
        """


rule MR2SLS:
    output:
        'results/mr_2sls/{tissue}_PR.tab.gz',
        'results/mr_2sls/{tissue}_UP.tab.gz'
    log:
        "logs/MR/{tissue}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        RNAtype = '|'.join(tissues)
    shell:
        """
        (python scripts/MR2SLS.py {wildcards.tissue}) &> {log}
        """

rule collectMR:
    input:
        expand('results/mr_2sls/{tissue}_PR.tab.gz', tissue = tissues)