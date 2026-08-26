
N_PermutationChunks = 2
ChunkNumbers = range(0, 1+N_PermutationChunks) 

rule sQTLs_perm:
    '''
        already ranknorm normalized, do not need the normal option 
    '''
    message: 'Map QTL using permutation pass'
    input: 
        phenoPrep = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/{datasource}/{group}/separateNoise/done',
        vcf = '/project/yangili1/cfbuenabadn/SpliFi/code/results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/{datasource}/{group}/separateNoise/{chrom}_CovMatrix.txt',
    output: temp('results/sqtl/{datasource}/{group}/cis_{window}/perm/chunks/{chrom}.{QTLTools_chunk_n}.txt')
    log: 'logs/NominalQTL_{datasource}_{group}_{window}_{chrom}.{QTLTools_chunk_n}.log'
    params:
        cis_window = lambda w: int(w.window),
        pheno = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter.qqnorm_{chrom}.gz',
        NChunks = N_PermutationChunks
    resources: mem_mb = 54000
    shell:
        '''
        module unload gsl && module load gsl/2.5
        {config[QTLtools]} cis \
            --std-err \
            --seed 123 \
            --permute 1000 \
            --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} \
            --vcf {input.vcf} --bed {params.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        '''


# top sQTL's nominal pass stats, used for beta-beta plot
rule NominalQTL:
    message: '### Map QTL using nominal pass'
    input: 
        phenoPrep = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/{datasource}/{group}/separateNoise/done',
        vcf = '/project/yangili1/cfbuenabadn/SpliFi/code/results/geno/{datasource}/{group}/{chrom}.vcf.gz',
        cov = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/{datasource}/{group}/separateNoise/{chrom}_CovMatrix.txt',
    output: temp('results/sqtl/{datasource}/{group}/cis_{window}/nom/chunks/{chrom}.{QTLTools_chunk_n}.txt')
    log: 'logs/NominalQTL_{datasource}_{group}_{window}_{chrom}.{QTLTools_chunk_n}.log'
    params:
        cis_window = lambda w: int(w.window) + 10000, # add 10kb to the window
        pheno = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter.qqnorm_{chrom}.gz',
        NChunks = N_PermutationChunks
    resources: mem_mb = 54000
    shell:
        '''
        module unload gsl && module load gsl/2.5
        {config[QTLtools]} cis \
            --std-err \
            --seed 123 \
            --nominal 1 \
            --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} \
            --vcf {input.vcf} --bed {params.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        '''
        
rule Gather_QTLtools_cis_pass:
    input:
        expand('results/sqtl/{{datasource}}/{{group}}/cis_{{window}}/{{pass}}/chunks/{{chrom}}.{QTLTools_chunk_n}.txt', QTLTools_chunk_n=ChunkNumbers)
    output:
       temp('results/sqtl/{datasource}/{group}/cis_{window}/{pass}/{chrom}.temp.txt'),
    log:
        "logs/Gather_QTLtools_cis_pass/sqtl/{datasource}/{group}/cis_{window}/{pass}/{chrom}.log"
    shell:
        """
        (cat {input} - > {output}) &> {log}
        """

rule tabixNominalPassQTLResults:
    """
    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
    for easy access with tabix
    """
    input:
        'results/sqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.temp.txt',
    params:
        sort_temp = '-T ' + config['scratch'][:-1]
    output:
        txt = 'results/sqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.txt.gz',
        tbi = 'results/sqtl/{datasource}/{group}/cis_{window}/nom/{chrom}.txt.gz.tbi',
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "logs/tabixNominalPassQTLResults/{datasource}/{group}/cis_{window}/nom/{chrom}.log"
    shadow: "shallow"
    shell:
        """
        (cat <(cat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(cat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k9,9 -k10,10n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
        tabix -b 10 -e10 -s9 {output.txt} &>> {log}
        """

rule collectQTLs:
    input:
        expand('results/sqtl/GTEx/{group}/cis_100000/nom/{chrom}.txt.gz', 
        group = tissues,
        chrom = ['chr' + str(x) for x in range(1, 23)])


rule MergeAndAnnotatePerm_sQTLs:
    """
    This rule merges sQTL permutation pass, annotates junctions by gene overlaps, 
    and calculates qvalue again on the whole run.
    """
    input:
        perm = expand(
        'results/sqtl/GTEx/{{tissue}}/cis_100000/perm/{chrom}.temp.txt',
        chrom = ['chr' + str(x) for x in range(1, 23)]
        ),
        gtf = "annotations/gencode.v26.GRCh38.genes.csv"
    output:
        'results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz'
    log:
        "logs/MergeAnnotateAndAddQ/{tissue}.log"
    resources: mem_mb = 24000
    shell:
        """
        (/software/R-4.1.0-el7-x86_64/bin/Rscript scripts/annotate_sQTL_Perm.R {wildcards.tissue} {input.gtf} {output}) &> {log}
        """
    
rule collect_perm_sQTLs:
    input:
        expand('results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz',
        tissue = tissues)





##### Ru sQTLs ####



rule PhenotypePCsRu:
    input:
        'AD_QTLs/carlos_processed/chr10.qqnorm.bed.gz',
    output:
        'AD_QTLs/carlos_processed/chr10.qqnorm.bed.pca',
    log:
        "logs/pca.log"
    resources:
        mem_mb = 24000
    conda:
        "../envs/r_essentials.yml"
    shell:
        """
        Rscript scripts/PermuteAndPCA.R {input} {output} &> {log}
        """
        
def getMAF(wildcards):
    if wildcards.MAF == 'MAF05':
        return '0.05'
    elif wildcards.MAF == 'ALL':
        return '1.0'
        
rule vcfMAF:
    input: 
        vcf = 'AD_QTLs/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.10.vcf.vcf',
    output: 
        vcf = 'AD_QTLs/carlos_processed/chr10.{MAF}.vcf.gz',
        tbi = 'AD_QTLs/carlos_processed/chr10.{MAF}.vcf.gz.tbi'
    log: 'logs/vcfmaf.{MAF}.log'
    params: getMAF
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    wildcard_constraints: MAF='ALL|MAF05'
    shell:
        '''
        bcftools view \
            --threads {threads} \
            -i "AF >= {params}" \
            -Oz -o {output.vcf} \
            {input.vcf} &> {log}
        
        bcftools index --threads {threads} --tbi {output.vcf} &>> {log}
        '''
        
rule make_standard_vcf:
    input:
        'AD_QTLs/carlos_processed/chr10.{MAF}.vcf.gz',
    output:
        'AD_QTLs/carlos_processed/chr10.{MAF}.standard.vcf',
    log:
        'logs/make_standard_vcf.{MAF}.log'
    wildcard_constraints: MAF='ALL|MAF05'
    shell:
        """
        python scripts/prepare_Ru_vcf.py {input} {output} &> {log}
        """
        
rule vcfChange:
    input: 
        vcf = 'AD_QTLs/carlos_processed/chr10.{MAF}.standard.vcf',
    output: 
        vcf = 'AD_QTLs/carlos_processed/chr10.{MAF}.standard.vcf.gz',
        tbi = 'AD_QTLs/carlos_processed/chr10.{MAF}.standard.vcf.gz.tbi',
    wildcard_constraints: MAF='ALL|MAF05'
    log: 'logs/vcfchange.{MAF}.log'
    threads: 4
    resources: time=2000, mem_mb=15000, cpu=4
    shell:
        '''
        (cat {input} | bgzip /dev/stdin -c > {output.vcf}) &> {log}
        bcftools index --threads {threads} --tbi {output.vcf} &>> {log}
        '''
        
rule SortQTLtoolsPhenotypeTable:
    input:
        'AD_QTLs/carlos_processed/chr10.qqnorm.bed.gz',
    output:
        bed = 'AD_QTLs/carlos_processed/chr10.qqnorm.sorted.bed.gz',
        tbi = 'AD_QTLs/carlos_processed/chr10.qqnorm.sorted.bed.gz.tbi',
    log:
        "logs/tabixru.log"
    resources:
        mem_mb = 24000
    shell:
        """
        (zcat {input} | bedtools sort -header -i - | bgzip /dev/stdin -c > {output.bed}) &> {log}
        (tabix -p bed {output.bed}) &>> {log}
        """


rule QTLTools_Ru:
    input: 
        pheno = 'AD_QTLs/carlos_processed/chr10.qqnorm.sorted.bed.gz',
        tbi = 'AD_QTLs/carlos_processed/chr10.qqnorm.sorted.bed.gz.tbi',
        vcf = 'AD_QTLs/carlos_processed/chr10.{MAF}.standard.vcf.gz',
        vcftbi = 'AD_QTLs/carlos_processed/chr10.{MAF}.standard.vcf.gz.tbi',
        cov = 'AD_QTLs/carlos_processed/chr10.qqnorm.bed.pca',
    output: temp('AD_QTLs/results/{MAF}.{QTLTools_chunk_n}.txt')
    log: 'logs/Ru.{QTLTools_chunk_n}.{MAF}.log'
    wildcard_constraints: MAF='ALL|MAF05'
    params:
        cis_window = 100000,
        NChunks = 10
    resources: mem_mb = 12000
    shell:
        """
        module unload gsl && module load gsl/2.5
        {config[QTLtools]} cis \
            --std-err \
            --seed 123 \
            --permute 1000 \
            --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} \
            --vcf {input.vcf} --bed {input.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} &> {log}
        """

RU_ChunkNumbers = range(0, 1+10)

rule Gather_QTLtools_Ru:
    input:
        expand('AD_QTLs/results/{{MAF}}.{QTLTools_chunk_n}.txt', QTLTools_chunk_n=RU_ChunkNumbers )
    output:
        "AD_QTLs/results/{MAF}.PermutationPass.txt.gz"
    wildcard_constraints: MAF='ALL|MAF05'
    log:
        "logs/Gather_QTLtools_{MAF}..log"
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """
        
        
rule AddQValueToPermutationPass:
    input:
        "AD_QTLs/results/{MAF}.PermutationPass.txt.gz"
    output:
        table = "AD_QTLs/results/{MAF}.PermutationPass.FDR_Added.txt.gz",
    log:
        "logs/{MAF}.Ru.Q.log"
    wildcard_constraints: MAF='ALL|MAF05'
    conda:
        "../envs/r_essentials.yml"
    priority:
        10
    shell:
        """
        Rscript scripts/AddQvalueToQTLtoolsOutput.Ru.R {input} {output} &> {log}
        """
        
        
rule collectRu:
    input:
        expand('AD_QTLs/results/{MAF}.PermutationPass.FDR_Added.txt.gz', MAF = ['ALL', 'MAF05'])


## -----------------------------------------------------------------------------
##   Per-tissue sQTL summary tables consumed by Figure 4
##
##   Both were previously built by cells in ../analysis/QTL_analysis.ipynb whose
##   to_csv() calls were commented out, written to the notebook's own directory
##   and moved into analysis_files/ by hand -- so neither had a traceable origin.
##   See code_report.txt.
##
##   Note: get_perm_counts() also reads the leafcutter2 noise tables from the
##   separate SpliFi project
##     /project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/
##   which is outside this workflow, so it is not declared as an input.
## -----------------------------------------------------------------------------

rule MakeTotalsQTLTable:
    """Per-tissue counts of significant sQTLs by cluster type (Figure 4a)."""
    input:
        perm = expand('results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz',
                      tissue = tissues)
    output:
        'analysis_files/Total_sQTLs.tsv.gz'
    params:
        tissues = ' '.join(tissues)
    resources:
        mem_mb = 24000
    log:
        'logs/MakeTotalsQTLTable.log'
    shell:
        """
        (python scripts/MakeTotalsQTLTable.py {output} {params.tissues}) &> {log}
        """


rule MakesQTLStatsTable:
    """Per-tissue sQTL-vs-eQTL effect-size correlations (Figure 4c, Supp. Table 9).

    Slow: one tabix query against the eQTL nominal pass per significant sQTL,
    over 49 tissues x 3 sQTL classes.
    """
    input:
        perm = expand('results/sqtl/GTEx/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz',
                      tissue = tissues),
        eqtl_nom = expand('results/eqtl/GTEx/{tissue}/cis_100000/nom/{chrom}.txt.gz',
                          tissue = tissues,
                          chrom = ['chr' + str(x) for x in range(1, 23)])
    output:
        'analysis_files/sQTL_stats.tsv.gz'
    params:
        tissues = ' '.join(tissues)
    resources:
        mem_mb = 24000
    log:
        'logs/MakesQTLStatsTable.log'
    shell:
        """
        (python scripts/MakesQTLStatsTable.py {output} {params.tissues}) &> {log}
        """
