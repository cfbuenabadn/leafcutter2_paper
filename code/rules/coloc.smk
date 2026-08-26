
def GetQTLtoolsWindowFlag(wildcards):
    if wildcards.FeatureCoordinatesRedefinedFor == '.ForGWASColoc':
        return "--window 1"
    else:
        return "--window 10000"
        
def GetPhenotypeInput(wildcards):
    if wildcards.pheno == 'leafcutter':
        return '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/{tissue}/separateNoise/leafcutter.qqnorm_{chrom}.gz'
    else:
        return '/project/yangili1/cfbuenabadn/SpliFi/code/results/eqtl/GTEx/{tissue}/qqnorm.sorted.{chrom}.bed.gz'
        


rule MakePhenotypeTableToColocFeaturesWithGWASLoci:
    """
    To colocalize molQTLs with gwas we need summary stats (betas and se) for
    every snp in the same window as a gwas locus. Here, I will take a phenotype
    table, find the features that intersect the gwas locus (1MB window centered
    on lead snp, defined in the input file), and output a phenotype table with
    every intersection with coordinates redefined to be the 1MB gwas locus
    window. Running QTLtools on this will help me get the necessary summary
    stats in those windows
    """
    input:
        bed = GetPhenotypeInput,
        fai = "resources/Annotations/GRCh38.primary_assembly.genome.fa.fai",
        loci = "resources/gwas/LeadSnpWindows.bed"
    output:
        bed = "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.gz",
        tbi = "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.gz.tbi"
    resources:
        mem_mb = 48000
    wildcard_constraints:
        pheno = 'leafcutter|expression',
    log:
        "/scratch/midway2/cnajar/logs/MakePhenotypeTableToColocFeaturesWithGWASLoci/{tissue}.{chrom}.{pheno}.log"
    shell:
        """
        (cat <(zcat {input.bed} | head -1) <(  bedtools intersect  -wo -a {input.bed} -b {input.loci} -sorted | awk -F'\\t' -v OFS='\\t' '{{$4=$4":"$(NF-1); $5=$(NF-1); $2=$(NF-3); $3=$(NF-2); print $0}}' | rev | cut -f 6- | rev ) | tr ' ' '\\t' | bedtools sort -i - -header | bgzip /dev/stdin -c > {output.bed} ) &> {log};
        tabix -p bed {output.bed}
        """
        
        
        
        
        
rule PhenotypePCs:
    """
    QTLtools format expression PCs as covariates
    including the number of PCs that explain more
    variance then when the phenotype table is
    permuted
    """
    input:
        "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.gz"
    output:
        "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.pca"
    log:
        "/scratch/midway2/cnajar/logs/QTLs/PhenotypePCs/{tissue}.{chrom}.{pheno}.log"
    resources:
        mem_mb = 24000
    wildcard_constraints:
        pheno = 'leafcutter|expression',
        tissue = '|'.join(tissues),
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)])
    shell:
        """
        Rscript scripts/PermuteAndPCA.R {input} {output} &> {log}
        """
        
    
    



N_PermutationChunks_coloc = 3
ChunkNumbers_coloc = range(0, 1+N_PermutationChunks_coloc) 


rule coloc_QTLs:
    '''
        already ranknorm normalized, do not need the normal option 
    '''
    message: 'Map QTL using permutation pass'
    input: 
        pheno = "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.gz",
        vcf = '/project/yangili1/cfbuenabadn/SpliFi/code/results/geno/GTEx/{tissue}/{chrom}.vcf.gz',
        cov = "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.pca",
    output: temp('results/coloc/qtls/{tissue}/{pheno}.PermutationPass_chunks/{chrom}.{QTLTools_chunk_n}.txt')
    log: '/scratch/midway2/cnajar/logs/PermutationQTL_GTEx_{tissue}_{pheno}_{chrom}.{QTLTools_chunk_n}.log'
    params:
        cis_window = 1,
        NChunks = N_PermutationChunks_coloc
    resources: mem_mb = 58000
    shell:
        '''
        module unload gsl && module load gsl/2.5
        {config[QTLtools]} cis \
            --std-err \
            --seed 123 \
            --permute 1000 \
            --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} \
            --vcf {input.vcf} --bed {input.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        '''
        
rule Gather_QTLtools_cis_pass_coloc:
    input:
        expand('results/coloc/qtls/{{tissue}}/{{pheno}}.PermutationPass_chunks/{{chrom}}.{QTLTools_chunk_n}.txt', 
        QTLTools_chunk_n=ChunkNumbers_coloc)
    output:
       'results/coloc/qtls/{tissue}/{pheno}.PermutationPass/{chrom}.txt.gz',
    log:
        "/scratch/midway2/cnajar/logs/Gather_QTLtools_cis_pass/qtls/{tissue}/{pheno}.PermutationPass/{chrom}.log"
    shell:
        """
        (cat {input} | gzip - > {output}) &> {log}
        """

rule GatherLeafcutter_coloc:
    input:
        expand('results/coloc/qtls/{tissue}/leafcutter.PermutationPass/{chrom}.txt.gz', tissue=tissues,
        chrom = ['chr' + str(x) for x in range(1, 23)])

        
rule GatherLeafcutter_coloc_tissues:
    input:
        expand('results/coloc/qtls/{tissue}/{pheno}.PermutationPass/{chrom}.txt.gz', 
        pheno = ['leafcutter', 'expression'],
        tissue=['Lung', 'Adipose-Visceral_Omentum', 'Muscle-Skeletal', 'Breast-MammaryTissue'],
        chrom = ['chr' + str(x) for x in range(1, 23)])
        

rule MergeAndAnnotatePerm_QTLs_coloc:
    """
    This rule merges sQTL permutation pass, annotates junctions by gene overlaps, 
    and calculates qvalue again on the whole run.
    """
    input:
        perm = expand(
         'results/coloc/qtls/{{tissue}}/{{pheno}}.PermutationPass/{chrom}.txt.gz',
        chrom = ['chr' + str(x) for x in range(1, 23)]
        ),
    output:
        'results/coloc/qtls/{tissue}/{pheno}.PermutationPass.FDR_Added.txt.gz'
    log:
        "logs/MergeAnnotateAndAddQ/{tissue}.{pheno}.log"
    resources: mem_mb = 24000
    shell:
        """
        (/software/R-4.1.0-el7-x86_64/bin/Rscript scripts/AddQvalueToQTLtoolsOutput.R {wildcards.tissue} {wildcards.pheno} {output}) &> {log}
        """
    

rule SelectSignificantTraitsForGWASColoc:
    input:
        bed = "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.gz",
        qtls = 'results/coloc/qtls/{tissue}/{pheno}.PermutationPass.FDR_Added.txt.gz'
    output:
        bed = temp('results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.significant_only.qqnorm.bed'),
        bgz = 'results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.significant_only.qqnorm.bed.gz',
        tbi = 'results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.significant_only.qqnorm.bed.gz.tbi'
    log:
        '/scratch/midway3/cnajar/logs/QTLs/get_significant_traits_.ForGWASColoc.{tissue}/{chrom}.{pheno}.log'
    wildcard_constraints:
        pheno = 'leafcutter|expression',
        tissue = '|'.join(tissues),
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)])
    shell:
        """
        python scripts/filter_permutation_pass_forGWAScoloc.py {input.bed} {input.qtls} {output.bed} &> {log};
        (bgzip {output.bed} -c > {output.bgz}) &>> {log}
        tabix -p bed {output.bgz} &>> {log}
        """

def GetInputPhenoForNominalPass(wildcards):
    if wildcards.pheno == 'leafcutter':
        return "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.significant_only.qqnorm.bed.gz"
    else:
        return "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.gz"

rule coloc_QTLs_nominal:
    '''
        already ranknorm normalized, do not need the normal option 
    '''
    message: 'Map QTL using permutation pass'
    input: 
        pheno = GetInputPhenoForNominalPass, #"results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.significant_only.qqnorm.bed.gz",
        vcf = '/project/yangili1/cfbuenabadn/SpliFi/code/results/geno/GTEx/{tissue}/{chrom}.vcf.gz',
        cov = "results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.qqnorm.bed.pca",
    output: 'results/coloc/qtls/{tissue}/{pheno}.NominalPass_chunks/{chrom}.{QTLTools_chunk_n}.txt'
    log: '/scratch/midway2/cnajar/logs/NominalQTL_GTEx_{tissue}_{pheno}_{chrom}.{QTLTools_chunk_n}.log'
    params:
        cis_window = 1,
        NChunks = N_PermutationChunks_coloc
    resources: mem_mb = 58000
    shell:
        '''
        module unload gsl && module load gsl/2.5
        {config[QTLtools]} cis \
            --std-err \
            --seed 123 \
            --nominal 1 \
            --chunk {wildcards.QTLTools_chunk_n} {params.NChunks} \
            --vcf {input.vcf} --bed {input.pheno} --cov {input.cov}  --out {output} \
            --window {params.cis_window} &> {log}
        if [ ! -f {output} ]
        then
            touch {output}
        fi
        '''
        
rule Gather_QTLtools_nominal_pass_coloc:
    input:
        expand('results/coloc/qtls/{{tissue}}/{{pheno}}.NominalPass_chunks/{{chrom}}.{QTLTools_chunk_n}.txt', 
        QTLTools_chunk_n=ChunkNumbers_coloc)
    output:
       bed = temp('results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.gz'),
       #tbi = 'results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.gz.tbi',
    log:
        "/scratch/midway2/cnajar/logs/Gather_QTLtools_nominal_pass/qtls/{tissue}/{pheno}.PermutationPass/{chrom}.log"
    shell:
        """
        (cat {input} | gzip - > {output.bed}) &> {log};
        #(tabix -p bed {output.bed}) &> {log}
        """
        
rule tabixNominalPassQTLResultsForGWASColoc:
    """
    Convert QTLtools output to tab delimited bgzipped and tabix indexed files
    for easy access with tabix
    """
    input:
        "results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.gz"
    wildcard_constraints:
        pheno = 'leafcutter|expression',
        tissue = '|'.join(tissues),
        chrom = '|'.join(['chr' + str(x) for x in range(1, 23)])
    params:
        sort_temp = '-T /scratch/midway2/cnajar/'
        # sort_temp = ""
    output:
        txt = "results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz",
        tbi = "results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz.tbi"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "/scratch/midway2/cnajar/logs/tabixNominalPassQTLResults/{tissue}/{pheno}.NominalPass/{chrom}.log"
    shadow: "shallow"
    shell:
        """
        (cat <(zcat {input} | head -1 | perl -p -e 'printf("#") if $. ==1; s/ /\\t/g') <(zcat {input} | awk 'NR>1' |  perl -p -e 's/ /\\t/g' | sort {params.sort_temp} -k9,9 -k10,10n  ) | bgzip /dev/stdin -c > {output.txt}) &> {log}
        tabix -b 10 -e10 -s9 {output.txt} &>> {log}
        """

rule tabixStatsForColoc:
    input:
        "resources/gwas/StatsForColoc/{trait}.standardized.txt.gz"
    output:
        txt = "resources/gwas/StatsForColoc/{trait}.standardized.txt.tabix.gz",
        tbi = "resources/gwas/StatsForColoc/{trait}.standardized.txt.tabix.gz.tbi"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    log:
        "/scratch/midway2/cnajar/logs/tabixStatsForColoc/{trait}.log"
    shadow: "shallow"
    shell:
        """
        (zcat {input} | bgzip -c > {output.txt}) &> {log};
        (tabix -S 1 -s 2 -b 3 -e 3 -h {output.txt}) &> {log}
        """
        
rule GatherGWASTraits_tabix:
    input:
        expand("resources/gwas/StatsForColoc/{trait}.standardized.txt.tabix.gz", trait = gwas_traits)
        
rule GatherLeafcutter_nominal_coloc_tissues:
    input:
        expand('results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz', 
        pheno = ['leafcutter', 'expression'],
        tissue=['Lung', 'Adipose-Visceral_Omentum', 'Muscle-Skeletal', 'Breast-MammaryTissue'],
        chrom = ['chr' + str(x) for x in range(1, 23)])
        
rule GatherLeafcutter_nominal_coloc_all_tissues:
    input:
        expand('results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz', 
        pheno = ['leafcutter', 'expression'],
        tissue=tissues,
        chrom = ['chr' + str(x) for x in range(1, 23)])
        
rule GatherLeafcutter_nominal_coloc_all_tissues_expression:
    input:
        expand('results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz', 
        pheno = ['expression'],
        tissue=tissues,
        chrom = ['chr' + str(x) for x in range(1, 23)])
        
        
def getGWAS_region(wildcards):
    row = leadSNPs.loc[leadSNPs.gwas_loci == wildcards.gwas_loci]
    chrom = row.iloc[0].chrom
    start = str(row.iloc[0].start)
    end = str(row.iloc[0].end)
    region = chrom + ':' + start + '-' + end
    return region
    
rule ColocGWASLoci:
    output:
        'results/coloc/colocboost/data/{gwas_loci}.rds',
        'results/coloc/colocboost/results/{gwas_loci}.rds'
    resources:
        mem_mb = 58000
    log:
        "/scratch/midway2/cnajar/logs/colocboost/{gwas_loci}.log"
    params:
        region = getGWAS_region
    wildcard_constraints:
        gwas_loci = '|'.join(gwas_loci)
    shell:
        """
        module load gsl/2.7 && module load R/4.1.0;
        /software/R-4.1.0-el8-x86_64/bin/Rscript scripts/colocboost_summarystats.R {params.region} {wildcards.gwas_loci} &> {log}
        """

rule collect_colocboost:
    input:
        expand('results/coloc/colocboost/results/{gwas_loci}.rds', gwas_loci = gwas_loci[:10])
    
def getGWAS_region_hyprcoloc(wildcards):
    row = leadSNPs.loc[leadSNPs.gwas_loci == wildcards.gwas_loci]
    chrom = row.iloc[0].chrom
    start = str(int(row.iloc[0].start)+450000)
    end = str(int(row.iloc[0].end)-450000)
    region = chrom + ':' + start + '-' + end
    return region    
  
rule HyprcolocGWASLoci:
    output:
        rds = 'results/coloc/hyprcoloc_results/rds/{gwas_loci}.rds',
        tsv = 'results/coloc/hyprcoloc_results/tables/temp/{gwas_loci}.tsv'
    resources:
        mem_mb = 58000
    log:
        "/scratch/midway2/cnajar/logs/hyprcoloc/{gwas_loci}.log"
    params:
        region = getGWAS_region #_hyprcoloc
    wildcard_constraints:
        gwas_loci = '|'.join(gwas_loci)
    shell:
        """
        module load gsl/2.7 && module load R/4.1.0;
        /software/R-4.1.0-el8-x86_64/bin/Rscript scripts/hyprcoloc_gwas.R {params.region} {wildcards.gwas_loci} &> {log}
        """
    
rule collect_hyprcoloc_table:
    input:
        expand('results/coloc/hyprcoloc_results/tables/temp/{gwas_loci}.tsv', gwas_loci=gwas_loci)
    output:
        'results/coloc/hyprcoloc_results/tables/hyprcoloc_results.tsv.gz'
    resources:
        mem_mb = 58000
    log:
        "/scratch/midway2/cnajar/logs/hyprcoloc.log"
    shell:
        """
        (cat {input} | gzip - > {output}) %>% {log}
        """
    
rule filter_hyprcoloc_table:
    """Drop the hyprcoloc iterations that found no colocalization.

    hyprcoloc emits one row per iteration and writes the literal 'None' in the
    traits column when that iteration colocalized nothing; 532,242 of the
    669,148 rows are such rows. Everything downstream (Figure 4d, and the
    UP_colocs / lfc_colocs / expr_colocs flags derived in
    analysis/Figure4/Figure4_helpers.py) works from the colocalizing subset.

    This step previously happened interactively and was never recorded, leaving
    hyprcoloc_results_filtered.tsv.gz with no traceable origin. The awk below
    reproduces the existing file byte for byte.
    """
    input:
        'results/coloc/hyprcoloc_results/tables/hyprcoloc_results.tsv.gz'
    output:
        'results/coloc/hyprcoloc_results/tables/hyprcoloc_results_filtered.tsv.gz'
    log:
        'logs/filter_hyprcoloc_table.log'
    shell:
        """
        (zcat {input} | awk -F'\t' 'NR==1 || $2 != "None"' | gzip -c > {output}) &> {log}
        """


rule collect_hyprcoloc:
    input:
        expand('results/coloc/hyprcoloc_results/rds/{gwas_loci}.rds', gwas_loci = gwas_loci)
    


rule HyprcolocGWASLoci_allTissues:
    output:
        rds = 'results/coloc/hyprcoloc_results_allTissues/rds/{gwas_loci}.rds',
        tsv = 'results/coloc/hyprcoloc_results_allTissues/tables/temp/{gwas_loci}.tsv'
    resources:
        mem_mb = 58000
    log:
        "/scratch/midway2/cnajar/logs/hyprcoloc_all/{gwas_loci}.log"
    params:
        region = getGWAS_region #_hyprcoloc
    wildcard_constraints:
        gwas_loci = '|'.join(gwas_loci)
    shell:
        """
        module load gsl/2.7 && module load R/4.1.0;
        /software/R-4.1.0-el8-x86_64/bin/Rscript scripts/hyprcoloc_gwas_allTissues.R {params.region} {wildcards.gwas_loci} &> {log}
        """
    

    
rule getClusterAnnotations:
    input:
        expand("/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/{{tissue}}/separateNoise/leafcutter.qqnorm_{chrom}.gz",
        chrom = ['chr' + str(i) for i in range(1, 23)])
    output:
        'results/data/cluster_annotations/tmp/{tissue}.txt'
    log:
        'logs/collect_cluster_annots/{tissue}.log'
    resources:
        mem_mb = 12000
    wildcard_constraints:
        tissue = '|'.join(tissues)
    shell:
        """
        for I in {input};
        do
          (zcat $I  | tail -n+2 - | awk -F'\\t' '{{print $4}}' - | awk -F':' '{{print $0, $4, $5}}' OFS='\\t' - >> {output}) &>> {log}
        done
        """
        
rule collectClusterAnnotations:
    input:
        expand('results/data/cluster_annotations/tmp/{tissue}.txt', tissue=tissues)
        
        
        
        
rule HyprcolocGWASLociPairwise:
    output:
        rds = 'results/coloc/hyprcoloc_results/pairwise_rds/{gwas_loci}.rds',
        tsv = 'results/coloc/hyprcoloc_results/tables/temp_pairwise/{gwas_loci}.tsv'
    resources:
        mem_mb = 58000
    log:
        "/scratch/midway2/cnajar/logs/hyprcoloc/{gwas_loci}.pairwise.log"
    params:
        region = getGWAS_region #_hyprcoloc
    wildcard_constraints:
        gwas_loci = '|'.join(gwas_loci)
    shell:
        """
        module load gsl/2.7 && module load R/4.1.0;
        /software/R-4.1.0-el8-x86_64/bin/Rscript scripts/hyprcoloc_pairwise.R {params.region} {wildcards.gwas_loci} &> {log}
        """
        
rule collect_hyprcoloc_pairwise:
    input:
        expand('results/coloc/hyprcoloc_results/pairwise_rds/{gwas_loci}.rds', gwas_loci = gwas_loci)
        
rule collect_hyprcoloc_allTissues:
    input:
        expand('results/coloc/hyprcoloc_results_allTissues/rds/{gwas_loci}.rds', gwas_loci = gwas_loci)