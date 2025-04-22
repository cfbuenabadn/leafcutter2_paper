
rule liftover_Gwas_stats:
    """
    Convert hg19 bed of summary stats to hg38_summarystats. For convenience with downstream rules, make sure the input bed has Pvalues in column4. Other summary stats can be in later columns
    """
    input:
        bed = "resources/gwas/temp/GRCh37/{trait}.bed",
        chain = "resources/Annotations/hg19ToHg38.over.chain.gz"
    output:
        "resources/gwas/temp/GRCh38/{trait}.bed"
    shadow: "shallow"
    wildcard_constraints:
        trait = '|'.join(list(gwas_df.loc[(gwas_df.source == 'Leafcutter2') & (gwas_df.assembly == 'GRCh37')].trait))
    log:
        "logs/crossmap.{trait}.log"
    shell:
        """
        (CrossMap.py bed {input.chain} {input.bed} {output}) &> {log}
        """
        


rule StandardizeSummaryStatistics:
    input:
        "resources/gwas/temp/GRCh38/{trait}.bed"
    output:
        "resources/gwas/temp/GRCh38/{trait}.beta_se.bed"
    log:
        'logs/gwas/standardize.summary_stats.{trait}.log'
    resources:
        mem_mb = 36000
    params:
        'beta_se'
    wildcard_constraints:
        trait = '|'.join(list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))
    shell:
        """
        Rscript scripts/StandardizeGwasStats.R {input} {params} {output} &> {log}
        """

#rule standardize_new_gwas:
#    input:
#        expand("resources/gwas/temp/GRCh38/{trait}.beta_se.bed", trait = list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))




rule SortCompressAndIndex:
    input:
        "resources/gwas/temp/GRCh38/{trait}.beta_se.bed"
    output:
        bed = temp("resources/gwas/hg38_summary_stats/{trait}.bed"),
        bedgz = "resources/gwas/hg38_summary_stats/{trait}.bed.gz",
        tbi = "resources/gwas/hg38_summary_stats/{trait}.bed.gz.tbi"
    log:
        'logs/gwas/standardize.sort_compress_index.{trait}.log'
    resources:
        mem_mb = 58000
    wildcard_constraints:
        trait = '|'.join(list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))
    shell:
        """
        (head -n1 {input} > {output.bed}) &> {log};
        (tail -n+2 {input} | sort -k 1,1 -k2,2n - >> {output.bed}) &> {log};
        (bgzip -c {output.bed} > {output.bedgz}) &> {log};
        (tabix -p bed {output.bedgz}) &> {log};
        """

rule standardize_new_gwas:
    input:
        expand("resources/gwas/hg38_summary_stats/{trait}.bed.gz", trait = list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))



def GetWindowSize(wildcards):
    if wildcards.WindowSize == 'Window1M':
        return 500000
    elif wildcards.WindowSize == 'Window100K':
        return 50000
    elif wildcards.WindowSize == 'Window10K':
        return 5000

rule GetGWAS_LeadSnpWindows:
    """
    output bed of 1MB window surrounding lead genome-wide signficant autosomal
    snps for each gwas. Exclude blacklistregions (ie MHC). This rule only works
    for summary stats in bed format with Pvalue in column4. Exit status 1 if
    output is empty
    """
    input:
        summarystats = "resources/gwas/hg38_summary_stats/{trait}.bed.gz",
        blacklistregions = "resources/Annotations/MHC.hg38.bed",
        chromsizes = "resources/Annotations/GRCh38.primary_assembly.genome.fa.fai"
    output:
        signif_loci = "resources/gwas/leadSnps/{trait}.bed"
    params:
        PvalThreshold = "5e-8"
    log:
        "logs/GetGWAS_LeadSnpWindows_NonGwasCatalog/{trait}.log"
    wildcard_constraints:
        trait = '|'.join(list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))
    resources:
        mem_mb = 58000
    shell:
        """
        (python scripts/GetGWASLeadVariantWindowsFromBed.py {input.summarystats} /dev/stdout {params.PvalThreshold} | awk -F'\\t' -v OFS='\\t' '$1~/chr[0-9]+/ {{print $1, $2, $2, $1"_"$2"_N_N_{wildcards.trait}" }}' | bedtools slop -i - -g {input.chromsizes} -b 500000 | bedtools sort -i - | bedtools intersect -a - -b {input.blacklistregions} -wa -sorted -v > {output} ) &> {log}
        [[ -s {output.signif_loci} ]]
        """

rule collect_grch38:
    input:
        expand("resources/gwas/leadSnps/{trait}.bed", trait = list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))



########################################

rule ConcatGwasLeadSnpWindows:
    input:
        expand(
            "resources/gwas/leadSnps/{trait}.bed", trait=gwas_traits
        ),
    output:
        "resources/gwas/LeadSnpWindows.bed"
    log:
        "logs/ConcatGwasLeadSnpWindows.log"
    shell:
        """
        cat {input} | bedtools sort -i - > {output}
        """


rule ChangeVarID:
    input:
        "resources/gwas/hg38_summary_stats/{trait}.bed.gz"
    output:
        bed = "resources/gwas/hg38_summary_stats/{trait}.var_id_renamed.bed.gz",
        tbi = "resources/gwas/hg38_summary_stats/{trait}.var_id_renamed.bed.gz.tbi",
    wildcard_constraints:
        trait = '|'.join(list(gwas_df.loc[(gwas_df.source == 'Leafcutter2')].trait))
    resources:
        mem_mb = 58000
    log: 'logs/change_var_name/{trait}.log'
    shell:
        """
        (zcat {input} | awk 'BEGIN{{OFS="\\t"}} NR==1 {{print $0; next}} {{split($6, id_parts, "_"); $6 = id_parts[1] ":" id_parts[2] ":" id_parts[4] ":" id_parts[3]; print $0 }}' - | bgzip -c > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &> {log}
        """
    

def GetSummaryStatsInput(wildcards):
    if wildcards.trait == 'Visceral_adipose_tissue_measurement':
        return "resources/gwas/hg38_summary_stats/{trait}.bed.gz"
    else:
        return "resources/gwas/hg38_summary_stats/{trait}.var_id_renamed.bed.gz"
        
def GetSummaryStatsTbiInput(wildcards):
    if wildcards.trait == 'Visceral_adipose_tissue_measurement':
        return "resources/gwas/hg38_summary_stats/{trait}.bed.gz.tbi"
    else:
        return "resources/gwas/hg38_summary_stats/{trait}.var_id_renamed.bed.gz.tbi"

rule GwasBedStatsAtWindows:
    """
    For coloc, gather summary stats in window centered on lead SNPs
    """
    input:
        signif_loci = "resources/gwas/leadSnps/{trait}.bed",
        bed = GetSummaryStatsInput, #"resources/gwas/hg38_summary_stats/{trait}.bed.gz",
        tbi = GetSummaryStatsTbiInput #"resources/gwas/hg38_summary_stats/{trait}.bed.gz.tbi",
    log:
        "logs/GwasBedStatsAtWindows/{trait}.log"
    output:
        stats = "resources/gwas/StatsForColoc/{trait}.standardized.txt.gz"
    wildcard_constraints:
        accession = '|'.join(gwas_new_traits)
    shell:
        """
        (tabix -h -R {input.signif_loci} {input.bed} | sort -k1,1 -k2,2n | bedtools intersect -sorted -a - -b {input.signif_loci} -wo | awk -F'\\t' -v OFS='\\t'  '{{print $12, $1, $2, $7, $8, $6}}'  | awk -F'[:\\t]' 'BEGIN{{print "loci\\tchrom\\tstart\\tbeta\\tSE\\tA1\\tA2"}} {{print $1, $2, $3, $4, $5, $8, $9}}' OFS='\\t' | gzip - > {output.stats} ) &> {log}
        """

rule CollectLeadSnps:
    input:
        expand("resources/gwas/StatsForColoc/{trait}.standardized.txt.gz", trait=gwas_traits),
        "resources/gwas/LeadSnpWindows.bed"


rule InstallHyprcoloc:
    """
    hyprcoloc r package is not on conda. This rule installs it on the conda environment. Here is a command to recreate the conda environment with dependencies (without hyprcoloc)
    mamba create --name r_hyprcoloc -c r r-rmpfr r-iterpc r-tidyverse r-devtools r-pheatmap r-rcppeigen r-essentials
    For reasons I don't understand, conda won't export this environment to yaml
    with `conda export`. So I manually created an environment with the command
    above, then ran `conda list -e` and manually indented lines to conform to
    yaml to create the conda-compatible yaml file specified.
    """
    output:
        touch("hyprcoloc/hyprcoloc_installed.touchfile")
    log:
        "logs/InstallHyprcoloc.log"
    resources:
        mem_mb = 58000
    shell:
        """
        module unload R/4.1.0 && module unload openblas && module unload gsl;
        module load gsl/2.5 && module load openblas/0.3.13 && module load R/4.1.0;
        /software/R-4.1.0-el8-x86_64/bin/Rscript -e 'remotes::install_github("jrs95/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = FALSE, dependencies=F); install.packages("R.utils", repos = "http://cran.us.r-project.org")' &> {log}
        """