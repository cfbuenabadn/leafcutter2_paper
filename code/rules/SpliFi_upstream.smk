## =============================================================================
##  UPSTREAM RULES FROM THE SpliFi PROJECT -- COMMENTED OUT ON PURPOSE
## =============================================================================
##
##  Every rule below is a verbatim copy of a rule in
##      /project/yangili1/cfbuenabadn/SpliFi/code/workflow/rules/
##  which is where the GTEx splicing phenotypes, expression phenotypes and
##  genotypes used by this paper were actually produced. They are reproduced
##  here, commented out, so that the provenance of those files is recorded in
##  this repository and so that the chain could be re-run here if the SpliFi
##  results directory ever became unavailable.
##
##  NOTHING HERE IS ACTIVE. Rules in this repository read the SpliFi outputs by
##  absolute path (/project/yangili1/cfbuenabadn/SpliFi/code/results/...). Each
##  of those rules also carries a commented-out `input:` line pointing at the
##  local path that the rules below would produce, so switching over is a matter
##  of uncommenting this file, its `include:` in the Snakefile, and those lines.
##
##  ---------------------------------------------------------------------------
##  UNDECLARED OUTPUTS -- the reason this file is worth having
##  ---------------------------------------------------------------------------
##  Three of these rules create files inside their `shell:` block that are never
##  declared in `output:`. Snakemake therefore cannot see them, will not rebuild
##  them if deleted, and reports the parent rule as up to date. This repository
##  consumes exactly those files. Where that happens the copy below has the
##  missing outputs ADDED (marked "# ADDED"), which is the one deliberate
##  difference from the SpliFi original:
##
##    Leafcutter2Gtex   declares  leafcutter_perind.counts.noise_by_intron.gz
##                      also writes  leafcutter_perind_numers.counts.noise_by_intron.gz  <- we read this
##                                   leafcutter_perind_numers.counts.gz                  <- and this
##    PreparePhenoBed   declares  done, leafcutter_names.txt
##                      also writes  separateNoise/leafcutter.qqnorm_chr{1..22}.gz       <- we read this
##                                   separateNoise/leafcutter.phen_chr{1..22}.gz         <- and this
##    sortGTExPhenotype declares  qqnorm.sorted.bed.gz
##                      also writes  qqnorm.sorted.chr{1..22}.bed.gz                     <- we read this
##
##  GTExPhenoPC and MakeCovarianceMatrix declare their outputs correctly.
##
##  ---------------------------------------------------------------------------
##  THE CHAIN, as consumed by this repository
##  ---------------------------------------------------------------------------
##    GTEx junctions
##      -> ClusterIntronsGtex            leafcutter_refined_noisy
##      -> Leafcutter2Gtex               leafcutter_perind*.counts.noise_by_intron.gz
##      -> PreparePhenoBed               separateNoise/leafcutter.qqnorm_{chrom}.gz
##      -> MakeCovarianceMatrix          separateNoise/{chrom}_CovMatrix.txt
##                                       [+ GenotypePCA -> geno/{chrom}.pca]
##    GTEx expression
##      -> PreareGTExEQTLPhenotype       eqtl/qqnorm.bed.gz
##      -> sortGTExPhenotype             eqtl/qqnorm.sorted.{chrom}.bed.gz
##      -> GTExPhenoPC                   eqtl/qqnorm.sorted.{chrom}.pca
##    GTEx WGS
##      -> ExtractGenotypeVCF            geno/{tissue}/{chrom}.vcf.gz
##    All-tissue clustering (differential splicing)
##      -> ClusterIntronsGtexAllTissues  ds/GTEx/all49tissues_refined_noisy
##
##  Paths are left exactly as they are in SpliFi and are relative to
##  SpliFi/code/. Run from there, or rewrite the prefixes, if this is ever
##  re-enabled.
## =============================================================================

## ---------------------------------------------------------------------------
## ClusterIntronsGtex   [SpliFi rules/gtex.smk]
## Pre-clusters GTEx introns; produces leafcutter_refined_noisy, the -c input to Leafcutter2Gtex.
## ---------------------------------------------------------------------------
# rule ClusterIntronsGtex:
#     message: '### Make intron clusters using any GTEx samples'
#     input:
#         junc_files_flag = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted/done',
#     output:
#         pooled   = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_pooled',
#         clusters = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_refined_noisy',
#         lowusage = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_lowusage_introns', # intermediate
#         refined  = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_refined' # intermediate
#     params:
#         run_dir    = 'results/pheno/noisy/GTEx/{tissue}',
#         out_prefix = 'leafcutter',
#         junc_files = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted',
#         py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
#     log: 'logs/ClusterIntronsGtex/{tissue}.log'
#     shell:
#         '''
#         python {params.py_script} \
#             -r {params.run_dir} \
#             -o {params.out_prefix} \
#             -j <(realpath {params.junc_files}/*.tsv.gz) &> {log}
#         ls -lah {output.pooled} {output.clusters} {output.lowusage} {output.refined} &>> {log}
#
#         '''

## ---------------------------------------------------------------------------
## Leafcutter2Gtex   [SpliFi rules/gtex.smk]
## Runs leafcutter2 per tissue. Source of the noise_by_intron tables this paper reads.
## ---------------------------------------------------------------------------
# rule Leafcutter2Gtex:
#     message:'### Run leafcutter2 on GTEx samples'
#     input:
#         junc_files_flag = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted/done',
#         junc_files = 'resources/GTEx/juncs/groupped_juncs/{tissue}/converted',
#         intron_clusters = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_refined_noisy',
#     output:
#         perind_noise_by_intron = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_perind.counts.noise_by_intron.gz',
#         # ADDED -- written by the same leafcutter2 call but undeclared upstream:
#         perind_numers_noise_by_intron = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_perind_numers.counts.noise_by_intron.gz',
#         perind_numers = 'results/pheno/noisy/GTEx/{tissue}/leafcutter_perind_numers.counts.gz'
#     params:
#         py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter2_regtools.py',
#         run_dir    = 'results/pheno/noisy/GTEx/{tissue}',
#         out_prefix = 'leafcutter', # note do not include parent dir
#         gtf = config['annotation']['gtf']['v43'],
#         genome = config['genome38'],
#         pre_clustered = '-c results/pheno/noisy/GTEx/{tissue}/leafcutter_refined_noisy',
#         max_juncs = 1000, # maximum number of introns per gene
#         other_params = '-k' # not keeping constitutive introns
#     threads: 1
#     resources: cpu=1, time=2100, mem_mb=25000
#     log: 'logs/Leafcutter2Gtex/{tissue}.log'
#     shell:
#         '''
#         python {params.py_script} \
#             -j <(realpath {input.junc_files}/*.tsv.gz) \
#             -r {params.run_dir} \
#             -o {params.out_prefix} \
#             -A {params.gtf} \
#             -G {params.genome} \
#             --max_juncs {params.max_juncs} \
#             {params.pre_clustered} {params.other_params} &> {log}
#
#         ls {output.perind_noise_by_intron} &>> {log}
#
#
#         '''

## ---------------------------------------------------------------------------
## PreparePhenoBed   [SpliFi rules/qtl.smk]
## Turns the leafcutter2 counts into QTLtools phenotype BEDs, split per chromosome.
## ---------------------------------------------------------------------------
# rule PreparePhenoBed:
#     message: '### Prepare phenotype bed file for qtltools, using *.counts.noise.gz'
#     input: getPreparePhenoBedInput
#     output:
#         flag = touch('results/pheno/noisy/{datasource}/{group}/{phenType}/done'),
#         samples = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter_names.txt',
#         # ADDED -- written by the bash loop in shell: but undeclared upstream:
#         qqnorm = expand('results/pheno/noisy/{{datasource}}/{{group}}/{{phenType}}/leafcutter.qqnorm_{chrom}.gz',
#                         chrom = ['chr' + str(i) for i in range(1, 23)]),
#         phen = expand('results/pheno/noisy/{{datasource}}/{{group}}/{{phenType}}/leafcutter.phen_{chrom}.gz',
#                       chrom = ['chr' + str(i) for i in range(1, 23)])
#     params:
#         pyscript = 'workflow/scripts/preparePheno.py',
#         outPrefix = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter',
#         vcfSamples = getVcfIndiv ,# individual ids in vcf file
#     log: 'logs/PreparePhenoBed_{datasource}_{group}_{phenType}.log'
#     shell:
#         '''
#         python {params.pyscript} {input} --sampleFile {params.vcfSamples} --outPrefix {params.outPrefix} &> {log}
#
#         for i in {{1..22}}; do
#             bgzip -f {params.outPrefix}.phen_chr${{i}} 2>> {log}
#             bgzip -f {params.outPrefix}.qqnorm_chr${{i}} 2>> {log}
#             tabix -f -p bed {params.outPrefix}.qqnorm_chr${{i}}.gz 2>> {log}
#         done
#
#         ls -l {output.samples} &>> {log}
#
#
#         '''
#
#
# def getExtractGenotypeInput(wildcards):
#     if wildcards.datasource == 'GTEx':
#         vcf = config['VCF']['GTEx']['HG38_v7']
#     elif wildcards.datasource == 'Geuvadis':
#         vcf_dir = config['VCF']['Geuvadis']['HG38_1kg_b38']
#         chrom = wildcards.chrom
#         vcf =  f'{vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz'
#     else:
#         print('Error. Invalid vcf file path. Exiting...')
#         exit(0)
#     return vcf
#
# def getExtractGenotypeParams(wildcards):
#     if wildcards.datasource == 'GTEx':
#         min_MAF, chrom = 0.05, wildcards.chrom
#         expr = f'AF >= {min_MAF} '
#     elif wildcards.datasource == 'Geuvadis':
#         min_MAF, max_HWE = 0.05, 1e-3
#         expr = f'AF >= {min_MAF} && HWE < {max_HWE}'
#     else:
#         print('Error. Invalid vcf file path. Exiting...')
#         exit(0)
#     return expr

## ---------------------------------------------------------------------------
## MakeCovarianceMatrix   [SpliFi rules/qtl.smk]
## Builds the sQTL covariate matrix (11 phenotype PCs + 5 genotype PCs).
## ---------------------------------------------------------------------------
# rule MakeCovarianceMatrix:
#     message: '### Make covariance matrix with fixed PCs'
#     input:
#         samples = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter_names.txt',  # for changes
#         GenoPCs = 'results/geno/{datasource}/{group}/{chrom}.pca',
#     output: 'results/pheno/noisy/{datasource}/{group}/{phenType}/{chrom}_CovMatrix.txt'
#     params:
#         n_phenoPCs = 11, # index starts from header row
#         n_genoPCs = 5, # index starts from header row
#         PhenoPCs = 'results/pheno/noisy/{datasource}/{group}/{phenType}/leafcutter.PCs',
#         fake = "fake",
#     log: 'logs/MakeCovarianceMatrix_{datasource}_{group}_{phenType}_{chrom}.log'
#     run:
#         fout = open(output[0], 'w')
#         with open(input.samples) as f:
#             samples = f.readlines()
#             samples = [x.strip() for x in samples]
#             fout.write('\t'.join(['id'] + samples) + '\n')
#
#         with open(input.GenoPCs) as f: # append genotype PCs
#             genoHeader = f.readline().strip().split()
#             print("Getting genotype PCs...")
#             if not all(a == b for a,b in zip(samples, genoHeader[1:])):
#                 print(f"Samples in genotype do not match!\nsamples:{samples}\nGenotype:{genoHeader[1:]}")
#                 exit("Samples genotype do not match! Exiting...")
#             i = 1
#             for ln in f.readlines():
#                 if i > params.n_genoPCs:
#                     break
#                 ln = ln.strip().split()
#                 ln[0] = f'genoPC:{ln[0]}'
#                 fout.write('\t'.join(ln) + '\n')
#                 i += 1
#
#         with open(params.PhenoPCs) as f: # append phenotype PCs
#             print("Getting phenotype PCs...")
#             phenoHeader = f.readline().strip().split()
#             if not all(a == b for a,b in zip(samples, phenoHeader[1:])):
#                 exit("Samples in phenotype do not match! Exiting...")
#             i = 1
#             for ln in f.readlines():
#                 if i > params.n_phenoPCs:
#                     break
#                 ln = ln.strip().split()
#                 ln[0] = f'phenoPC:{ln[0]}'
#                 fout.write('\t'.join(ln) + '\n')
#                 i += 1
#
#         print(f"Done.\nWrote to covariance matrix file: {output[0]}")
#         fout.close()

## ---------------------------------------------------------------------------
## ExtractGenotypeVCF   [SpliFi rules/qtl.smk]
## Per-tissue, per-chromosome genotypes from the GTEx WGS VCF.
## ---------------------------------------------------------------------------
# rule ExtractGenotypeVCF:
#     '''
#         Procedures:
#             1.  extract vcf per chromosome for samples in phenotype
#             2.  tabix vcf
#     '''
#     message: '### Extract genotype vcf for phenotype samples'
#     input:
#         vcf = getExtractGenotypeInput,
#         sample_file = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter_names.txt'
#     output:
#         vcf = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz',
#         tbi = 'results/geno/{datasource}/{group}/{chrom}.vcf.gz.tbi'
#     params:
#         expr = getExtractGenotypeParams
#     log: 'logs/ExtractGenotypeVCF_{datasource}_{group}_{chrom}.log'
#     threads: 4
#     resources: time=2000, mem_mb=15000, cpu=4
#     group: 'geno'
#     shell:
#         '''
#         bcftools view \
#             --threads {threads} \
#             --samples-file {input.sample_file} \
#             -i "{params.expr}" \
#             -Oz -o {output.vcf} \
#             {input.vcf} \
#             {wildcards.chrom} &> {log}
#
#         bcftools index --threads {threads} --tbi {output.vcf} &>> {log}
#         '''
#
# # sometimes qtltools fail, module unload gsl gcc and reload gcc/10.2.0 and gsl/2.2.1 works

## ---------------------------------------------------------------------------
## GenotypePCA   [SpliFi rules/qtl.smk]
## Genotype PCs, consumed by MakeCovarianceMatrix.
## ---------------------------------------------------------------------------
# rule GenotypePCA:
#     input:  'results/geno/{datasource}/{group}/{chrom}.vcf.gz'
#     output: 'results/geno/{datasource}/{group}/{chrom}.pca'
#     params:
#         out_prefix = 'results/geno/{datasource}/{group}/{chrom}'
#     log: 'logs/GenotypePCA_{datasource}_{group}_{chrom}.log'
#     group: 'geno'
#     shell:
#         '''
#             module unload gsl && module load gsl/2.5
#             QTLtools pca \
#             --seed 123 \
#             --maf 0.05 \
#             --vcf {input}  \
#             --out {params.out_prefix} \
#             --center \
#             --scale &> {log}
#         '''

## ---------------------------------------------------------------------------
## PreareGTExEQTLPhenotype   [SpliFi rules/eqtl.smk]
## GTEx expression phenotype BED (note the typo in the rule name is upstream).
## ---------------------------------------------------------------------------
# rule PreareGTExEQTLPhenotype:
#     input:
#         cnt = 'resources/{datasource}/expression/{group}_gene_reads.tsv.gz',
#         genelist = 'resources/{datasource}/ExpressedGeneList.txt',
#         samples = 'results/pheno/noisy/{datasource}/{group}/separateNoise/leafcutter_names.txt' # intersect of gtex pheno and gtex vcf avail. samples
#     output:
#         cpm = 'results/eqtl/{datasource}/{group}/cpm.bed.gz',
#         qq = 'results/eqtl/{datasource}/{group}/qqnorm.bed.gz',
#     params:
#         R_script = 'workflow/scripts/fromBen_PrepareGTExEQTLPheno.R',
#     log: 'logs/PrepareGTExEQTLPheno_{datasource}_{group}.log'
#     shell:
#         '''
#         Rscript {params.R_script} {input.cnt} {input.genelist} {input.samples} {output.cpm} {output.qq}
#
#         '''

## ---------------------------------------------------------------------------
## sortGTExPhenotype   [SpliFi rules/eqtl.smk]
## Sorts and indexes the expression BED, and splits it per chromosome.
## ---------------------------------------------------------------------------
# rule sortGTExPhenotype:
#     input: 'results/eqtl/{datasource}/{group}/qqnorm.bed.gz',
#     output:
#         bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz',
#         tbi = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.bed.gz.tbi',
#         # ADDED -- written by the bash loop in shell: but undeclared upstream:
#         per_chrom = expand('results/eqtl/{{datasource}}/{{group}}/qqnorm.sorted.{chrom}.bed.gz',
#                            chrom = ['chr' + str(i) for i in range(1, 23)]),
#     params:
#         bed_prefix = 'results/eqtl/{datasource}/{group}/qqnorm.sorted',
#     shell:
#         '''
#         bedtools sort -header -i {input} | bgzip -c > {output.bed}
#         tabix -p bed {output.bed}
#
#         # split file by chromosome
#         for i in {{1..22}}; do
#             tabix --print-header {output.bed} chr$i| bgzip -c > {params.bed_prefix}.chr$i.bed.gz
#             tabix -p bed {params.bed_prefix}.chr$i.bed.gz
#         done
#
#         '''

## ---------------------------------------------------------------------------
## GTExPhenoPC   [SpliFi rules/eqtl.smk]
## Expression phenotype PCs, consumed by MakeCovMatrixEqtl.
## ---------------------------------------------------------------------------
# rule GTExPhenoPC:
#     input:
#         bed = 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.bed.gz',
#     output: 'results/eqtl/{datasource}/{group}/qqnorm.sorted.{chrom}.pca'
#     params:
#         R_script = 'workflow/scripts/fromBen_PermuteAndPCA.R',
#     shell:
#         '''
#         Rscript {params.R_script} {input} {output}
#         '''

## ---------------------------------------------------------------------------
## ClusterIntronsGtexAllTissues   [SpliFi rules/ds-dge.smk]
## All-49-tissue clustering; produces all49tissues_refined_noisy used by the DS rules here.
## ---------------------------------------------------------------------------
# rule ClusterIntronsGtexAllTissues:
#     message: '### Make intron clusters using all (49) GTEx tissues'
#     output:
#         pooled   = 'results/ds/GTEx/all49tissues_pooled',
#         clusters = 'results/ds/GTEx/all49tissues_refined_noisy',
#         lowusage = 'results/ds/GTEx/all49tissues_lowusage_introns', # intermediate
#         refined  = 'results/ds/GTEx/all49tissues_refined' # intermediate
#     params:
#         run_dir    = 'results/ds/GTEx',
#         out_prefix = 'all49tissues',
#         junc_files = 'resources/GTEx/juncs/all49tissues',
#         py_script  = 'workflow/submodules/leafcutter2/scripts/leafcutter_make_clusters.py'
#     log: 'logs/ClusterIntronsGtexAllTissues/all49tissues.log'
#     resources: cpu=1, time=2100, mem_mb=35000
#     shell:
#         '''
#         python {params.py_script} \
#             -r {params.run_dir} \
#             -o {params.out_prefix} \
#             -j <(ls {params.junc_files}/*.tsv.gz) &> {log}
#
#         ls -lah {output.pooled} {output.clusters} {output.lowusage} {output.refined} &>> {log}
#         '''

