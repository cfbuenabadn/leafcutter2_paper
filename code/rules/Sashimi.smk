rule BedGraph2BigWig:
    input:
        bed = '/project2/yangili1/GTEx_v8/bedGraph/{tissue}/{sample}.bed.gz',
        chromSizes = 'resources/hg38_w_chrEBV.chrom.sizes'
    output:
        tmp = 'resources/GTEx/BigWig/{tissue}/{sample}.tmp.bed',
        tmp_sort = 'resources/GTEx/BigWig/{tissue}/{sample}.tmp.sorted.bed',
        bw = 'resources/GTEx/BigWig/{tissue}/{sample}.bw'
    threads: 8
    resources: cpu=8, mem_mb=24000, time=1200
    log:
        'logs/bed2bw.{tissue}.{sample}.log'
    shell:
        """
        (bgzip -f -d -c {input.bed} > {output.tmp}) &> {log};
        (bedtools sort -i {output.tmp} > {output.tmp_sort}) &>> {log};
        (scripts/bedGraphToBigWig {output.tmp_sort} {input.chromSizes} {output.bw}) &>> {log}
        """
    
      
      
rule CollectBigWig:
    input:
        expand('resources/GTEx/BigWig/Brain_Anterior_cingulate_cortex_BA24/{sample}.bw', sample = ba24_samples),
        expand('resources/GTEx/BigWig/Brain_Frontal_Cortex_BA9/{sample}.bw', sample = ba9_samples),
        expand('resources/GTEx/BigWig/Brain_Cortex/{sample}.bw', sample = bc_samples),
        expand('resources/GTEx/BigWig/Brain_Putamen_basal_ganglia/{sample}.bw', sample = bputamen_samples),
        expand('resources/GTEx/BigWig/Heart_Atrial_Appendage/{sample}.bw', sample = heart_samples),
        expand('resources/GTEx/BigWig/Lung/{sample}.bw', sample = lung_samples),
        expand('resources/GTEx/BigWig/Skin_Not_Sun_Exposed_Suprapubic/{sample}.bw', sample = skin_samples),
        expand('resources/GTEx/BigWig/Liver/{sample}.bw', sample = liver_samples),
        expand('resources/GTEx/BigWig/Muscle_Skeletal/{sample}.bw', sample = ms_samples),
        expand('resources/GTEx/BigWig/Whole_Blood/{sample}.bw', sample = wb_samples)


rule TabixPhenTables:
    input:
        "/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/{tissue}/separateNoise/leafcutter.phen_{chrom}.gz"
    output:
        bed = "results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.gz",
        tbi = "results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.gz.tbi"
    resources: mem_mb=24000
    log:
        "logs/sort_phenotypes_and_tabix/{tissue}.{chrom}.log"
    shell:
        """
        (zcat {input} | bedtools sort -header -i - | bgzip -c > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log}
        """
        
rule TabixCountTables:
    input:
        "/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/{tissue}/leafcutter_perind_numers.counts.gz"
    output:
        bed = "results/pheno/GTEx/{tissue}/leafcutter.counts.sorted.gz",
        tbi = "results/pheno/GTEx/{tissue}/leafcutter.counts.sorted.gz.tbi"
    log:
        "logs/sort_counts/{tissue}.log"
    resources: mem_mb=24000
    shell:
        """
        (zcat {input} | tail -n+2 - | awk -F':' '{{print $1, $2, $3, $4, $0}}' OFS='\\t' - | awk '{{gsub(/ /, "\\t"); print}}' - | bedtools sort -header -i - | bgzip -c > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log}
        """
        
rule CollectTabixPhenotypes:
    input:
        expand("results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.gz", tissue = tissues, 
        chrom = ['chr' + str(x) for x in range(1, 23)]),
        expand("results/pheno/GTEx/{tissue}/leafcutter.counts.sorted.gz", tissue=tissues)
        
        

rule PrepareBedForSashimiLinks:
    input:
        "results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.gz"
    output:
        temp("results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.for_links.tmp.bed.gz")
    resources: mem_mb=24000
    log:
        "logs/prepare_for_links/{tissue}.{chrom}.log"
    shell:
        """
        python scripts/prepare_bed_for_sashimi_links.py {input} {output} &> {log}
        """
        
rule TabixPhenTablesForLinks:
    input:
        "results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.for_links.tmp.bed.gz"
    output:
        bed = "results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.for_links.bed.gz",
        tbi = "results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.for_links.bed.gz.tbi"
    resources: mem_mb=24000
    log:
        "logs/sort_phenotypes_and_tabix/{tissue}.{chrom}.log"
    shell:
        """
        (zcat {input} | bedtools sort -header -i - | bgzip -c > {output.bed}) &> {log};
        (tabix -p bed {output.bed}) &>> {log}
        """        

rule CollectTabixPhenotypesForLinks:
    input:
        expand("results/pheno/GTEx/{tissue}/leafcutter.phen_{chrom}.sorted.for_links.bed.gz", tissue = tissues, 
        chrom = ['chr' + str(x) for x in range(1, 23)])
        