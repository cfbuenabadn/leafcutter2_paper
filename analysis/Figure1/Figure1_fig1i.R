#!/usr/bin/env Rscript
# Panel 1i: unproductive splicing against intron length and against gene
# expression, in the Geuvadis EUR panel of 373 lymphoblastoid cell lines.
#
# Ported from SpliFi/analysis/2024-01-20-unprod-splicing-vs-gene-features.ipynb.
# Writes two tidy tables into figure_data/; the Python notebook draws the ECDFs
# from them, so this only has to run when the underlying data changes.
#
#   fig1i_intron_length.tsv   one row per cluster: mean unproductive ratio and
#                             the length of its most-used productive intron
#   fig1i_expression.tsv      one row per gene: median unproductive ratio and
#                             mean RPKM across the 373 samples
suppressMessages({library(data.table); library(dplyr); library(stringr)})

SPLIFI <- '/project/yangili1/cfbuenabadn/SpliFi'
ANALYSIS <- file.path(SPLIFI, 'analysis')
# The notebook read gene coordinates from a CSV export of this GTF that no
# longer exists. gene_len is just end - start of protein_coding `gene` features,
# so the GTF itself gives identical values.
GTF <- '/project/yangili1/cdai/annotations/hg38/gencode.v43.primary_assembly.annotation.gtf.gz'
OUT <- 'figure_data'
dir.create(OUT, showWarnings = FALSE)

cts <- readRDS(file.path(SPLIFI, 'data/ExtractFractions/Geuvadis/EUR.numerators_constcounts.noise_by_intron.rds'))
stored <- readRDS(file.path(ANALYSIS, '2024-01-20-Geuvadis.EUR.storedcomputes.rds'))
datcols <- names(cts)[grepl('^(HG|NA)[0-9]', names(cts))]
cat('samples:', length(datcols), '\n')

# ---- intron length -------------------------------------------------------- #
# Clusters carrying both productive and unproductive introns; within each, the
# most-used productive intron supplies the length.
dt <- copy(cts)[str_detect(clu_type, 'PR\\,UP')]
dt[, intron_length := as.numeric(str_split(chrom, ':', simplify = TRUE)[, 3]) -
                      as.numeric(str_split(chrom, ':', simplify = TRUE)[, 2])]
dt[, meanCount := rowMeans(.SD), .SDcols = datcols]
pr <- dt[intron_type == 'PR'][order(-meanCount), .SD[1], by = clu][, .(clu, intron_length)]

ratio_clu <- as.data.table(stored$unprod_ratio_by_clu)
ratio_clu[, meanUnprod := rowMeans(.SD, na.rm = TRUE), .SDcols = datcols]
len <- merge(pr, ratio_clu[, .(clu, meanUnprod)], by = 'clu')
len <- len[is.finite(meanUnprod) & is.finite(intron_length)]
len[, intron_length_bin := cut(intron_length,
      breaks = c(0, 1000, 5000, 20000, 50000, Inf),
      labels = c('<1kb', '1-5kb', '5-20kb', '20-50kb', '>50kb'), right = FALSE)]
ct <- cor.test(len$meanUnprod, len$intron_length, method = 'spearman')
cat(sprintf('intron length: n=%d  rho=%.3f  P=%.3g\n', nrow(len), ct$estimate, ct$p.value))
write.table(len, file.path(OUT, 'fig1i_intron_length.tsv'), sep = '\t',
            quote = FALSE, row.names = FALSE)

# ---- expression ----------------------------------------------------------- #
# RPKM per gene: total junction reads over the genomic span of the gene.
gtf <- fread(cmd = paste('zcat', GTF, '| awk -F\'\\t\' \'$3=="gene"\''),
             header = FALSE, sep = '\t',
             select = c(4, 5, 9), col.names = c('start', 'end', 'attr'))
gtf <- gtf[grepl('gene_type "protein_coding"', attr)]
gtf[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
gtf[, gene_len := end - start]
gene_lens <- gtf[, .(gene_len = max(gene_len)), by = gene_name]

juncClass <- fread(file.path(SPLIFI, 'code/results/pheno/noisy/Geuvadis/EUR/wConst_junction_classifications.txt'))
jmap <- unique(juncClass[, .(Intron_coord, Gene_name)])
cts2 <- copy(cts)
cts2[, Intron_coord := paste0(str_split(chrom, ':', simplify = TRUE)[, 1], ':',
                              str_split(chrom, ':', simplify = TRUE)[, 2], '-',
                              str_split(chrom, ':', simplify = TRUE)[, 3])]
cts2 <- merge(cts2, jmap, by = 'Intron_coord', allow.cartesian = FALSE)
total_by_gene <- cts2[, lapply(.SD, sum), by = Gene_name, .SDcols = datcols]
total_by_gene <- merge(total_by_gene, gene_lens, by.x = 'Gene_name', by.y = 'gene_name')

rpkm <- as.data.table(lapply(total_by_gene[, ..datcols],
          function(x) 1e9 * as.double(x) / (sum(as.double(x)) * as.double(total_by_gene$gene_len))))
meanRPKM <- data.table(Gene_name = total_by_gene$Gene_name,
                       meanrpkm = rowMeans(as.matrix(rpkm), na.rm = TRUE))

ratio_gene <- as.data.table(stored$unprod_ratio_by_gene)
ratio_gene[, medUPratio := matrixStats::rowMedians(as.matrix(.SD), na.rm = TRUE), .SDcols = datcols]
expr <- merge(ratio_gene[, .(Gene_name, medUPratio)], meanRPKM, by = 'Gene_name')
expr <- expr[medUPratio > 0 & is.finite(meanrpkm)]
qs <- quantile(expr$meanrpkm, seq(0, 1, .2), na.rm = TRUE)
expr[, rpkm_bin := cut(meanrpkm, breaks = c(0, qs[2:5], Inf),
      labels = c('Q1 - lowly expressed', 'Q2', 'Q3', 'Q4', 'Q5 - highly expressed'))]
ct2 <- cor.test(expr$medUPratio, expr$meanrpkm, method = 'spearman')
cat(sprintf('expression  : n=%d  rho=%.3f  P=%.3g\n', nrow(expr), ct2$estimate, ct2$p.value))
write.table(expr, file.path(OUT, 'fig1i_expression.tsv'), sep = '\t',
            quote = FALSE, row.names = FALSE)
