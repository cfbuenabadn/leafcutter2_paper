library(tidyverse)
library(dplyr)
library(bedr)
library(readr)
library(glue)
library(qvalue)

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(data.table))
suppressMessages(library(bedtoolsr))
suppressMessages(library(GenomicRanges))

if (interactive()) {
  tissue <- "Liver"
  gtf.f <- "/project2/yangili1/cdai/annotations/hg38/gencode.v26.GRCh38.genes.csv"
  out.f <- '/project/yangili1/cfbuenabadn/tmp.tsv'


} else {
  args = commandArgs(trailingOnly = TRUE)
  tissue <- args[1]
  gtf.f <- args[2]
  out.f <- args[3]
}

# Function copied from Chao, slightly modified to fit QTLTools output

labelIntron_sqtl <- function(df) {
    df <- df[, 
          # get intron id (removing label)
        .(intron = str_sub(phe_id, 1, -4),
          # extract cluster id
          cluster = str_split(phe_id, ":") %>% map_chr(~paste(.x[1], .x[4], sep=":")),
          # extract intron type (label)
          itype = str_sub(phe_id, -2, -1)
    )]
    # get cluster type
    df <- df[, ctype := paste(sort(unique(itype)), sep="", collapse=","), by = cluster][]

    return(df)
}

annotate_perm_sqtls <- function(tissue, gtf, chrom){
    pfx <- '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/sqtl/GTEx/'
    sqtl.file <- glue::glue(pfx, tissue, '/cis_100000/perm/', chrom, '.temp.txt')
    sqtl <- fread(sqtl.file)
    
    sqtl.intron <- labelIntron_sqtl(sqtl[,1])
    sqtl.w.intron <- cbind(sqtl, (sqtl.intron[,'intron']))
    
    coords <- str_split(sqtl.intron$intron, ":", simplify = TRUE)  %>% as.data.table()
    setnames(coords, c("seqname", "start", "end", "cluster"))
    coords[, strand := str_sub(cluster, -1, -1)]
    coords[, cluster := NULL]
    sqtl.intron <- cbind(sqtl.intron, coords)
    sqtl.intron <- makeGRangesFromDataFrame(
      sqtl.intron, 
      keep.extra.columns = TRUE,
      ignore.strand = FALSE,
      starts.in.df.are.0based = TRUE)

    # use overlap to get gene_name and gene_id that intron belongs to
    olaps <- findOverlaps(sqtl.intron, gtf, type="within", select="all", ignore.strand=FALSE)
    
    an.intron <- sqtl.intron[olaps@from] # annotated introns
    mcols(an.intron) <- cbind(mcols(gtf[olaps@to]), mcols(an.intron))


    an.intron <- as.data.table(an.intron)[
      , .(gene_name, gene_id, rk = frank(gene_name)), 
      by = .(seqnames, start, end, strand, intron, cluster, itype, ctype)
      ][, .(gene_name, gene_id, maxrk = max(rk)), 
        by = .(seqnames, start, end, strand, intron, cluster, itype, ctype)
      ][maxrk == 1, -c("maxrk")]

    an.sqtl <- inner_join(an.intron[, -c("seqnames", "start", "end", "strand")], sqtl.w.intron, by = c("intron"))


    return(an.sqtl)
}


#gtf.f <- "/project/yangili1/cfbuenabadn/leafcutter2_paper/code/annotations/gencode.v26.GRCh38.genes.csv"
gtf <- fread(gtf.f) %>% 
  .[feature == 'gene' & gene_type == 'protein_coding',
    .(seqname, start, end, gene_name, gene_id, strand)] %>%
  unique()

gtf <- makeGRangesFromDataFrame(gtf,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE
)

chrom_list <- paste0('chr', 1:22)

an.sqtls.list <- lapply(chrom_list, function(chrom) {
  annotate_perm_sqtls(tissue, gtf, chrom)
})

# Combine results into a single dataframe
an.sqtls <- do.call(rbind, an.sqtls.list)

an.sqtls[,'q'] <- (an.sqtls %>% pull(adj_beta_pval) %>% qvalue(.))$qvalues %>% signif(., 5)

an.sqtls %>% write_tsv(., out.f)

