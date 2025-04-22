#!/usr/bin/env Rscript

suppressMessages(library(glue))

if (interactive()) {
    args = scan(what = character(), 
                text = "plotdata/ds_v_dge/Liver_v_Lung_data.rds 
                        plotdata/ds_v_dge/Brain-Cortex_v_Lung_data.rds"
                )
    rds_files <- args
    FDR = c(ds = 1e-3, dge = .05)
    dPSI = .2
} else {
    args = commandArgs(trailingOnly = TRUE)
    if (length(args) != 5) {
        stop("Usage: prepGTExHeatmapData.R < contrast_dir > < out_file > < FDR_ds > < FDR_dge > < dPSI >")
    }
    contrast_dir <- args[1]
    rds_files  <- dir(contrast_dir, '.+_v_.+.rds', full.names = TRUE)
    out_file  <- args[2]
    FDR = c(ds = as.numeric(args[3]), dge = as.numeric(args[4]))
    dPSI = as.numeric(args[5])
}

print(glue("# {Sys.time()} Prepare data input for GTEx UP PSI heatmap:"))
print(glue("With FDR_ds = {FDR['ds']}, FDR_dge = {FDR['dge']}, dPSI = {dPSI}"))


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


#--------------------------------------------------------------
# NOTE: Functions
#--------------------------------------------------------------

GetTopUpIntrons <- function(data, FDR_ds, FDR_dge, dPSI) {
    #' Select unproductive introns based on the following criteria:
    #' 1. ds FDR'
    #' 2. dge FDR
    #' 3. deltaPSI'
    #' resturs joined datatable of introns where deltapsi * log2FoldChange < 0

    dge <- data$dge
    ds <- data$ds

    # if a cluster has multiple UP introns, only select the best 1
    ds <- ds[itype == 'UP' & ctype == 'PR,UP'][, rk := rank(-abs(deltapsi), ties.method = "first"), by = cluster][rk ==1][, rk := NULL][]
    ds <- ds[`p.adjust` < FDR_ds & abs(deltapsi) > dPSI,]

    dge <- dge[padj < FDR_dge,]
    
    ds_excl_cols <- c('itype', 'ctype', 'df', 'p', 'p.adjust','logef', 'loglr', 'status')
    dge_excl_cols <- c('baseMean', 'lfcSE', 'stat', 'pvalue', 'padj')
    chosen <- inner_join(
        x = ds[, -ds_excl_cols, with = FALSE],
        y = dge[, -dge_excl_cols, with = FALSE],
        by = "gene_id",
        suffix = c("_ds", "_dge")
      ) %>%
      .[deltapsi * log2FoldChange < 0, ]

    return(chosen)
}

GetHeatmapMatrix <- function(dt, clusters, contrast) {
  tissues <- str_split(contrast, "_v_") %>% unlist()
  keep1 <- c("intron", "cluster", tissues[1], "itype", "ctype")
  keep2 <- c("intron", "cluster", tissues[2], "itype", "ctype")
  out <- list(
    dt[cluster %in% clusters, keep1, with = FALSE],
    dt[cluster %in% clusters, keep2, with = FALSE]
  )
  names(out) <- tissues

  return(out)
}

sumPSI <- function(dt, cn) {
  dt <- rename(dt, psi = {{ cn }})
  dt <- dt[itype == "UP"][, .(psi = sum(psi)), by = cluster] # sum unprod PSI by cluster
  names(dt) <- c("cluster", cn)
  return(as.data.table(dt))
}

#--------------------------------------------------------------
# NOTE: Load data
#--------------------------------------------------------------

print(glue("# {Sys.time()} Load ds and dge results:"))
contrast_ls <- basename(rds_files) %>%
    str_remove_all('_data\\.rds')
tissues <- str_split(contrast_ls, "_v_") %>% unlist %>% unique
names(rds_files) <- contrast_ls
data <- map(rds_files, readRDS)
print(glue("# {Sys.time()} Loaded {length(data)} datasets(contrasts)"))

print(glue("# {Sys.time()} Get a union set of clusters with UP introns:"))
# UNIONED set of clusters across all tissues, where
# each cluster has at least 1 unproductive intron passing selection criteria
chosen_ls <- map(data, GetTopUpIntrons, FDR_ds = FDR['ds'], FDR_dge = FDR['dge'], dPSI = dPSI)
chosen_clusters <- map(chosen_ls, ~ .x[, cluster]) %>% reduce(union)
print(glue("# {Sys.time()} Collected {length(chosen_clusters)} clusters (unioned) with UP introns"))

# COMMON set of clusters across all tissues, where
# each cluster has been selected in the unioned set
print(glue("# {Sys.time()} Get a common set of clusters with UP introns:"))
plot_clusters <- map(data, ~ .x$ds[, cluster]) %>% # get clusters for each pair
  map(., ~ intersect(.x, chosen_clusters)) %>% # intersect with chosen clusters
  reduce(intersect) # chosen clusters that are common across all tissue pairs
print(glue("# {Sys.time()} Collected {length(plot_clusters)} clusters (common) with UP introns"))

# ds result for common set of colusters
print(glue("# {Sys.time()} Get Heatmap data for common clusters:"))
plotdata <- imap(data, ~ GetHeatmapMatrix(.x$ds, plot_clusters, .y)) %>% 
    unlist(recursive = FALSE)

# get the names of stored plot datatable without duplicating dataset
plotdata.names <- names(plotdata) %>%
  str_split("\\.", simplify = T) %>%
  as.data.table() %>%
  .[, .(V1, rk = rank(V1)), by = V2] %>% # after split, V1=contrast, V2=tissue
  .[rk == 1] %>% # since a tissue can be used in multiple contrasts, only keep 1
  .[, .(V1, nm = paste(V1, ".", V2, sep = ""))] %>%
  .[, nm]

# keep non-duplicated datatable
plotdata <- plotdata[plotdata.names]
names(plotdata) <- str_split(names(plotdata), "\\.") %>% map_chr(~.[2])

# Summarize Unproductive PSI by cluster
plotdata <- imap(plotdata, ~sumPSI(.x, .y))

# ensure same order of clusters
plot_clusters.dt  <- data.table(cluster = naturalsort::naturalsort(plot_clusters))
plotdata <- map(
                plotdata,
                \(df) left_join(plot_clusters.dt, df, by = 'cluster')
)

print(glue("# {Sys.time()} Prepared heatmap data for {length(plotdata)} tissues"))
out_rds = list(
           unioned_chosen_clus = chosen_clusters,
           shared_chosen_clus = plot_clusters,
           shared_UP_psi = plotdata
)

print(glue("# {Sys.time()} Save data to {out_file}"))
write_rds(out_rds, out_file)


