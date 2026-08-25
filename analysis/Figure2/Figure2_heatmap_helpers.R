# ---------------------------------------------------------------------------
# Data loading and processing for fig2_heatmap.
#
# Created with Figure code cleaner.
# Source notebooks: ../DS_Heatmap2.ipynb and ../DS_Heatmap.ipynb (the plotting
# cell is byte-identical in both; the duplicate was collapsed).
#
# The source notebook builds the matrix twice -- once per unproductive intron
# and once averaged per gene (its commented-out fig2_heatmap_aver.pdf) -- but
# the two give the same figure, so only the per-intron matrix is built here.
# ---------------------------------------------------------------------------

.libPaths(c('/project/yangili1/cfbuenabadn/R/x86_64-pc-linux-gnu-library/4.1',
            '/software/R-4.1.0-el8-x86_64/lib64/R/library'))

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# --- Paths (absolute: this notebook sits one directory deeper than the
# --- originals, where '../code/...' would no longer resolve) ----------------

PSI_TABLE <- paste0('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/',
                    'analysis_files/GTEx.psi.tsv.gz')
RDS_DIR <- paste0('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/',
                  'results/ds_v_dge_confounder/rds_files/')

# Column indices of the 50 GTEx tissues in GTEx.psi.tsv.gz
TISSUE_COLS <- 7:56

# Display names for those 50 columns, in file order (source cell `colnames(X)`).
TISSUE_DISPLAY_NAMES <- c(
    'Artery - Aorta', 'Brain - Cortex', 'Brain - Nucleus accumbens (basal ganglia)', 'Esophagus - Muscularis', 
    'Adipose - Visceral (omentum)', 'Adrenal gland', 'Artery - Coronary', 'Bladder', 'Brain - Amygdala', 
    'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar hemisphere', 
    'Brain - Cerebellum', 'Brain - Hypothalamus', 'Brain - Putamen (basal ganglia)', 'Brain - Spinalcord (cervical c-1)', 
    'Brain - Substantia nigra', 'Cells - Cultured fibroblasts', 'Esophagus - Gastroesophageal junction', 
    'Esophagus - Mucosa', 'Heart - Atrial appendage', 
    'Heart - Left ventricle', 'Kidney - Cortex', 'Lung', 'Minor salivary gland', 'Ovary', 'Pancreas', 'Pituitary', 
    'Prostate', 'Small intestine - Terminal ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 
    'Adipose - Subcutaneous', 'Artery - Tibial', 'Brain - Frontal cortex (BA9)', 'Brain - Hippocampus', 
    'Breast - Mammary tissue', 'Cells - EBV-transformed lymphocytes', 'Colon - Sigmoid', 'Colon - Transverse', 
    'Liver', 'Muscle - Skeletal', 'Nerve - Tibial', 'Skin - Not sun exposed (suprapubic)', 
    'Skin - Sun exposed (lower leg)', 'Whole blood'
)


# GTEx tissue colours, keyed by the display names above (source notebook).
gtex_colors <- list(
  "Adipose - Subcutaneous" = list(
    "tissue_abbrv" = "ADPSBQ", 
    "tissue_color_hex" = "FFA54F", 
    "tissue_color_rgb" = "255,165,79"
  ), 
  "Adipose - Visceral (omentum)" = list(
    "tissue_abbrv" = "ADPVSC", 
    "tissue_color_hex" = "EE9A00", 
    "tissue_color_rgb" = "238,154,0"
  ), 
  "Adrenal gland" = list(
    "tissue_abbrv" = "ADRNLG", 
    "tissue_color_hex" = "8FBC8F", 
    "tissue_color_rgb" = "143,188,143"
  ), 
  "Artery - Aorta" = list(
    "tissue_abbrv" = "ARTAORT", 
    "tissue_color_hex" = "8B1C62", 
    "tissue_color_rgb" = "139,28,98"
  ), 
  "Artery - Coronary" = list(
    "tissue_abbrv" = "ARTCRN", 
    "tissue_color_hex" = "EE6A50", 
    "tissue_color_rgb" = "238,106,80"
  ), 
  "Artery - Femoral" = list(
    "tissue_abbrv" = "ARTFMR", 
    "tissue_color_hex" = "FF4500", 
    "tissue_color_rgb" = "255,69,0"
  ), 
  "Artery - Tibial" = list(
    "tissue_abbrv" = "ARTTBL", 
    "tissue_color_hex" = "FF0000", 
    "tissue_color_rgb" = "255,0,0"
  ), 
  "Bladder" = list(
    "tissue_abbrv" = "BLDDER", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Brain - Amygdala" = list(
    "tissue_abbrv" = "BRNAMY", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Anterior cingulate cortex (BA24)" = list(
    "tissue_abbrv" = "BRNACC", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Caudate (basal ganglia)" = list(
    "tissue_abbrv" = "BRNCDT", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Cerebellar hemisphere" = list(
    "tissue_abbrv" = "BRNCHB", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Cerebellum" = list(
    "tissue_abbrv" = "BRNCHA", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Cortex" = list(
    "tissue_abbrv" = "BRNCTXA", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Frontal cortex (BA9)" = list(
    "tissue_abbrv" = "BRNCTXB", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Hippocampus" = list(
    "tissue_abbrv" = "BRNHPP", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Hypothalamus" = list(
    "tissue_abbrv" = "BRNHPT", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Nucleus accumbens (basal ganglia)" = list(
    "tissue_abbrv" = "BRNNCC", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Putamen (basal ganglia)" = list(
    "tissue_abbrv" = "BRNPTM", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Spinalcord (cervical c-1)" = list(
    "tissue_abbrv" = "BRNSPC", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Brain - Substantia nigra" = list(
    "tissue_abbrv" = "BRNSNG", 
    "tissue_color_hex" = "EEEE00", 
    "tissue_color_rgb" = "238,238,0"
  ), 
  "Breast - Mammary tissue" = list(
    "tissue_abbrv" = "BREAST", 
    "tissue_color_hex" = "00CDCD", 
    "tissue_color_rgb" = "0,205,205"
  ), 
  "Cells - EBV-transformed lymphocytes" = list(
    "tissue_abbrv" = "LCL", 
    "tissue_color_hex" = "EE82EE", 
    "tissue_color_rgb" = "238,130,238"
  ), 
  "Cells - Cultured fibroblasts" = list(
    "tissue_abbrv" = "FIBRBLS", 
    "tissue_color_hex" = "9AC0CD", 
    "tissue_color_rgb" = "154,192,205"
  ), 
  "Cervix - Ectocervix" = list(
    "tissue_abbrv" = "CVXECT", 
    "tissue_color_hex" = "EED5D2", 
    "tissue_color_rgb" = "238,213,210"
  ), 
  "Cervix - Endocervix" = list(
    "tissue_abbrv" = "CVSEND", 
    "tissue_color_hex" = "EED5D2", 
    "tissue_color_rgb" = "238,213,210"
  ), 
  "Colon - Sigmoid" = list(
    "tissue_abbrv" = "CLNSGM", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Colon - Transverse" = list(
    "tissue_abbrv" = "CLNTRN", 
    "tissue_color_hex" = "EEC591", 
    "tissue_color_rgb" = "238,197,145"
  ), 
  "Esophagus - Gastroesophageal junction" = list(
    "tissue_abbrv" = "ESPGEJ", 
    "tissue_color_hex" = "8B7355", 
    "tissue_color_rgb" = "139,115,85"
  ), 
  "Esophagus - Mucosa" = list(
    "tissue_abbrv" = "ESPMCS", 
    "tissue_color_hex" = "8B7355", 
    "tissue_color_rgb" = "139,115,85"
  ), 
  "Esophagus - Muscularis" = list(
    "tissue_abbrv" = "ESPMSL", 
    "tissue_color_hex" = "CDAA7D", 
    "tissue_color_rgb" = "205,170,125"
  ), 
  "Fallopian tube" = list(
    "tissue_abbrv" = "FLLPNT", 
    "tissue_color_hex" = "EED5D2", 
    "tissue_color_rgb" = "238,213,210"
  ), 
  "Heart - Atrial appendage" = list(
    "tissue_abbrv" = "HRTAA", 
    "tissue_color_hex" = "B452CD", 
    "tissue_color_rgb" = "180,82,205"
  ), 
  "Heart - Left ventricle" = list(
    "tissue_abbrv" = "HRTLV", 
    "tissue_color_hex" = "7A378B", 
    "tissue_color_rgb" = "122,55,139"
  ), 
  "Kidney - Cortex" = list(
    "tissue_abbrv" = "KDNCTX", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Kidney - Medulla" = list(
    "tissue_abbrv" = "KDNMDL", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Liver" = list(
    "tissue_abbrv" = "LIVER", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Lung" = list(
    "tissue_abbrv" = "LUNG", 
    "tissue_color_hex" = "9ACD32", 
    "tissue_color_rgb" = "154,205,50"
  ), 
  "Minor salivary gland" = list(
    "tissue_abbrv" = "SLVRYG", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Muscle - Skeletal" = list(
    "tissue_abbrv" = "MSCLSK", 
    "tissue_color_hex" = "7A67EE", 
    "tissue_color_rgb" = "122,103,238"
  ), 
  "Nerve - Tibial" = list(
    "tissue_abbrv" = "NERVET", 
    "tissue_color_hex" = "FFD700", 
    "tissue_color_rgb" = "255,215,0"
  ), 
  "Ovary" = list(
    "tissue_abbrv" = "OVARY", 
    "tissue_color_hex" = "FFB6C1", 
    "tissue_color_rgb" = "255,182,193"
  ), 
  "Pancreas" = list(
    "tissue_abbrv" = "PNCREAS", 
    "tissue_color_hex" = "CD9B1D", 
    "tissue_color_rgb" = "205,155,29"
  ), 
  "Pituitary" = list(
    "tissue_abbrv" = "PTTARY", 
    "tissue_color_hex" = "B4EEB4", 
    "tissue_color_rgb" = "180,238,180"
  ), 
  "Prostate" = list(
    "tissue_abbrv" = "PRSTTE", 
    "tissue_color_hex" = "D9D9D9", 
    "tissue_color_rgb" = "217,217,217"
  ), 
  "Skin - Not sun exposed (suprapubic)" = list(
    "tissue_abbrv" = "SKINNS", 
    "tissue_color_hex" = "3A5FCD", 
    "tissue_color_rgb" = "58,95,205"
  ), 
  "Skin - Sun exposed (lower leg)" = list(
    "tissue_abbrv" = "SKINS", 
    "tissue_color_hex" = "1E90FF", 
    "tissue_color_rgb" = "30,144,255"
  ), 
  "Small intestine - Terminal ileum" = list(
    "tissue_abbrv" = "SNTTRM", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Spleen" = list(
    "tissue_abbrv" = "SPLEEN", 
    "tissue_color_hex" = "CDB79E", 
    "tissue_color_rgb" = "205,183,158"
  ), 
  "Stomach" = list(
    "tissue_abbrv" = "STMACH", 
    "tissue_color_hex" = "FFD39B", 
    "tissue_color_rgb" = "255,211,155"
  ), 
  "Testis" = list(
    "tissue_abbrv" = "TESTIS", 
    "tissue_color_hex" = "A6A6A6", 
    "tissue_color_rgb" = "166,166,166"
  ), 
  "Thyroid" = list(
    "tissue_abbrv" = "THYROID", 
    "tissue_color_hex" = "008B45", 
    "tissue_color_rgb" = "0,139,69"
  ), 
  "Uterus" = list(
    "tissue_abbrv" = "UTERUS", 
    "tissue_color_hex" = "EED5D2", 
    "tissue_color_rgb" = "238,213,210"
  ), 
  "Vagina" = list(
    "tissue_abbrv" = "VAGINA", 
    "tissue_color_hex" = "EED5D2", 
    "tissue_color_rgb" = "238,213,210"
  ), 
  "Whole blood" = list(
    "tissue_abbrv" = "WHLBLD", 
    "tissue_color_hex" = "FF00FF", 
    "tissue_color_rgb" = "255,0,255"
  )
)

# ---------------------------------------------------------------------------
# Loading / processing
# ---------------------------------------------------------------------------

load_psi <- function(psi_table = PSI_TABLE) {
  # The source notebook also read GTEx.cluster_counts.tsv.gz to build `result`
  # (a per-cluster PSI used for a `cluster_psi >= 0.5` filter). That filter
  # produces `cluster_list_`, which never reaches the heatmap, so the 37 MB
  # counts table is not read here.
  read.table(psi_table, sep = '\t', header = TRUE)
}


get_sig_clusters <- function(comparison, rds_dir = RDS_DIR) {
  # Clusters with a significant unproductive splicing change in a tissue pair,
  # restricted to genes that are also differentially expressed in that pair.
  rds <- readRDS(glue::glue(rds_dir, '{comparison}.rds'))

  sig_genes <- rds$dge %>%
    filter(padj <= 1e-2, abs(log2FoldChange) >= 2) %>%
    pull(gene_id) %>% unique()

  rds$ds %>%
    filter(abs(deltapsi) >= 0.2, p.adjust <= 1e-2, itype == 'UP', ctype == 'PR,UP',
           (gene_id %in% sig_genes)) %>%
    pull(cluster) %>% unique()
}


collect_sig_clusters <- function(rds_dir = RDS_DIR) {
  # Union of significant clusters over every tissue pair (Bladder excluded).
  comparisons <- sub("\\.rds$", "", list.files(rds_dir))
  comparisons <- comparisons[!grepl("Bladder", comparisons)]

  cluster_list <- c()
  pb <- txtProgressBar(min = 0, max = length(comparisons), initial = 0)
  stepi <- 0
  for (comparison in comparisons) {
    cluster_list <- union(cluster_list, get_sig_clusters(comparison, rds_dir))
    setTxtProgressBar(pb, stepi)
    stepi <- stepi + 1
  }
  close(pb)

  cluster_list
}


make_ds_df <- function(cluster_list) {
  ds_df <- cluster_list %>% table() %>% as.data.frame()
  colnames(ds_df) <- c('cluster', 'freq')
  ds_df
}


make_X_introns <- function(psi, ds_df, display_names = TISSUE_DISPLAY_NAMES) {
  # fig2_heatmap: one row per unproductive intron of a selected PR,UP cluster.
  up_psi <- psi %>%
    filter(cluster %in% (ds_df %>% pull(cluster)), itype == 'UP', ctype == 'PR,UP')
  up_psi <- up_psi[rowSums(is.na(up_psi)) <= 0, ]

  tissues <- (up_psi %>% colnames())[TISSUE_COLS]
  X <- up_psi[rowSums(is.na(up_psi)) <= 0, tissues]
  rownames(X) <- up_psi[rowSums(is.na(up_psi)) <= 0, ]$intron

  colnames(X) <- display_names
  X
}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

PLOT_READY_VARS <- c('heatmap_X_introns')


run_all <- function(data_dir = 'figure_data') {
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  psi <- load_psi()
  cluster_list <- collect_sig_clusters()
  ds_df <- make_ds_df(cluster_list)

  heatmap_X_introns <- make_X_introns(psi, ds_df)

  data <- list(heatmap_X_introns = heatmap_X_introns)

  for (name in PLOT_READY_VARS) {
    saveRDS(data[[name]], file.path(data_dir, paste0(name, '.rds')))
  }

  data
}


load_plot_data <- function(data_dir = 'figure_data') {
  data <- list()
  for (name in PLOT_READY_VARS) {
    data[[name]] <- readRDS(file.path(data_dir, paste0(name, '.rds')))
  }
  data
}
