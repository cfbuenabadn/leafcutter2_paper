# Data loading for the R panels of Figure 5.
#
# Created with Figure code cleaner.
# Source notebook: ../Fig5_R.ipynb
#
# Figures served by this file:
#   * Fig5A-left   QQ plot, u-sQTL vs p-sQTL against eQTL
#   * Fig5A-right  QQ plot, u-sQTL vs p-sQTL against pQTL
#   * Fig5C_wide   PTWAS gene-index plot
#
# Only data loading and processing live here; everything that draws goes in
# Figure5_R_plot_helpers.R.

BASE <- "/project/yangili1/cfbuenabadn/leafcutter2_paper"
RU_PLOTS <- file.path(BASE, "code/Ru_plots")
PLOTS_DIR <- file.path(BASE, "analysis/Figure5/plots")

QQ_EQTL_RDS  <- file.path(RU_PLOTS, "leafcutter2_fig4a_qqplot_eQTL.rds")
QQ_PQTL_RDS  <- file.path(RU_PLOTS, "leafcutter2_fig4a_qqplot_pQTL.rds")
FIG5B_RDS    <- file.path(RU_PLOTS, "leafcutter2_fig4b.ptwas_new.rds")

# Shared p-sQTL / u-sQTL colours (Fig5_R.ipynb)
UP_PALETTE <- c("#1f77b4", "#d62728")

# Fig5C_wide colours: no colocalization / the two colocalization classes
FIG5B_PALETTE <- c("#1f77b4", "#d62728", "#7f7f7f")


load_qq_pvalues <- function(rds_path, pro_slot, unpro_slot) {
  # The p-value vectors behind one QQ panel, named as they appear in the legend.
  rds <- readRDS(rds_path)
  list(`p-sQTL` = rds[[pro_slot]][, "pvalue"],
       `u-sQTL` = rds[[unpro_slot]][, "pvalue"])
}


load_Fig5B <- function(rds_path = FIG5B_RDS) {
  # PTWAS results for Fig5C_wide.
  #
  # Fig5_R.ipynb recodes NA colocalization as the string "no colocalization"
  # in the cell between the Fig5C and Fig5C_wide plots, so Fig5C_wide is drawn
  # from the recoded table. Without this the NA points drop out of the colour
  # scale and the labelled subset changes.
  suppressPackageStartupMessages(library(dplyr))
  Fig5B <- readRDS(rds_path)
  Fig5B %>% mutate(coloc = ifelse(is.na(coloc), "no colocalization", coloc))
}


run_all <- function(data_dir = "figure_data") {
  # Load every plot-ready object, cache it, and return it.
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

  Fig5A_left  <- load_qq_pvalues(QQ_EQTL_RDS, "eqtl_lf2_pro",  "eqtl_lf2_unpro")
  Fig5A_right  <- load_qq_pvalues(QQ_PQTL_RDS, "pQTL_lf2_pro",  "pQTL_lf2_unpro")
  Fig5B    <- load_Fig5B()

  saveRDS(Fig5A_left, file.path(data_dir, "Fig5A_left.rds"))
  saveRDS(Fig5A_right, file.path(data_dir, "Fig5A_right.rds"))
  saveRDS(Fig5B,   file.path(data_dir, "Fig5B.rds"))

  list(Fig5A_left = Fig5A_left, Fig5A_right = Fig5A_right, Fig5B = Fig5B)
}


load_plot_data <- function(data_dir = "figure_data") {
  list(Fig5A_left = readRDS(file.path(data_dir, "Fig5A_left.rds")),
       Fig5A_right = readRDS(file.path(data_dir, "Fig5A_right.rds")),
       Fig5B   = readRDS(file.path(data_dir, "Fig5B.rds")))
}
