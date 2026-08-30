#!/usr/bin/env Rscript
# Source Data tables for the R-drawn panels: Fig. 2c and 2c_confounders, and the
# three Figure 5 panels built from RDS. Companion to make_source_data.py.
#
# Fig. 2c already had source-data TSVs written by save_heatmap_source_data();
# they are copied to source_data/ as CSV so every panel is in one place and one
# format.

base <- dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
if (is.na(base) || base == '') base <- 'analysis'

out_dir <- function(fig) {
  d <- file.path(base, fig, 'source_data'); dir.create(d, showWarnings = FALSE, recursive = TRUE); d
}
report <- function(df, fig, panel) {
  p <- file.path(out_dir(fig), paste0(panel, '.csv'))
  write.csv(df, p, row.names = FALSE, na = 'NA')
  cat(sprintf('   %-44s %6d rows x %2d cols\n', paste0(panel, '.csv'), nrow(df), ncol(df)))
}

cat('Figure2 (R panels)\n')
for (nm in c('fig2c', 'fig2c_confounders')) {
  f <- file.path(base, 'Figure2', 'figure_data', paste0(nm, '_source_data.tsv'))
  if (file.exists(f)) report(read.delim(f, check.names = FALSE), 'Figure2', nm)
  ft <- file.path(base, 'Figure2', 'figure_data', paste0(nm, '_source_data.tissues.tsv'))
  if (file.exists(ft)) report(read.delim(ft, check.names = FALSE), 'Figure2',
                              paste0(nm, '_tissue_groups'))
}

cat('Figure5 (R panels)\n')
# 5a -- the nominal P values behind each QQ curve, one row per test
for (side in c('left', 'right')) {
  f <- file.path(base, 'Figure5', 'figure_data', sprintf('Fig5A_%s.rds', side))
  if (!file.exists(f)) next
  o <- readRDS(f)
  df <- do.call(rbind, lapply(names(o), function(k)
    data.frame(sqtl_class = k, nominal_p = as.numeric(o[[k]]), stringsAsFactors = FALSE)))
  # what the QQ plot actually draws, per class: observed vs expected -log10 P
  df <- do.call(rbind, lapply(split(df, df$sqtl_class), function(d) {
    d <- d[order(d$nominal_p), ]
    d$observed_minus_log10_p <- -log10(d$nominal_p)
    d$expected_minus_log10_p <- -log10(ppoints(nrow(d)))
    d
  }))
  report(df, 'Figure5', sprintf('Fig5A-%s', side))
}

f <- file.path(base, 'Figure5', 'figure_data', 'Fig5B.rds')
if (file.exists(f)) report(as.data.frame(readRDS(f)), 'Figure5', 'Fig5B')
