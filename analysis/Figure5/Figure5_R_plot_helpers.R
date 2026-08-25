# Plotting functions for the R panels of Figure 5.
#
# Created with Figure code cleaner.
# Source notebook: ../Fig5_R.ipynb
#
# Every function takes plot-ready data (as produced by Figure5_R_helpers.R) and
# returns a ggplot object. No data loading, no computation.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrastr)
  library(foreach)
  library(ggrepel)
})


multiqq <- function(pvalues) {
  # QQ plot of several p-value vectors against a common uniform expectation.
  #
  # NOTE: the expected quantiles come from runif(), so this is stochastic --
  # Fig5_R.ipynb seeds it before each call (set.seed(0) for Fig5A-left,
  # set.seed(3) for Fig5A-right). Keep those seeds to reproduce the panels.
  punif <- -log10(runif(max(sapply(pvalues, length))))
  df <- do.call(rbind, foreach(i = seq_len(length(pvalues))) %do% {
    df <- as.data.frame(qqplot(punif[1:length(pvalues[[i]])],
                               -log10(pvalues[[i]]), plot.it = FALSE))
    df$group <- names(pvalues)[i]
    df
  })
  df$group <- factor(df$group, names(pvalues))
  ggplot(df, aes(x, y, col = group)) +
    geom_point_rast() +
    geom_abline(intercept = 0, slope = 1) +
    theme_bw(base_size = 18) +
    xlab("Expected -log10(p)") +
    ylab("Observed -log10(p)")
}


plot_qq_panel <- function(pvalues, palette) {
  # Fig5A-left / Fig5A-right. Seed before calling (see multiqq).
  multiqq(pvalues) +
    theme(legend.title = element_blank(),
          legend.position = c(0.25, 0.6)) +
    scale_color_manual(values = palette) +
    theme_classic() + geom_point_rast(size = 2) +
    theme_classic(base_size = 8) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = c(0.25, 0.6),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 8),
      axis.ticks.length = unit(0.15, "cm"),
      aspect.ratio = 1
    )
}


plot_Fig5B <- function(Fig5B_results, palette) {
  # Fig5C_wide. Expects the table with NA colocalization already recoded to
  # "no colocalization" (Figure5_R_helpers::load_Fig5B does this).
  threshold_high <- -log10(0.05 / nrow(Fig5B_results))
  threshold_low <- -log10(0.05)

  ggplot(Fig5B_results, aes(x = seq_along(PVAL), y = -log10(PVAL))) +
    geom_point(aes(color = coloc, shape = `Brain regions`, fill = coloc), size = 6) +
    geom_label_repel(
      data = Fig5B_results %>% as.data.frame() %>%
        mutate(gene_index = row_number()) %>%
        filter(-log10(PVAL) > threshold_low &
                 (-log10(PVAL) > threshold_high | coloc != "no colocalization")),
      aes(x = gene_index, y = -log10(PVAL), label = gene_name, color = coloc),
      box.padding = 0.7,
      point.padding = 1,
      size = 7,
      fontface = "italic",   # gene symbols are set in italics
      max.overlaps = Inf, force = 40, force_pull = 20
    ) +
    geom_hline(yintercept = threshold_high, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = threshold_low, color = "grey", linetype = "dashed") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(8, 17, 16, 15, 2, 0, 5)) +
    labs(
      x = "Gene Index",
      y = "-log10(p)",
      color = "Colocalization",
      fill = "Colocalization",
      shape = "Brain regions"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 26),
      axis.title = element_text(size = 26),
      axis.text = element_text(size = 22),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 0.75,
      axis.ticks.length = unit(0.25, "cm")
    )
}


save_panel <- function(plot, name, width, height, plots_dir = PLOTS_DIR, dpi = 300) {
  # Write one panel as svg + pdf + png, matching Fig5_R.ipynb's ggsave calls.
  for (ext in c("svg", "pdf", "png")) {
    ggsave(file.path(plots_dir, paste0(name, ".", ext)), plot,
           width = width, height = height, units = "in", dpi = dpi)
  }
  invisible(file.path(plots_dir, paste0(name, c(".svg", ".pdf", ".png"))))
}
