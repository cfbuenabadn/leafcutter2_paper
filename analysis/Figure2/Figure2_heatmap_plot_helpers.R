# ---------------------------------------------------------------------------
# Plotting for fig2_heatmap.
#
# Created with Figure code cleaner.
# Source notebooks: ../DS_Heatmap2.ipynb (cell with the commented-out
# fig2_heatmap pdf() calls) and ../DS_Heatmap.ipynb, which carries the same
# cell byte-for-byte.
#
# Takes a plot-ready matrix from Figure2_heatmap_helpers.R and draws it. No
# loading, no processing.
# ---------------------------------------------------------------------------

.libPaths(c('/project/yangili1/cfbuenabadn/R/x86_64-pc-linux-gnu-library/4.1',
            '/software/R-4.1.0-el8-x86_64/lib64/R/library'))

suppressMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(circlize)
  library(grid)
})


# Diverging colour scale for the PSI z-scores (source notebook).
HEATMAP_BREAKS <- seq(-4, 4, 2)
HEATMAP_COLORS <- c("#053061", "#2166AC", "#F7F7F7", "#B2182B", "#67001F")

make_colfunc <- function(breaks = HEATMAP_BREAKS, colors = HEATMAP_COLORS) {
  circlize::colorRamp2(breaks = breaks, colors = colors)
}


plot_fig2_heatmap <- function(X, gtex_colors, seed = 2, use_raster = TRUE,
                              raster_quality = 2, colfunc = NULL) {
  # Rows are z-scored across tissues, then clustered into 6 row groups and
  # 3 column groups by k-means -- hence set.seed(), which fixes both.
  #
  # use_raster = TRUE reproduces the source notebook: the heatmap body is
  # rasterized inside the PDF. Pass use_raster = FALSE for a fully vector body
  # (much larger file, since every cell becomes its own rectangle).
  if (is.null(colfunc)) colfunc <- make_colfunc()

  column_colors <- setNames(
    sapply(colnames(X), function(tissue) paste0("#", gtex_colors[[tissue]]$tissue_color_hex)),
    colnames(X)
  )

  column_ha <- HeatmapAnnotation(
    Tissue = colnames(X),
    col = list(Tissue = column_colors),
    show_legend = FALSE
  )

  set.seed(seed)
  ht <- X %>%
    as.matrix() %>%
    t() %>%
    scale() %>%
    t() %>%
    Heatmap(
      row_title = glue::glue("Unproductive splicing events (N = {nrow(X)})"),
      col = colfunc,
      row_km = 6, row_gap = unit(0.5, "mm"),
      column_km = 3, column_gap = unit(0.5, "mm"),
      show_parent_dend_line = FALSE,
      show_row_names = FALSE, clustering_method_rows = "complete", show_row_dend = FALSE,
      clustering_method_columns = "complete", column_dend_height = unit(1.5, "in"),
      show_column_dend = TRUE,
      heatmap_legend_param = list(title = "PSI (z-score)"),
      use_raster = use_raster, raster_quality = raster_quality,
      bottom_annotation = column_ha,
      width = unit(15, "cm"),
      height = unit(12, "cm"), column_names_gp = gpar(fontsize = 10)
    )

  draw(ht)
}


save_fig2_heatmap <- function(ht, path, width = 12, height = 15, raster_quality = 2) {
  # The source notebook's pdf() / draw() / dev.off() block, as a function.
  pdf(path, width = width, height = height)
  draw(ht, raster_quality = raster_quality)
  dev.off()
  invisible(path)
}
