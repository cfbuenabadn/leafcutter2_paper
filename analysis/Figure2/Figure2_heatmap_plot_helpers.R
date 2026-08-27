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


save_heatmap_source_data <- function(ht, X, prefix) {
  # Source-data tables for a drawn fig2c heatmap.
  #
  # `ht` must be the object returned by plot_fig2_heatmap(), i.e. already drawn
  # -- the k-means row and column groups do not exist until draw() runs, so this
  # is the only point at which they can be recovered. Writes two files:
  #
  #   <prefix>.tsv          one row per intron: gene, k-means row group, the
  #                         position it occupies in the drawn figure, and the
  #                         z-scored PSI actually plotted in each tissue
  #   <prefix>.tissues.tsv  one row per tissue: k-means column group and drawn
  #                         position
  #
  # The z-scores are recomputed with the same t/scale/t as plot_fig2_heatmap;
  # scale() is deterministic, so they match the figure exactly.
  ro <- ComplexHeatmap::row_order(ht)
  co <- ComplexHeatmap::column_order(ht)
  if (!is.list(ro)) ro <- list(`1` = ro)
  if (!is.list(co)) co <- list(`1` = co)

  Z <- X %>% as.matrix() %>% t() %>% scale() %>% t()

  row_group <- rep(NA_character_, nrow(X))
  row_pos <- rep(NA_integer_, nrow(X))
  pos <- 1L
  for (g in names(ro)) {
    for (i in ro[[g]]) { row_group[i] <- g; row_pos[i] <- pos; pos <- pos + 1L }
  }

  col_group <- rep(NA_character_, ncol(X))
  col_pos <- rep(NA_integer_, ncol(X))
  pos <- 1L
  for (g in names(co)) {
    for (j in co[[g]]) { col_group[j] <- g; col_pos[j] <- pos; pos <- pos + 1L }
  }

  ann <- attr(X, 'row_annotation')
  out <- data.frame(intron = rownames(X), stringsAsFactors = FALSE)
  if (!is.null(ann)) {
    out$cluster <- ann$cluster
    out$gene_name <- ann$gene_name
    out$gene_id <- ann$gene_id
  }
  out$row_group <- row_group
  out$position_in_figure <- row_pos
  out <- cbind(out, as.data.frame(Z, check.names = FALSE))
  out <- out[order(out$position_in_figure), ]

  tissues <- data.frame(tissue = colnames(X),
                        column_group = col_group,
                        position_in_figure = col_pos,
                        stringsAsFactors = FALSE)
  tissues <- tissues[order(tissues$position_in_figure), ]

  # Base R rather than readr::write_tsv: readr pulls in vroom, whose compiled
  # library needs a newer libstdc++ than some environments provide. A TSV writer
  # should not be able to fail for that reason. na = 'NA' matches write_tsv().
  write_tsv_base <- function(df, path) {
    write.table(df, path, sep = '\t', quote = FALSE, row.names = FALSE,
                col.names = TRUE, na = 'NA')
  }
  write_tsv_base(out, paste0(prefix, '.tsv'))
  write_tsv_base(tissues, paste0(prefix, '.tissues.tsv'))
  cat(sprintf("  %s.tsv          %d introns x %d tissues, %d row groups\n",
              basename(prefix), nrow(out), ncol(Z), length(ro)))
  cat(sprintf("  %s.tissues.tsv  %d tissues, %d column groups\n",
              basename(prefix), nrow(tissues), length(co)))
  invisible(list(rows = out, tissues = tissues))
}


write_up_clusters <- function(ht, X, path, group_labels = NULL, group_patterns = NULL) {
  # UP_clusters.tsv -- the input GO_analysis.ipynb reads.
  #
  # One row per unproductive intron, in the matrix's own row order, carrying the
  # z-scored PSI actually plotted plus the k-means group the intron falls in.
  # `ht` must already be drawn: the k-means groups do not exist before draw().
  #
  # Groups are indexed by their top-to-bottom order in the drawn figure, which
  # is reproducible; the raw k-means labels are arbitrary and are not used.
  #
  # `group_labels`, if given, renames those positions -- element i is the label
  # for the i-th group from the top. Use it to write the manuscript's group
  # numbering (I-VI), which is assigned by tissue pattern and so does not follow
  # figure order. `group_patterns` adds a human-readable `pattern` column
  # alongside, so the file says what each group is without a lookup elsewhere.
  # With neither, `Cluster` is simply the figure position, which is what the
  # archived published file used.
  ro <- ComplexHeatmap::row_order(ht)
  if (!is.list(ro)) ro <- list(`1` = ro)

  # row_order() returns the groups already in top-to-bottom figure order, and
  # each element holds that group's row indices in drawn order.
  row_group <- rep(NA_integer_, nrow(X))
  for (k in seq_along(ro)) row_group[ro[[k]]] <- k

  Z <- X %>% as.matrix() %>% t() %>% scale() %>% t()
  ann <- attr(X, 'row_annotation')

  out <- data.frame(intron = rownames(X),
                    cluster = ann$cluster,
                    itype = 'UP',        # make_X_introns filters to these,
                    ctype = 'PR,UP',     # so they are constant by construction
                    gene_name = ann$gene_name,
                    gene_id = ann$gene_id,
                    stringsAsFactors = FALSE, check.names = FALSE)
  # 17 significant digits so the doubles round-trip exactly
  Zc <- as.data.frame(lapply(as.data.frame(Z, check.names = FALSE),
                             function(v) vapply(v, function(x) sprintf('%.17g', x), '')),
                      check.names = FALSE, stringsAsFactors = FALSE)
  colnames(Zc) <- colnames(Z)
  out <- cbind(out, Zc)
  if (!is.null(group_labels)) {
    stopifnot(length(group_labels) == length(ro))
    out$Cluster <- group_labels[row_group]
  } else {
    out$Cluster <- row_group
  }
  if (!is.null(group_patterns)) {
    stopifnot(length(group_patterns) == length(ro))
    out$pattern <- group_patterns[row_group]
  }

  write.table(out, path, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
  invisible(path)
}
