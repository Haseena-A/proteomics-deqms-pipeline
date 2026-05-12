# =============================================================================
# R/visualize.R
# QC and results visualisation for proteomics DE analysis
#
# All ggplot2 plots use theme_bw() + custom styling (ggprism optional).
# Heatmaps use ComplexHeatmap::Heatmap().
# Both PNG and SVG are saved for every plot.
#
# Main entry point:
#   generate_all_plots(fit_deqms, deqms_results, protein_matrix,
#                      metadata, output_dir, raw_matrix)
# =============================================================================

# ── Package loading ──────────────────────────────────────────────────────────

.load_viz_deps <- function() {
  pkgs_cran <- c("ggplot2", "ggrepel", "circlize")
  pkgs_bioc <- c("ComplexHeatmap")

  for (p in pkgs_cran) .ensure_pkg(p, bioc = FALSE)
  for (p in pkgs_bioc) .ensure_pkg(p, bioc = TRUE)

  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(ComplexHeatmap)
    library(circlize)
  })

  # ggprism is optional — graceful fallback to theme_bw
  if (requireNamespace("ggprism", quietly = TRUE)) {
    library(ggprism)
    .base_theme <<- ggprism::theme_prism(base_size = 12)
  } else {
    .base_theme <<- ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92")
      )
  }

  # matrixStats for colMedians
  .ensure_pkg("matrixStats", bioc = FALSE)
  library(matrixStats)
}

# Shared colour palette for conditions (up to 8)
.COND_COLS <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                "#FF7F00", "#A65628", "#F781BF", "#999999")

# =============================================================================
# File-saving helpers
# =============================================================================

#' Save a ggplot as PNG + SVG
#' @keywords internal
.save_plot <- function(plot, base_path, width = 8, height = 6, dpi = 300) {
  png_path <- sub("\\.(svg|png)$", ".png", base_path)
  ggplot2::ggsave(png_path, plot = plot,
                  width = width, height = height, dpi = dpi)
  message("   Saved: ", png_path)

  svg_path <- sub("\\.(svg|png)$", ".svg", base_path)
  tryCatch(
    ggplot2::ggsave(svg_path, plot = plot, width = width, height = height),
    error = function(e) {
      tryCatch({
        grDevices::svg(svg_path, width = width, height = height)
        print(plot)
        grDevices::dev.off()
        message("   Saved: ", svg_path)
      }, error = function(e2) message("   (SVG export skipped)"))
    }
  )
}

#' Save a ComplexHeatmap as PNG + SVG
#' @keywords internal
.save_heatmap <- function(ht, base_path, width = 10, height = 8, dpi = 300) {
  png_path <- sub("\\.(svg|png)$", ".png", base_path)
  grDevices::png(png_path, width = width, height = height,
                 units = "in", res = dpi)
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
  message("   Saved: ", png_path)

  svg_path <- sub("\\.(svg|png)$", ".svg", base_path)
  tryCatch({
    grDevices::svg(svg_path, width = width, height = height)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
    message("   Saved: ", svg_path)
  }, error = function(e) message("   (SVG export skipped for heatmap)"))
}

# =============================================================================
# Individual plot functions
# =============================================================================

#' Intensity distribution boxplot (before vs after normalization)
#'
#' @param raw_matrix  Pre-normalisation log2 matrix (proteins × samples)
#' @param norm_matrix Post-normalisation log2 matrix
#' @param metadata    Sample metadata with 'condition' column
#' @param output_dir  Output directory
#' @param width,height Plot dimensions in inches
#' @export
plot_intensity_distribution <- function(raw_matrix, norm_matrix, metadata,
                                        output_dir = "results",
                                        width = 11, height = 6) {
  message("\n   Plotting intensity distribution...")

  .reshape_for_boxplot <- function(mat, stage_label) {
    data.frame(
      sample    = rep(colnames(mat), each = nrow(mat)),
      intensity = as.vector(mat),
      stage     = stage_label,
      stringsAsFactors = FALSE
    )
  }

  df <- rbind(
    .reshape_for_boxplot(raw_matrix,  "Before normalisation"),
    .reshape_for_boxplot(norm_matrix, "After normalisation")
  )
  df$stage     <- factor(df$stage,
                         levels = c("Before normalisation", "After normalisation"))
  df$condition <- metadata[df$sample, "condition"]

  n_conds <- length(levels(metadata$condition))
  cond_cols <- setNames(.COND_COLS[seq_len(n_conds)], levels(metadata$condition))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = intensity,
                                         fill = condition)) +
    ggplot2::geom_boxplot(outlier.size = 0.4, alpha = 0.75, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = cond_cols) +
    ggplot2::facet_wrap(~stage, scales = "free_y") +
    .base_theme +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text   = ggplot2::element_text(face = "bold", size = 12)
    ) +
    ggplot2::labs(
      title = "Protein Intensity Distribution",
      x     = "Sample",
      y     = "Log\u2082 Intensity",
      fill  = "Condition"
    )

  .save_plot(p, file.path(output_dir, "intensity_distribution.png"),
             width = width, height = height)
}

# -----------------------------------------------------------------------------
#' Missing value pattern heatmap
#'
#' @param intensity_matrix Log2 matrix (may contain NA)
#' @param metadata Sample metadata
#' @param output_dir Output directory
#' @param width,height Plot dimensions
#' @export
plot_missing_values <- function(intensity_matrix, metadata,
                                output_dir = "results",
                                width = 10, height = 8) {
  message("\n   Plotting missing value heatmap...")

  present_mat <- ifelse(is.na(intensity_matrix), 0L, 1L)
  has_missing <- rowSums(present_mat == 0L) > 0

  if (!any(has_missing)) {
    message("   No missing values — skipping heatmap")
    return(invisible(NULL))
  }

  # Limit to ≤200 most-missing proteins
  n_miss_per_prot <- rowSums(present_mat == 0L)
  top_idx <- head(order(n_miss_per_prot, decreasing = TRUE),
                  min(200L, sum(has_missing)))
  plot_mat <- present_mat[top_idx, , drop = FALSE]

  n_conds  <- length(levels(metadata$condition))
  cond_col <- setNames(.COND_COLS[seq_len(n_conds)], levels(metadata$condition))

  col_anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition = metadata[colnames(plot_mat), "condition"],
    col       = list(Condition = cond_col),
    annotation_name_side = "left"
  )

  col_fun <- circlize::colorRamp2(c(0, 1), c("#2C3E50", "#ECF0F1"))

  ht <- ComplexHeatmap::Heatmap(
    plot_mat,
    name                = "Present",
    col                 = col_fun,
    top_annotation      = col_anno,
    show_row_names      = FALSE,
    column_title        = "Missing Value Pattern",
    column_title_gp     = grid::gpar(fontsize = 14, fontface = "bold"),
    cluster_rows        = TRUE,
    cluster_columns     = TRUE,
    heatmap_legend_param = list(labels = c("Missing", "Present"), at = c(0, 1))
  )

  .save_heatmap(ht, file.path(output_dir, "missing_values_heatmap.png"),
                width = width, height = height)
}

# -----------------------------------------------------------------------------
#' PCA of normalised protein abundances
#'
#' @param protein_matrix Normalised log2 matrix
#' @param metadata Sample metadata
#' @param output_dir Output directory
#' @param label_samples Show sample labels (default TRUE)
#' @param width,height Plot dimensions
#' @export
plot_pca <- function(protein_matrix, metadata,
                     output_dir    = "results",
                     label_samples = TRUE,
                     width = 7, height = 6) {
  message("\n   Plotting PCA...")

  mat  <- protein_matrix[stats::complete.cases(protein_matrix), ]
  pca  <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  vexp <- summary(pca)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    sample    = rownames(pca$x),
    condition = metadata[rownames(pca$x), "condition"]
  )

  n_conds  <- length(levels(metadata$condition))
  cond_col <- setNames(.COND_COLS[seq_len(n_conds)], levels(metadata$condition))

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2,
                                              colour = condition)) +
    ggplot2::geom_point(size = 4.5, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = cond_col) +
    .base_theme +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                       face = "bold", size = 14)) +
    ggplot2::labs(
      title  = "PCA of Protein Abundances",
      x      = sprintf("PC1 (%.1f%% variance)", vexp[1]),
      y      = sprintf("PC2 (%.1f%% variance)", vexp[2]),
      colour = "Condition"
    )

  if (label_samples) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample),
                                      size = 3, max.overlaps = 20,
                                      colour = "grey30")
  }

  .save_plot(p, file.path(output_dir, "pca_plot.png"), width = width, height = height)
}

# -----------------------------------------------------------------------------
#' Sample–sample Pearson correlation heatmap
#'
#' @param protein_matrix Normalised log2 matrix
#' @param metadata Sample metadata
#' @param output_dir Output directory
#' @param width,height Plot dimensions
#' @export
plot_sample_correlation <- function(protein_matrix, metadata,
                                    output_dir = "results",
                                    width = 8, height = 7) {
  message("\n   Plotting sample correlation heatmap...")

  cor_mat <- stats::cor(protein_matrix, use = "pairwise.complete.obs")

  n_conds  <- length(levels(metadata$condition))
  cond_col <- setNames(.COND_COLS[seq_len(n_conds)], levels(metadata$condition))

  col_anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition = metadata[colnames(cor_mat), "condition"],
    col       = list(Condition = cond_col),
    annotation_name_side = "left"
  )

  rng     <- range(cor_mat, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(
    c(rng[1], mean(c(rng[1], 1)), 1),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  ht <- ComplexHeatmap::Heatmap(
    cor_mat,
    name                = "Pearson r",
    col                 = col_fun,
    top_annotation      = col_anno,
    column_title        = "Sample-to-Sample Correlation (Pearson)",
    column_title_gp     = grid::gpar(fontsize = 13, fontface = "bold"),
    row_names_gp        = grid::gpar(fontsize = 8),
    column_names_gp     = grid::gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y,
                      gp = grid::gpar(fontsize = 7))
    }
  )

  .save_heatmap(ht, file.path(output_dir, "sample_correlation_heatmap.png"),
                width = width, height = height)
}

# -----------------------------------------------------------------------------
#' Volcano plot with labeled top hits
#'
#' @param deqms_results DEqMS results data.frame
#' @param output_dir Output directory
#' @param alpha BH-adjusted p-value threshold
#' @param lfc_threshold |log2FC| threshold
#' @param label_top Number of top hits to label
#' @param width,height Plot dimensions
#' @export
plot_volcano <- function(deqms_results,
                         output_dir    = "results",
                         alpha         = 0.05,
                         lfc_threshold = 0.58,
                         label_top     = 12,
                         width = 8, height = 6) {
  message("\n   Plotting volcano plot...")

  df <- deqms_results
  df$neg_log10p <- -log10(df$sca.adj.pval + .Machine$double.eps)
  df$sig <- "Not significant"
  df$sig[df$sca.adj.pval < alpha &  df$logFC >  lfc_threshold] <- "Up"
  df$sig[df$sca.adj.pval < alpha &  df$logFC < -lfc_threshold] <- "Down"
  df$sig <- factor(df$sig, levels = c("Up", "Down", "Not significant"))

  top_hits <- head(df[df$sig != "Not significant", ][
    order(df[df$sig != "Not significant", "sca.adj.pval"]), ], label_top)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = neg_log10p,
                                         colour = sig)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.6, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = c("Up" = "#D62728", "Down" = "#1F77B4",
                 "Not significant" = "#AAAAAA")
    ) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    .base_theme +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                       face = "bold", size = 14)) +
    ggplot2::labs(
      title  = "Volcano Plot (DEqMS)",
      x      = "Log\u2082 Fold Change",
      y      = "-Log\u2081\u2080 Adjusted P-value",
      colour = "Significance"
    )

  if (nrow(top_hits) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data        = top_hits,
      ggplot2::aes(label = protein),
      size        = 2.8,
      max.overlaps = 25,
      colour      = "black",
      fontface    = "italic",
      box.padding = 0.35
    )
  }

  .save_plot(p, file.path(output_dir, "volcano_plot.png"),
             width = width, height = height)
}

# -----------------------------------------------------------------------------
#' MA plot (log2FC vs average expression)
#'
#' @param deqms_results DEqMS results data.frame
#' @param output_dir Output directory
#' @param alpha BH-adjusted p-value threshold
#' @param label_top Number of top hits to label
#' @param width,height Plot dimensions
#' @export
plot_ma <- function(deqms_results,
                    output_dir = "results",
                    alpha      = 0.05,
                    label_top  = 12,
                    width = 8, height = 6) {
  message("\n   Plotting MA plot...")

  df  <- deqms_results
  df$significant <- !is.na(df$sca.adj.pval) & df$sca.adj.pval < alpha

  top_hits <- head(df[df$significant, ][order(df[df$significant, "sca.adj.pval"]), ],
                   label_top)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = AveExpr, y = logFC,
                                         colour = significant)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.6, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#D62728", "FALSE" = "#AAAAAA"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                        colour = "grey30", linewidth = 0.5) +
    .base_theme +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                       face = "bold", size = 14)) +
    ggplot2::labs(
      title  = "MA Plot (DEqMS)",
      x      = "Average Log\u2082 Expression",
      y      = "Log\u2082 Fold Change",
      colour = ""
    )

  if (nrow(top_hits) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data        = top_hits,
      ggplot2::aes(label = protein),
      size        = 2.8, max.overlaps = 25,
      colour      = "black", fontface = "italic"
    )
  }

  .save_plot(p, file.path(output_dir, "ma_plot.png"),
             width = width, height = height)
}

# -----------------------------------------------------------------------------
#' DEqMS diagnostic: protein variance vs PSM count
#'
#' @param fit_deqms DEqMS fit object
#' @param output_dir Output directory
#' @param width,height Plot dimensions
#' @export
plot_variance_psm <- function(fit_deqms,
                              output_dir = "results",
                              width = 8, height = 6) {
  message("\n   Plotting variance vs PSM count (DEqMS diagnostic)...")

  df <- data.frame(
    psm_count    = fit_deqms$count,
    log_variance = log2(fit_deqms$sigma^2)
  )
  df <- df[!is.na(df$psm_count) & !is.na(df$log_variance), ]

  df$psm_bin <- factor(
    pmin(df$psm_count, 20L),
    levels = sort(unique(pmin(df$psm_count, 20L)))
  )
  levels(df$psm_bin)[levels(df$psm_bin) == "20"] <- "20+"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = psm_bin, y = log_variance)) +
    ggplot2::geom_boxplot(fill = "#377EB8", alpha = 0.65,
                          outlier.size = 0.5, na.rm = TRUE) +
    .base_theme +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(size = 8)
    ) +
    ggplot2::labs(
      title = "Protein Variance vs PSM Count (DEqMS Diagnostic)",
      x     = "PSM Count",
      y     = "Log\u2082 Variance"
    )

  .save_plot(p, file.path(output_dir, "variance_psm_plot.png"),
             width = width, height = height)
}

# =============================================================================
# Master function
# =============================================================================

#' Generate all QC and results plots
#'
#' Runs all seven plot functions and saves output to `output_dir`.
#'
#' @param fit_deqms      DEqMS fit object from \code{run_de_analysis()}
#' @param deqms_results  DEqMS results data.frame
#' @param protein_matrix Normalised protein matrix
#' @param metadata       Sample metadata
#' @param output_dir     Output directory (default: "results")
#' @param raw_matrix     Pre-normalisation matrix for intensity distribution plot
#' @param padj_threshold Significance threshold for volcano/MA (default: 0.05)
#' @param lfc_threshold  |log2FC| threshold (default: 0.58)
#'
#' @return Invisible NULL (all plots saved to disk)
#' @export
generate_all_plots <- function(fit_deqms,
                               deqms_results,
                               protein_matrix,
                               metadata,
                               output_dir    = "results",
                               raw_matrix    = NULL,
                               padj_threshold = 0.05,
                               lfc_threshold  = 0.58) {

  .load_viz_deps()

  message("\n=== Generating Proteomics QC and Results Plots ===")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 1. Intensity distribution (needs raw_matrix)
  if (!is.null(raw_matrix)) {
    plot_intensity_distribution(raw_matrix, protein_matrix, metadata,
                                output_dir = output_dir)
  } else {
    message("\n   (Skipping intensity distribution — raw_matrix not provided)")
  }

  # 2. Missing value heatmap
  plot_mat <- if (!is.null(raw_matrix)) raw_matrix else protein_matrix
  plot_missing_values(plot_mat, metadata, output_dir = output_dir)

  # 3. PCA
  plot_pca(protein_matrix, metadata, output_dir = output_dir)

  # 4. Sample correlation
  plot_sample_correlation(protein_matrix, metadata, output_dir = output_dir)

  # 5. Volcano
  plot_volcano(deqms_results, output_dir = output_dir,
               alpha = padj_threshold, lfc_threshold = lfc_threshold)

  # 6. MA
  plot_ma(deqms_results, output_dir = output_dir, alpha = padj_threshold)

  # 7. Variance vs PSM
  plot_variance_psm(fit_deqms, output_dir = output_dir)

  message("\n\u2713 All plots generated — saved to: ", output_dir, "\n")
  invisible(NULL)
}
