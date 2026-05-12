# =============================================================================
# R/export.R
# Export proteomics DE results: CSVs, RDS, Markdown report, PDF report
#
# Main entry point:
#   export_all(fit_deqms, deqms_results, protein_matrix, metadata, ...)
# =============================================================================

#' Export all proteomics DE results and generate reports
#'
#' Saves: all_results.csv, significant_results.csv, top100_proteins.csv,
#' normalized_protein_matrix.csv, psm_counts.csv, analysis_object.rds,
#' analysis_report.md, and (optionally) analysis_report.pdf.
#'
#' @param fit_deqms       DEqMS fit object from \code{run_de_analysis()}
#' @param deqms_results   DEqMS results data.frame
#' @param protein_matrix  Normalised protein intensity matrix
#' @param metadata        Sample metadata data.frame
#' @param psm_counts      Named integer vector of PSM counts (optional)
#' @param comparison_name Contrast string, e.g. "miR372-ctrl"
#' @param output_dir      Output directory (default: "results")
#' @param padj_threshold  Significance threshold (default: 0.05)
#' @param lfc_threshold   |log2FC| threshold (default: 0.58)
#'
#' @return Invisible list of output file paths
#' @export
export_all <- function(fit_deqms,
                       deqms_results,
                       protein_matrix,
                       metadata,
                       psm_counts      = NULL,
                       comparison_name = "comparison",
                       output_dir      = "results",
                       padj_threshold  = 0.05,
                       lfc_threshold   = 0.58) {

  message("\n=== Exporting Proteomics DE Results ===\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created directory: ", output_dir, "\n")
  }

  # Extract PSM counts from fit object if not provided
  if (is.null(psm_counts) && !is.null(fit_deqms$count)) {
    psm_counts <- fit_deqms$count
  }

  out_files <- list()

  # ── 1. All DEqMS results ────────────────────────────────────────────────
  message("1. All DEqMS results...")
  f <- file.path(output_dir, "all_results.csv")
  write.csv(deqms_results, f, row.names = FALSE)
  out_files$all_results <- f
  message("   Saved: all_results.csv  (", nrow(deqms_results), " proteins)\n")

  # ── 2. Significant results ──────────────────────────────────────────────
  message("2. Significant results...")
  sig <- deqms_results[
    !is.na(deqms_results$sca.adj.pval) &
    deqms_results$sca.adj.pval < padj_threshold &
    abs(deqms_results$logFC)   > lfc_threshold, ]
  sig <- sig[order(sig$sca.adj.pval), ]
  f <- file.path(output_dir, "significant_results.csv")
  write.csv(sig, f, row.names = FALSE)
  out_files$significant_results <- f
  n_up   <- sum(sig$logFC > 0)
  n_down <- sum(sig$logFC < 0)
  message("   Saved: significant_results.csv  (", nrow(sig),
          " proteins: ", n_up, " up / ", n_down, " down)")
  message("   Thresholds: sca.adj.pval < ", padj_threshold,
          ", |logFC| > ", lfc_threshold, "\n")

  # ── 3. Top 100 ──────────────────────────────────────────────────────────
  message("3. Top 100 proteins...")
  top100 <- head(deqms_results[order(deqms_results$sca.adj.pval), ], 100)
  f <- file.path(output_dir, "top100_proteins.csv")
  write.csv(top100, f, row.names = FALSE)
  out_files$top100 <- f
  message("   Saved: top100_proteins.csv\n")

  # ── 4. Normalised matrix ────────────────────────────────────────────────
  message("4. Normalised protein matrix...")
  f <- file.path(output_dir, "normalized_protein_matrix.csv")
  write.csv(protein_matrix, f, row.names = TRUE)
  out_files$normalized_matrix <- f
  message("   Saved: normalized_protein_matrix.csv  (",
          nrow(protein_matrix), " \u00d7 ", ncol(protein_matrix), ")\n")

  # ── 5. PSM counts ───────────────────────────────────────────────────────
  if (!is.null(psm_counts)) {
    message("5. PSM counts...")
    psm_df <- data.frame(
      protein   = names(psm_counts),
      psm_count = as.integer(psm_counts)
    )
    psm_df <- psm_df[order(psm_df$psm_count, decreasing = TRUE), ]
    f <- file.path(output_dir, "psm_counts.csv")
    write.csv(psm_df, f, row.names = FALSE)
    out_files$psm_counts <- f
    message("   Saved: psm_counts.csv  (", nrow(psm_df), " proteins)\n")
  }

  # ── 6. Analysis object (RDS) ────────────────────────────────────────────
  message("6. Saving analysis object (RDS)...")
  analysis_obj <- list(
    fit_deqms       = fit_deqms,
    deqms_results   = deqms_results,
    protein_matrix  = protein_matrix,
    metadata        = metadata,
    psm_counts      = psm_counts,
    comparison_name = comparison_name,
    thresholds      = list(padj = padj_threshold, lfc = lfc_threshold),
    created         = Sys.time()
  )
  f <- file.path(output_dir, "analysis_object.rds")
  saveRDS(analysis_obj, f)
  out_files$rds <- f
  message("   Saved: analysis_object.rds")
  message("   Load with: obj <- readRDS('", f, "')\n")

  # ── 7. Markdown report ──────────────────────────────────────────────────
  message("7. Generating Markdown report...")
  .write_markdown_report(
    deqms_results, metadata, comparison_name, output_dir,
    padj_threshold, lfc_threshold
  )
  out_files$md_report <- file.path(output_dir, "analysis_report.md")
  message("   Saved: analysis_report.md\n")

  # ── 8. PDF report (optional) ────────────────────────────────────────────
  message("8. PDF report (optional — requires rmarkdown + LaTeX)...")
  pdf_path <- tryCatch(
    .generate_pdf_report(
      deqms_results, metadata, comparison_name,
      output_dir, padj_threshold, lfc_threshold
    ),
    error = function(e) {
      message("   PDF skipped: ", conditionMessage(e))
      message("   (Markdown report still available)\n")
      NULL
    }
  )
  if (!is.null(pdf_path)) {
    out_files$pdf_report <- pdf_path
    message("   Saved: analysis_report.pdf\n")
  }

  # ── Summary ──────────────────────────────────────────────────────────────
  message("\n=== Export Complete ===")
  message("Output directory: ", output_dir)
  all_out <- list.files(output_dir, pattern = "\\.(csv|rds|md|pdf)$")
  for (f in all_out) message("  \u2022 ", f)
  message("")

  invisible(out_files)
}

# =============================================================================
# Markdown report
# =============================================================================

#' @keywords internal
.write_markdown_report <- function(deqms_results, metadata,
                                   comparison_name, output_dir,
                                   padj_threshold, lfc_threshold) {
  sig <- deqms_results[
    !is.na(deqms_results$sca.adj.pval) &
    deqms_results$sca.adj.pval < padj_threshold &
    abs(deqms_results$logFC)   > lfc_threshold, ]

  n_total <- nrow(deqms_results)
  n_sig   <- nrow(sig)
  n_up    <- sum(sig$logFC > 0)
  n_down  <- sum(sig$logFC < 0)

  top10 <- head(deqms_results[order(deqms_results$sca.adj.pval), ], 10)

  cond_tab <- as.data.frame(table(metadata$condition))
  cond_lines <- paste0("| ", cond_tab$Var1, " | ", cond_tab$Freq, " |",
                       collapse = "\n")

  # Build top-10 table rows
  top_rows <- apply(top10, 1, function(r) {
    sprintf("| %s | %.3f | %.2e | %.2e | %s |",
            r["protein"],
            as.numeric(r["logFC"]),
            as.numeric(r["sca.adj.pval"]),
            as.numeric(r["adj.P.Val"]),
            r["count"])
  })

  lines <- c(
    "# Proteomics Differential Expression Report",
    "",
    paste("**Date:**", format(Sys.Date(), "%B %d, %Y")),
    paste("**Comparison:**", comparison_name),
    "",
    "---",
    "",
    "## Summary",
    "",
    paste0("| Metric | Value |"),
    paste0("|--------|-------|"),
    paste0("| Proteins tested | ", format(n_total, big.mark = ","), " |"),
    paste0("| Significant (DEqMS) | **", n_sig, "** |"),
    paste0("| &nbsp;&nbsp;Upregulated | ", n_up, " |"),
    paste0("| &nbsp;&nbsp;Downregulated | ", n_down, " |"),
    paste0("| padj threshold | ", padj_threshold, " |"),
    paste0("| |log2FC| threshold | ", lfc_threshold, " |"),
    paste0("| Total samples | ", nrow(metadata), " |"),
    "",
    "### Samples per condition",
    "",
    "| Condition | n |",
    "|-----------|---|",
    cond_lines,
    "",
    "---",
    "",
    "## Methods",
    "",
    "Analysis performed using the **limma + DEqMS** pipeline:",
    "",
    "1. **PSM → protein aggregation** via `medianSweeping()` (DEqMS package)",
    "2. **Missing value filtering** — proteins with >50% missing in all conditions removed",
    paste0("3. **Imputation** — MinProb (draws from bottom 1% of observed distribution; ",
           "appropriate for MNAR data in MS proteomics)"),
    "4. **Normalization** — median centering (per-sample median subtracted)",
    "5. **limma linear model** — `lmFit()` + `contrasts.fit()` + `eBayes()`",
    "6. **DEqMS variance correction** — `spectraCounteBayes()`: PSM-count-aware empirical Bayes",
    "",
    "---",
    "",
    "## Top 10 Differentially Expressed Proteins",
    "",
    "| Protein | log2FC | DEqMS adj.p | limma adj.p | PSM count |",
    "|---------|--------|-------------|-------------|-----------|",
    top_rows,
    "",
    "---",
    "",
    "## Output Files",
    "",
    "| File | Description |",
    "|------|-------------|",
    "| `all_results.csv` | Full DEqMS results for all proteins |",
    "| `significant_results.csv` | Significant hits only |",
    "| `top100_proteins.csv` | Top 100 by DEqMS adjusted p-value |",
    "| `normalized_protein_matrix.csv` | Normalised log2 intensities |",
    "| `psm_counts.csv` | PSM counts per protein |",
    "| `analysis_object.rds` | Full R analysis object for downstream use |",
    "| `volcano_plot.png/svg` | Volcano plot |",
    "| `ma_plot.png/svg` | MA plot |",
    "| `pca_plot.png/svg` | PCA of samples |",
    "| `sample_correlation_heatmap.png/svg` | Sample correlation |",
    "| `missing_values_heatmap.png/svg` | Missing value pattern |",
    "| `variance_psm_plot.png/svg` | DEqMS diagnostic |",
    "",
    "---",
    "",
    "## References",
    "",
    "1. Zhu Y, et al. DEqMS: A Method for Accurate Variance Estimation in Differential",
    "   Protein Expression Analysis. *Mol Cell Proteomics.* 2020;19(6):1047-1057.",
    "   <https://doi.org/10.1074/mcp.TIR119.001646>",
    "",
    "2. Ritchie ME, et al. limma powers differential expression analyses for",
    "   RNA-sequencing and microarray studies. *Nucleic Acids Res.* 2015;43(7):e47.",
    "   <https://doi.org/10.1093/nar/gkv007>",
    ""
  )

  writeLines(lines, file.path(output_dir, "analysis_report.md"))
}

# =============================================================================
# PDF report (optional)
# =============================================================================

#' @keywords internal
.generate_pdf_report <- function(deqms_results, metadata,
                                 comparison_name, output_dir,
                                 padj_threshold, lfc_threshold) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown not installed. Run: install.packages('rmarkdown')")
  }

  # Check for LaTeX
  has_latex <- FALSE
  if (requireNamespace("tinytex", quietly = TRUE)) {
    has_latex <- tinytex::is_tinytex() || nchar(Sys.which("xelatex")) > 0
  }
  if (!has_latex) {
    has_latex <- nchar(Sys.which("pdflatex")) > 0 ||
                 nchar(Sys.which("xelatex"))  > 0
  }
  if (!has_latex) {
    stop("No LaTeX installation found. Install with: tinytex::install_tinytex()")
  }

  n_total <- nrow(deqms_results)
  n_sig   <- sum(deqms_results$sca.adj.pval < padj_threshold &
                 abs(deqms_results$logFC)   > lfc_threshold, na.rm = TRUE)
  n_up    <- sum(deqms_results$sca.adj.pval < padj_threshold &
                 deqms_results$logFC        > lfc_threshold, na.rm = TRUE)
  n_down  <- sum(deqms_results$sca.adj.pval < padj_threshold &
                 deqms_results$logFC        < -lfc_threshold, na.rm = TRUE)

  top20 <- head(deqms_results[order(deqms_results$sca.adj.pval), ], 20)

  # Build figure section
  fig_map <- c(
    "volcano_plot.png"              = "Volcano plot — significant proteins highlighted.",
    "ma_plot.png"                   = "MA plot — log2FC vs average expression.",
    "pca_plot.png"                  = "PCA of normalised protein abundances.",
    "intensity_distribution.png"    = "Intensity distributions before/after normalization.",
    "sample_correlation_heatmap.png" = "Sample-to-sample Pearson correlation.",
    "missing_values_heatmap.png"    = "Missing value pattern across samples.",
    "variance_psm_plot.png"         = "Protein variance vs PSM count (DEqMS diagnostic)."
  )

  fig_chunks <- ""
  for (nm in names(fig_map)) {
    fpath <- normalizePath(file.path(output_dir, nm), mustWork = FALSE)
    if (file.exists(fpath)) {
      fig_chunks <- paste0(
        fig_chunks,
        '\n```{r, out.width="100%", fig.cap="', fig_map[[nm]], '"}\n',
        'knitr::include_graphics("', fpath, '")\n```\n'
      )
    }
  }

  # Serialise top20 table inline (avoid file dependency in Rmd)
  top20_json <- jsonlite::toJSON(top20[, c("protein", "logFC", "sca.adj.pval",
                                            "adj.P.Val", "count")],
                                 digits = 6)

  rmd <- paste0(
'---
title: "Proteomics Differential Expression Report"
subtitle: "limma + DEqMS Pipeline"
date: "', format(Sys.Date(), "%B %d, %Y"), '"
output:
  pdf_document:
    toc: true
    toc_depth: 2
    number_sections: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
```

# Summary

| Metric | Value |
|--------|-------|
| Comparison | ', comparison_name, ' |
| Proteins tested | ', format(n_total, big.mark = ","), ' |
| Significant (DEqMS) | **', n_sig, '** |
| &nbsp;&nbsp;Upregulated | ', n_up, ' |
| &nbsp;&nbsp;Downregulated | ', n_down, ' |
| padj threshold | ', padj_threshold, ' |
| log2FC threshold | ', sprintf("%.2f", lfc_threshold), ' |
| Samples | ', nrow(metadata), ' |

# Methods

1. **PSM → protein**: `medianSweeping()` (PSM log2-transform, median sweep, protein median)
2. **Filtering**: proteins with >50% missing in all conditions removed
3. **Imputation**: MinProb (drawn from low-intensity distribution; MNAR assumption)
4. **Normalization**: median centering
5. **limma**: `lmFit` + `contrasts.fit` + `eBayes`
6. **DEqMS**: `spectraCounteBayes` — PSM-count-aware empirical Bayes prior variance

# Results

## Top 20 Differentially Expressed Proteins

```{r top-proteins}
tbl <- jsonlite::fromJSON(\'', top20_json, '\')
knitr::kable(tbl,
  col.names = c("Protein", "log2FC", "DEqMS adj.p", "limma adj.p", "PSMs"),
  digits     = c(0, 3, 3, 3, 0),
  caption    = "Top 20 proteins by DEqMS adjusted p-value",
  format     = "latex", booktabs = TRUE)
```

## Figures
', fig_chunks, '

# Conclusions

- **', n_sig, '** proteins were significantly differentially expressed
  (DEqMS adj.p < ', padj_threshold, ', |log2FC| > ', sprintf("%.2f", lfc_threshold), ')
  in the comparison **', comparison_name, '**.
- Of these, **', n_up, '** were upregulated and **', n_down, '** downregulated.

# References

1. Zhu Y et al. DEqMS. *Mol Cell Proteomics.* 2020;19(6):1047–1057.
2. Ritchie ME et al. limma. *Nucleic Acids Res.* 2015;43(7):e47.
'
  )

  rmd_path <- file.path(output_dir, "analysis_report.Rmd")
  writeLines(rmd, rmd_path)

  pdf_path <- file.path(output_dir, "analysis_report.pdf")

  rmarkdown::render(
    rmd_path,
    output_file = basename(pdf_path),
    output_dir  = output_dir,
    quiet       = TRUE
  )

  # Tidy up
  for (f in c(rmd_path,
              list.files(output_dir, pattern = "\\.(tex|log|aux)$",
                         full.names = TRUE))) {
    if (file.exists(f)) file.remove(f)
  }

  pdf_path
}
