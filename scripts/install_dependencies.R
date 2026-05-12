# =============================================================================
# scripts/install_dependencies.R
# Install all required and optional packages for the pipeline
#
# USAGE:
#   source("scripts/install_dependencies.R")
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("\n=== Installing proteomics pipeline dependencies ===\n\n")

# ── BiocManager ──────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# ── Required CRAN packages ───────────────────────────────────────────────────
cran_required <- c("matrixStats", "ggplot2", "ggrepel", "circlize")
cat("Installing required CRAN packages:", paste(cran_required, collapse = ", "), "\n")
install.packages(cran_required[!sapply(cran_required, requireNamespace, quietly = TRUE)])

# ── Required Bioconductor packages ───────────────────────────────────────────
bioc_required <- c("limma", "DEqMS", "ComplexHeatmap", "ExperimentHub")
cat("Installing required Bioconductor packages:",
    paste(bioc_required, collapse = ", "), "\n")
BiocManager::install(
  bioc_required[!sapply(bioc_required, requireNamespace, quietly = TRUE)],
  update = FALSE, ask = FALSE
)

# ── Optional packages ────────────────────────────────────────────────────────
cat("\nInstalling optional packages (graceful if unavailable)...\n")

# ggprism: prettier plots (falls back to theme_bw if missing)
tryCatch(
  install.packages("ggprism"),
  error = function(e) message("  ggprism: skipped (", conditionMessage(e), ")")
)

# svglite: higher-quality SVG export
tryCatch(
  install.packages("svglite"),
  error = function(e) message("  svglite: skipped")
)

# impute: kNN imputation alternative
tryCatch(
  BiocManager::install("impute", update = FALSE, ask = FALSE),
  error = function(e) message("  impute: skipped")
)

# rmarkdown + tinytex: PDF report generation
tryCatch({
  install.packages(c("rmarkdown", "knitr", "jsonlite"))
  cat("  rmarkdown installed.\n")
  cat("  For PDF reports, install LaTeX:\n")
  cat("    tinytex::install_tinytex()\n")
}, error = function(e) message("  rmarkdown: skipped"))

# ── Verification ─────────────────────────────────────────────────────────────
cat("\n=== Checking installation ===\n")
all_pkgs <- c(cran_required, bioc_required)
status <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
for (i in seq_along(status)) {
  mark <- if (status[i]) "\u2713" else "\u2717"
  cat(sprintf("  %s  %s\n", mark, names(status)[i]))
}

if (all(status)) {
  cat("\n\u2713 All required packages installed successfully!\n")
  cat("  Run: source('scripts/run_example.R')\n\n")
} else {
  failed <- names(status[!status])
  cat("\n\u2717 Some packages failed to install:", paste(failed, collapse = ", "), "\n")
  cat("  Try installing manually:\n")
  cat("    BiocManager::install(c(", paste0('"', failed, '"', collapse = ", "), "))\n\n")
}
