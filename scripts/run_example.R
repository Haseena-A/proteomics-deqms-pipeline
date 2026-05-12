# =============================================================================
# scripts/run_example.R
# Quick demo — full pipeline on the built-in TMT10plex A431/miRNA dataset
#
# Dataset: A431 cells transfected with miR-191, miR-372, miR-519, or ctrl
#          10 TMT channels, ~2,600 proteins, PXD004163
#
# Runtime: ~3–5 min (first run downloads ~50 MB from ExperimentHub)
#
# USAGE:
#   source("scripts/run_example.R")
# =============================================================================

# ── Load pipeline functions ──────────────────────────────────────────────────
source("R/load_data.R")
source("R/workflow.R")
source("R/visualize.R")
source("R/export.R")

# ── Configuration ────────────────────────────────────────────────────────────
COMPARISON_NAME      <- "miR372-ctrl"   # miR-372 vs control
PADJ_THRESHOLD       <- 0.05
LFC_THRESHOLD        <- 0.58            # log2(1.5) ≈ 1.5-fold
IMPUTATION_METHOD    <- "MinProb"
NORMALIZATION_METHOD <- "median"
OUTPUT_DIR           <- "results"

# ── Run ──────────────────────────────────────────────────────────────────────
cat("\n")
cat(strrep("=", 60), "\n")
cat("  Proteomics DE Pipeline — Example Run\n")
cat("  Dataset : A431 TMT10plex (PXD004163)\n")
cat("  Contrast: ", COMPARISON_NAME, "\n")
cat(strrep("=", 60), "\n\n")

# 1. Load example data
example <- load_example_data()
psm_data <- example$psm_data
metadata <- example$metadata

# 2. Run differential expression
results <- run_de_analysis(
  psm_data             = psm_data,
  metadata             = metadata,
  comparison_name      = COMPARISON_NAME,
  padj_threshold       = PADJ_THRESHOLD,
  lfc_threshold        = LFC_THRESHOLD,
  imputation_method    = IMPUTATION_METHOD,
  normalization_method = NORMALIZATION_METHOD
)

# Unpack for interactive use
protein_matrix <- results$protein_matrix
raw_matrix     <- results$raw_matrix
fit_deqms      <- results$fit_deqms
deqms_results  <- results$deqms_results
psm_counts     <- results$psm_counts

# 3. Generate all plots
generate_all_plots(
  fit_deqms      = fit_deqms,
  deqms_results  = deqms_results,
  protein_matrix = protein_matrix,
  metadata       = metadata,
  output_dir     = OUTPUT_DIR,
  raw_matrix     = raw_matrix,
  padj_threshold = PADJ_THRESHOLD,
  lfc_threshold  = LFC_THRESHOLD
)

# 4. Export results, RDS, and report
export_all(
  fit_deqms       = fit_deqms,
  deqms_results   = deqms_results,
  protein_matrix  = protein_matrix,
  metadata        = metadata,
  psm_counts      = psm_counts,
  comparison_name = COMPARISON_NAME,
  output_dir      = OUTPUT_DIR,
  padj_threshold  = PADJ_THRESHOLD,
  lfc_threshold   = LFC_THRESHOLD
)

cat(strrep("=", 60), "\n")
cat("  Example run complete!\n")
cat("  Results saved to: ", OUTPUT_DIR, "\n\n")
cat("  Quick look at top hits:\n")
cat(strrep("-", 60), "\n")
print(head(deqms_results[, c("protein", "logFC", "sca.adj.pval", "count")], 10))
cat(strrep("=", 60), "\n\n")

# ── Interactive exploration ──────────────────────────────────────────────────
# After running, explore results in R:
#
# # Top significant proteins
# sig <- deqms_results[deqms_results$sca.adj.pval < 0.05 &
#                      abs(deqms_results$logFC) > 0.58, ]
# nrow(sig)
#
# # Reload the full analysis object later
# obj <- readRDS("results/analysis_object.rds")
# obj$deqms_results
