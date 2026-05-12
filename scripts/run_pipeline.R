# =============================================================================
# scripts/run_pipeline.R
# Main entry point — run the full proteomics DE pipeline on your own data
#
# USAGE:
#   1. Edit the CONFIGURATION section below
#   2. source("scripts/run_pipeline.R")
#
# For a quick demo with the built-in dataset, run:
#   source("scripts/run_example.R")
# =============================================================================

# ── 0. Load pipeline functions ───────────────────────────────────────────────
source("R/load_data.R")
source("R/workflow.R")
source("R/visualize.R")
source("R/export.R")

# =============================================================================
# CONFIGURATION — edit these before running
# =============================================================================

# --- Input files (choose ONE option) ---

# Option A: PSM-level data (recommended for DEqMS)
INPUT_MODE      <- "psm"                       # "psm" or "protein"
PSM_FILE        <- "data/my_psms.csv"          # path to your PSM CSV
METADATA_FILE   <- "data/metadata.csv"         # path to your metadata CSV

# Option B: Pre-aggregated protein matrix
# INPUT_MODE      <- "protein"
# MATRIX_FILE     <- "data/protein_matrix.csv"
# METADATA_FILE   <- "data/metadata.csv"
# PSM_COUNTS_FILE <- NULL                      # optional: "data/psm_counts.csv"
# LOG2_TRANSFORM  <- FALSE                     # set TRUE if matrix is on raw scale

# --- Contrast specification ---
# Must be valid limma contrast syntax using exact condition names from metadata
# Examples:
#   "treatment-ctrl"       (simple two-group)
#   "groupB-groupA"
#   "drugHigh-drugLow"
COMPARISON_NAME <- "treatment-ctrl"

# --- Statistical thresholds ---
PADJ_THRESHOLD  <- 0.05    # BH-adjusted p-value cutoff
LFC_THRESHOLD   <- 0.58    # |log2FC| cutoff (0.58 ≈ 1.5-fold change)

# --- Processing options ---
IMPUTATION_METHOD    <- "MinProb"  # "MinProb" (MNAR) or "kNN" (MCAR)
NORMALIZATION_METHOD <- "median"   # "median", "quantile", or "none"

# --- Output ---
OUTPUT_DIR <- "results"

# =============================================================================
# PIPELINE — do not edit below this line
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("  Proteomics DE Pipeline — limma + DEqMS\n")
cat(strrep("=", 60), "\n\n")

# ── Step 1: Load data ────────────────────────────────────────────────────────
if (INPUT_MODE == "psm") {

  cat("Loading PSM-level data...\n")
  data_list  <- load_user_psm_data(PSM_FILE, METADATA_FILE)
  validated  <- validate_input_data(
    intensity_matrix = data_list$psm_data[, -(1:2)],  # check intensity cols
    metadata         = data_list$metadata
  )
  psm_data   <- data_list$psm_data
  metadata   <- data_list$metadata
  protein_matrix <- NULL
  psm_counts <- NULL

} else if (INPUT_MODE == "protein") {

  cat("Loading pre-aggregated protein matrix...\n")
  data_list  <- load_user_protein_data(
    matrix_file     = MATRIX_FILE,
    metadata_file   = METADATA_FILE,
    psm_file        = PSM_COUNTS_FILE,
    log2_transform  = LOG2_TRANSFORM
  )
  validated  <- validate_input_data(
    intensity_matrix = data_list$protein_matrix,
    metadata         = data_list$metadata,
    psm_counts       = data_list$psm_counts
  )
  psm_data       <- NULL
  protein_matrix <- validated$intensity_matrix
  metadata       <- validated$metadata
  psm_counts     <- validated$psm_counts

} else {
  stop("INPUT_MODE must be 'psm' or 'protein'")
}

# ── Step 2: Run DE analysis ──────────────────────────────────────────────────
results <- run_de_analysis(
  psm_data             = psm_data,
  protein_matrix       = protein_matrix,
  metadata             = if (INPUT_MODE == "psm") metadata else validated$metadata,
  psm_counts           = psm_counts,
  comparison_name      = COMPARISON_NAME,
  padj_threshold       = PADJ_THRESHOLD,
  lfc_threshold        = LFC_THRESHOLD,
  imputation_method    = IMPUTATION_METHOD,
  normalization_method = NORMALIZATION_METHOD
)

# Unpack results into the global environment for interactive use
protein_matrix  <- results$protein_matrix
raw_matrix      <- results$raw_matrix
fit_deqms       <- results$fit_deqms
deqms_results   <- results$deqms_results
psm_counts      <- results$psm_counts

# ── Step 3: Generate plots ───────────────────────────────────────────────────
generate_all_plots(
  fit_deqms      = fit_deqms,
  deqms_results  = deqms_results,
  protein_matrix = protein_matrix,
  metadata       = results$metadata %||% metadata,
  output_dir     = OUTPUT_DIR,
  raw_matrix     = raw_matrix,
  padj_threshold = PADJ_THRESHOLD,
  lfc_threshold  = LFC_THRESHOLD
)

# ── Step 4: Export results ───────────────────────────────────────────────────
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
cat("  Pipeline complete! Results in:", OUTPUT_DIR, "\n")
cat(strrep("=", 60), "\n\n")

# ── Null-coalescing operator for internal use ────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b
