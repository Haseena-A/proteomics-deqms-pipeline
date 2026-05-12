# =============================================================================
# R/workflow.R
# Core differential expression analysis: PSM aggregation → DEqMS
#
# Inputs  (set in calling env or passed to run_de_analysis()):
#   psm_data            PSM-level data.frame  OR
#   protein_matrix      Pre-aggregated log2 matrix (proteins × samples)
#   metadata            Sample metadata with 'condition' factor column
#   psm_counts          Named integer vector of PSM counts (optional)
#
# Outputs (returned by run_de_analysis()):
#   protein_matrix      Log2-normalised protein matrix
#   raw_matrix          Pre-normalisation log2 matrix (for QC plots)
#   fit_deqms           DEqMS fit object (augmented limma MArrayLM)
#   deqms_results       DEqMS results data.frame (sorted by sca.adj.pval)
#   psm_counts          Named PSM count vector
#   comparison_name     Contrast string
# =============================================================================

#' Run the complete proteomics DE analysis
#'
#' Executes the full pipeline from raw PSM data (or a pre-aggregated matrix)
#' through DEqMS differential expression testing.
#'
#' @param psm_data          PSM-level data.frame (gene col + intensity cols).
#'                          Provide either this OR `protein_matrix`.
#' @param protein_matrix    Pre-aggregated log2 protein intensity matrix.
#' @param metadata          Sample metadata data.frame with 'condition' column.
#' @param psm_counts        Named integer vector of PSM counts per protein.
#'                          Extracted automatically from `psm_data` when available.
#' @param comparison_name   Contrast string in limma syntax, e.g. "treatment-ctrl".
#' @param padj_threshold    BH-adjusted p-value cutoff for summary output.
#' @param lfc_threshold     |log2FC| cutoff for summary (default 0.58 = 1.5-fold).
#' @param imputation_method "MinProb" (MNAR, default) or "kNN" (MCAR).
#' @param normalization_method "median" (default), "quantile", or "none".
#'
#' @return Named list with: protein_matrix, raw_matrix, fit_deqms,
#'         deqms_results, psm_counts, comparison_name
#'
#' @examples
#' \dontrun{
#' data <- load_example_data()
#' results <- run_de_analysis(
#'   psm_data        = data$psm_data,
#'   metadata        = data$metadata,
#'   comparison_name = "miR372-ctrl"
#' )
#' }
#' @export
run_de_analysis <- function(psm_data            = NULL,
                            protein_matrix      = NULL,
                            metadata,
                            psm_counts          = NULL,
                            comparison_name     = "treatment-ctrl",
                            padj_threshold      = 0.05,
                            lfc_threshold       = 0.58,
                            imputation_method   = "MinProb",
                            normalization_method = "median") {

  # ── Dependency check ──────────────────────────────────────────────────────
  .check_workflow_deps()

  suppressPackageStartupMessages({
    library(limma)
    library(DEqMS)
    library(matrixStats)
  })

  message("\n=== Proteomics DE Analysis (limma + DEqMS) ===\n")
  message("  Comparison        : ", comparison_name)
  message("  Imputation        : ", imputation_method)
  message("  Normalization     : ", normalization_method)
  message("  padj threshold    : ", padj_threshold)
  message("  |logFC| threshold : ", lfc_threshold, "\n")

  if (is.null(psm_data) && is.null(protein_matrix)) {
    stop("Provide either 'psm_data' (PSM-level) or 'protein_matrix' (pre-aggregated).")
  }

  # ── Step 1: PSM → Protein aggregation ─────────────────────────────────────
  if (!is.null(psm_data) && is.null(protein_matrix)) {
    agg  <- .aggregate_psms(psm_data, metadata)
    protein_matrix <- agg$protein_matrix
    psm_counts     <- agg$psm_counts
  } else {
    message("Step 1: Using pre-aggregated protein matrix (",
            nrow(protein_matrix), " proteins)")
    if (is.null(psm_counts)) {
      message("  WARNING: No PSM counts — DEqMS will use uniform weights\n")
      psm_counts <- rep(1L, nrow(protein_matrix))
      names(psm_counts) <- rownames(protein_matrix)
    }
  }

  # ── Step 2: Missing value assessment & filtering ───────────────────────────
  message("\nStep 2: Missing value assessment & filtering")
  result_filter <- .filter_missing(protein_matrix, metadata, psm_counts)
  protein_matrix <- result_filter$protein_matrix
  psm_counts     <- result_filter$psm_counts

  # ── Step 3: Save raw matrix for QC ────────────────────────────────────────
  raw_matrix <- protein_matrix

  # ── Step 4: Imputation ────────────────────────────────────────────────────
  message("\nStep 3: Imputation (", imputation_method, ")")
  protein_matrix <- .impute(protein_matrix, method = imputation_method)

  # ── Step 5: Normalization ─────────────────────────────────────────────────
  message("\nStep 4: Normalization (", normalization_method, ")")
  protein_matrix <- .normalize(protein_matrix, method = normalization_method)

  # ── Step 6: limma model ───────────────────────────────────────────────────
  message("\nStep 5: limma linear model")
  fit2 <- .fit_limma(protein_matrix, metadata, comparison_name)

  # ── Step 7: DEqMS variance correction ────────────────────────────────────
  message("\nStep 6: DEqMS PSM-count-aware variance correction")
  fit_deqms <- .apply_deqms(fit2, psm_counts)

  # ── Step 8: Extract results ───────────────────────────────────────────────
  message("\nStep 7: Extracting results")
  deqms_results <- .extract_results(fit_deqms)

  # ── Summary ───────────────────────────────────────────────────────────────
  .print_summary(deqms_results, comparison_name, padj_threshold, lfc_threshold)

  list(
    protein_matrix  = protein_matrix,
    raw_matrix      = raw_matrix,
    fit_deqms       = fit_deqms,
    deqms_results   = deqms_results,
    psm_counts      = psm_counts,
    comparison_name = comparison_name
  )
}

# =============================================================================
# Step implementations
# =============================================================================

#' @keywords internal
.aggregate_psms <- function(psm_data, metadata) {
  message("Step 1: Aggregating PSMs → proteins (medianSweeping)")

  if (!"gene" %in% colnames(psm_data)) {
    stop("PSM data must have a 'gene' column.")
  }

  gene_col    <- which(colnames(psm_data) == "gene")
  sample_cols <- which(colnames(psm_data) %in% rownames(metadata))
  if (length(sample_cols) == 0) {
    # Fallback: assume intensity columns start at column 3
    sample_cols <- seq(3, ncol(psm_data))
    message("  WARNING: no sample names matched metadata rownames; ",
            "assuming intensity columns are 3:", ncol(psm_data))
  }

  # Log2-transform; set log2(0) → NA
  dat.log <- psm_data
  dat.log[, sample_cols] <- log2(psm_data[, sample_cols])
  dat.log[, sample_cols][is.infinite(as.matrix(dat.log[, sample_cols]))] <- NA

  # medianSweeping aggregation
  prot_mat   <- as.matrix(DEqMS::medianSweeping(dat.log, group_col = gene_col))

  # PSM counts per protein
  cnt        <- as.data.frame(table(psm_data$gene))
  psm_counts <- setNames(as.integer(cnt$Freq), as.character(cnt$Var1))

  message("  PSMs aggregated  : ", nrow(psm_data), " → ", nrow(prot_mat), " proteins")
  message("  PSM count range  : ", min(psm_counts), " – ", max(psm_counts))

  list(protein_matrix = prot_mat, psm_counts = psm_counts)
}

#' @keywords internal
.filter_missing <- function(protein_matrix, metadata, psm_counts) {
  n_before <- nrow(protein_matrix)
  n_miss   <- sum(is.na(protein_matrix))
  pct_miss <- n_miss / length(protein_matrix) * 100

  message("  Total proteins  : ", n_before)
  message("  Missing values  : ", sprintf("%.1f%%", pct_miss))

  # Keep proteins with ≥50% non-missing in at least ONE condition
  keep <- rep(FALSE, nrow(protein_matrix))
  for (cond in levels(metadata$condition)) {
    cols <- rownames(metadata)[metadata$condition == cond]
    cols <- intersect(cols, colnames(protein_matrix))
    if (length(cols) == 0) next
    pct_present <- rowSums(!is.na(protein_matrix[, cols, drop = FALSE])) /
                   length(cols)
    keep <- keep | (pct_present >= 0.5)
  }

  protein_matrix <- protein_matrix[keep, , drop = FALSE]
  psm_counts     <- psm_counts[rownames(protein_matrix)]

  message("  After filtering : ", nrow(protein_matrix), " proteins retained (",
          n_before - sum(keep), " removed)")

  list(protein_matrix = protein_matrix, psm_counts = psm_counts)
}

#' @keywords internal
.impute <- function(protein_matrix, method = "MinProb") {
  n_miss <- sum(is.na(protein_matrix))

  if (n_miss == 0) {
    message("  No missing values — skipping imputation")
    return(protein_matrix)
  }

  if (method == "MinProb") {
    # Draw from a down-shifted normal at the 1st percentile
    set.seed(42L)
    for (j in seq_len(ncol(protein_matrix))) {
      na_idx <- is.na(protein_matrix[, j])
      if (!any(na_idx)) next
      obs <- protein_matrix[!na_idx, j]
      protein_matrix[na_idx, j] <- rnorm(
        n    = sum(na_idx),
        mean = quantile(obs, probs = 0.01, na.rm = TRUE),
        sd   = 0.3 * sd(obs, na.rm = TRUE)
      )
    }
    message("  MinProb: imputed ", n_miss, " values",
            " (drawn from low-intensity distribution)")

  } else if (method == "kNN") {
    if (requireNamespace("impute", quietly = TRUE)) {
      res            <- impute::impute.knn(protein_matrix, k = 10)
      protein_matrix <- res$data
      message("  kNN (k=10): imputed ", n_miss, " values")
    } else {
      message("  'impute' package not available — falling back to MinProb")
      return(.impute(protein_matrix, method = "MinProb"))
    }
  } else {
    stop("imputation_method must be 'MinProb' or 'kNN'")
  }

  protein_matrix
}

#' @keywords internal
.normalize <- function(protein_matrix, method = "median") {
  if (method == "median") {
    col_meds       <- matrixStats::colMedians(protein_matrix, na.rm = TRUE)
    protein_matrix <- sweep(protein_matrix, 2, col_meds, "-")
    message("  Median centering applied")

  } else if (method == "quantile") {
    protein_matrix <- limma::normalizeBetweenArrays(protein_matrix,
                                                     method = "quantile")
    message("  Quantile normalization applied")

  } else if (method == "none") {
    message("  No normalization applied")

  } else {
    stop("normalization_method must be 'median', 'quantile', or 'none'")
  }

  protein_matrix
}

#' @keywords internal
.fit_limma <- function(protein_matrix, metadata, comparison_name) {
  design <- model.matrix(~0 + condition, data = metadata)
  colnames(design) <- gsub("^condition", "", colnames(design))

  message("  Design groups   : ", paste(colnames(design), collapse = ", "))
  message("  Contrast        : ", comparison_name)

  fit1      <- limma::lmFit(protein_matrix, design)
  contrasts <- limma::makeContrasts(contrasts = comparison_name, levels = design)
  fit2      <- limma::contrasts.fit(fit1, contrasts = contrasts)
  fit2      <- limma::eBayes(fit2)

  fit2
}

#' @keywords internal
.apply_deqms <- function(fit2, psm_counts) {
  # Align PSM counts with fit rows
  counts <- psm_counts[rownames(fit2$coefficients)]
  counts[is.na(counts)] <- 1L   # default for missing
  fit2$count <- counts

  fit_deqms <- DEqMS::spectraCounteBayes(fit2)

  message("  PSM count range : ",
          min(fit_deqms$count, na.rm = TRUE), " – ",
          max(fit_deqms$count, na.rm = TRUE))

  fit_deqms
}

#' @keywords internal
.extract_results <- function(fit_deqms) {
  res         <- DEqMS::outputResult(fit_deqms, coef_col = 1)
  res$protein <- rownames(res)
  res[order(res$sca.adj.pval), ]
}

#' @keywords internal
.print_summary <- function(deqms_results, comparison_name,
                           padj_threshold, lfc_threshold) {
  n_sig  <- sum(deqms_results$sca.adj.pval < padj_threshold &
                abs(deqms_results$logFC)   > lfc_threshold, na.rm = TRUE)
  n_up   <- sum(deqms_results$sca.adj.pval < padj_threshold &
                deqms_results$logFC        > lfc_threshold, na.rm = TRUE)
  n_down <- sum(deqms_results$sca.adj.pval < padj_threshold &
                deqms_results$logFC        < -lfc_threshold, na.rm = TRUE)
  n_limma <- sum(deqms_results$adj.P.Val < padj_threshold &
                 abs(deqms_results$logFC) > lfc_threshold, na.rm = TRUE)

  message("\n\u2713 DE analysis complete!")
  message("  Proteins tested : ", nrow(deqms_results))
  message("  Comparison      : ", comparison_name)
  message("  Thresholds      : sca.adj.pval < ", padj_threshold,
          ", |logFC| > ", lfc_threshold)
  message("  Significant (DEqMS): ", n_sig,
          "  (", n_up, " up / ", n_down, " down)")
  message("  Significant (limma): ", n_limma)
  message("  Top hit         : ", deqms_results$protein[1],
          "  (logFC = ", sprintf("%.2f", deqms_results$logFC[1]),
          ", sca.adj.pval = ",
          sprintf("%.2e", deqms_results$sca.adj.pval[1]), ")\n")
}

#' @keywords internal
.check_workflow_deps <- function() {
  required <- list(
    list(pkg = "limma",       bioc = TRUE),
    list(pkg = "DEqMS",       bioc = TRUE),
    list(pkg = "matrixStats", bioc = FALSE)
  )
  for (r in required) {
    .ensure_pkg(r$pkg, bioc = r$bioc)
  }
}
