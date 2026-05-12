# =============================================================================
# R/load_data.R
# Load, validate, and prepare proteomics input data
#
# Exports:
#   load_example_data()       — TMT10plex A431/miRNA dataset from ExperimentHub
#   load_user_psm_data()      — Load PSM CSV from disk
#   load_user_protein_data()  — Load pre-aggregated protein matrix from disk
#   validate_input_data()     — Validate any intensity matrix + metadata
# =============================================================================

#' Install a package if not already available
#' @keywords internal
.ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      message("Installing Bioconductor package: ", pkg, " ...")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      message("Installing CRAN package: ", pkg, " ...")
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

# -----------------------------------------------------------------------------
#' Load the built-in TMT10plex example dataset
#'
#' Downloads PSM-level data (EH1663) from Bioconductor ExperimentHub.
#' Dataset: A431 human epidermoid carcinoma cells transfected with miRNAs.
#' Publication: PXD004163 (Ohman et al., MCP 2019).
#'
#' @return Named list:
#'   \describe{
#'     \item{psm_data}{PSM-level data.frame (gene column + 10 TMT channels)}
#'     \item{metadata}{Sample metadata (sample, condition) with condition factor}
#'     \item{description}{Dataset description string}
#'   }
#'
#' @examples
#' \dontrun{
#' data <- load_example_data()
#' psm_data <- data$psm_data
#' metadata  <- data$metadata
#' }
#' @export
load_example_data <- function() {
  .ensure_pkg("BiocManager")
  .ensure_pkg("ExperimentHub", bioc = TRUE)
  .ensure_pkg("DEqMS", bioc = TRUE)

  suppressPackageStartupMessages({
    library(ExperimentHub)
  })

  message("\nDownloading TMT10plex proteomics data from ExperimentHub (EH1663)...")
  options(timeout = 300)
  eh      <- ExperimentHub()
  dat.psm <- eh[["EH1663"]]

  # TMT channel → biological condition mapping
  sample_names <- colnames(dat.psm)[3:12]
  conditions   <- c(
    "ctrl",   "miR191", "miR372", "miR519",
    "ctrl",   "miR372", "miR519",
    "ctrl",   "miR191", "miR372"
  )

  metadata <- data.frame(
    sample    = sample_names,
    condition = factor(conditions,
                       levels = c("ctrl", "miR191", "miR372", "miR519")),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  # Sanity checks
  stopifnot(ncol(dat.psm) >= 12)
  stopifnot(all(sample_names %in% colnames(dat.psm)))

  n_psms     <- nrow(dat.psm)
  n_proteins <- length(unique(dat.psm$gene))
  n_samples  <- length(sample_names)

  message("\n\u2713 Example data loaded successfully")
  message("  PSMs     : ", n_psms)
  message("  Proteins : ", n_proteins)
  message("  Samples  : ", n_samples, " (TMT 10-plex)")
  message("  Conditions: ", paste(levels(metadata$condition), collapse = ", "))
  message("  Replicates per condition:")
  print(table(metadata$condition))

  list(
    psm_data    = dat.psm,
    metadata    = metadata,
    description = paste(
      "A431 human epidermoid carcinoma cells treated with miRNAs",
      "(TMT 10-plex, PXD004163, ExperimentHub EH1663)"
    )
  )
}

# -----------------------------------------------------------------------------
#' Load user-provided PSM-level data from a CSV/TSV file
#'
#' Expected format:
#'   - Column 1 : any row ID (ignored)
#'   - Column 2 : `gene` — protein/gene identifier (character)
#'   - Columns 3+ : raw (non-log) intensity values, one per sample
#'
#' @param psm_file  Path to PSM CSV/TSV file
#' @param metadata_file Path to metadata CSV file (columns: sample, condition)
#' @param sep       Field separator (default: auto-detect from extension)
#' @param gene_col  Name of the protein/gene column (default: "gene")
#'
#' @return Named list with `psm_data` and `metadata`
#'
#' @examples
#' \dontrun{
#' data <- load_user_psm_data("data/my_psms.csv", "data/metadata.csv")
#' }
#' @export
load_user_psm_data <- function(psm_file,
                               metadata_file,
                               sep      = NULL,
                               gene_col = "gene") {
  stopifnot(file.exists(psm_file),
            file.exists(metadata_file))

  sep <- .detect_sep(psm_file, sep)

  message("Reading PSM data from: ", psm_file)
  psm_data <- read.csv(psm_file, sep = sep, check.names = FALSE,
                       stringsAsFactors = FALSE)
  message("  Rows (PSMs): ", nrow(psm_data))
  message("  Columns    : ", ncol(psm_data))

  if (!gene_col %in% colnames(psm_data)) {
    stop("PSM file must contain a column named '", gene_col, "'. ",
         "Available columns: ", paste(colnames(psm_data), collapse = ", "))
  }

  # Standardise gene column name
  if (gene_col != "gene") {
    colnames(psm_data)[colnames(psm_data) == gene_col] <- "gene"
  }

  metadata <- .read_metadata(metadata_file)

  message("\n\u2713 PSM data loaded")
  message("  PSMs     : ", nrow(psm_data))
  message("  Proteins : ", length(unique(psm_data$gene)))
  message("  Conditions: ", paste(levels(metadata$condition), collapse = ", "))

  list(psm_data = psm_data, metadata = metadata)
}

# -----------------------------------------------------------------------------
#' Load user-provided pre-aggregated protein matrix
#'
#' The file should be a matrix of log2 intensities (proteins × samples).
#' Row names = protein identifiers, column names = sample names.
#'
#' @param matrix_file   Path to protein matrix CSV (rows = proteins, cols = samples)
#' @param metadata_file Path to metadata CSV (columns: sample, condition)
#' @param psm_file      Optional path to PSM counts CSV (columns: protein, psm_count)
#' @param sep           Field separator (default: auto-detect)
#' @param log2_transform Logical: apply log2 transform? (default FALSE — assumes already log2)
#'
#' @return Named list with `protein_matrix`, `metadata`, `psm_counts` (or NULL)
#' @export
load_user_protein_data <- function(matrix_file,
                                   metadata_file,
                                   psm_file       = NULL,
                                   sep            = NULL,
                                   log2_transform = FALSE) {
  stopifnot(file.exists(matrix_file), file.exists(metadata_file))

  sep <- .detect_sep(matrix_file, sep)

  message("Reading protein matrix from: ", matrix_file)
  mat <- read.csv(matrix_file, sep = sep, check.names = FALSE,
                  stringsAsFactors = FALSE, row.names = 1)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  if (log2_transform) {
    message("  Applying log2 transform...")
    mat[mat <= 0] <- NA
    mat <- log2(mat)
  }

  metadata <- .read_metadata(metadata_file)

  psm_counts <- NULL
  if (!is.null(psm_file) && file.exists(psm_file)) {
    psm_df     <- read.csv(psm_file, stringsAsFactors = FALSE)
    psm_counts <- setNames(as.integer(psm_df$psm_count), psm_df$protein)
    message("  PSM counts loaded for ", length(psm_counts), " proteins")
  } else {
    message("  No PSM counts file provided — DEqMS will use equal weights")
  }

  message("\n\u2713 Protein matrix loaded")
  message("  Proteins : ", nrow(mat))
  message("  Samples  : ", ncol(mat))
  message("  Missing  : ", sprintf("%.1f%%",
           sum(is.na(mat)) / length(mat) * 100))

  list(
    protein_matrix = mat,
    metadata       = metadata,
    psm_counts     = psm_counts
  )
}

# -----------------------------------------------------------------------------
#' Validate intensity matrix and metadata before running the pipeline
#'
#' Checks for: numeric matrix, condition column, sample alignment, all-NA rows,
#' minimum replicates, and PSM count name matching.
#'
#' @param intensity_matrix Numeric matrix (proteins × samples)
#' @param metadata         Sample metadata data.frame
#' @param psm_counts       Optional named numeric vector of PSM counts
#' @param condition_col    Column name for condition in metadata (default: "condition")
#'
#' @return List with validated `intensity_matrix`, `metadata`, `psm_counts`
#' @export
validate_input_data <- function(intensity_matrix,
                                metadata,
                                psm_counts    = NULL,
                                condition_col = "condition") {
  message("Validating input data...")

  # — Matrix checks —
  if (!is.matrix(intensity_matrix) && !is.data.frame(intensity_matrix)) {
    stop("intensity_matrix must be a matrix or data.frame")
  }
  intensity_matrix <- as.matrix(intensity_matrix)
  storage.mode(intensity_matrix) <- "numeric"

  # — Metadata checks —
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame")

  if (!condition_col %in% colnames(metadata)) {
    stop("metadata must contain column '", condition_col, "'. ",
         "Found: ", paste(colnames(metadata), collapse = ", "))
  }

  # — Sample alignment —
  if (!is.null(rownames(metadata)) &&
      all(colnames(intensity_matrix) %in% rownames(metadata))) {
    metadata <- metadata[colnames(intensity_matrix), , drop = FALSE]
  } else if (nrow(metadata) != ncol(intensity_matrix)) {
    stop("Samples in metadata (", nrow(metadata), ") do not match ",
         "columns in intensity_matrix (", ncol(intensity_matrix), ").")
  }

  # — Remove all-NA rows —
  all_na <- rowSums(!is.na(intensity_matrix)) == 0
  if (any(all_na)) {
    message("  Removing ", sum(all_na), " proteins with all-NA values")
    intensity_matrix <- intensity_matrix[!all_na, , drop = FALSE]
  }

  # — Condition factor —
  if (!is.factor(metadata[[condition_col]])) {
    metadata[[condition_col]] <- factor(metadata[[condition_col]])
  }

  # — Minimum replicates warning —
  n_per_cond <- table(metadata[[condition_col]])
  low <- names(n_per_cond[n_per_cond < 2])
  if (length(low) > 0) {
    warning("Conditions with fewer than 2 replicates: ",
            paste(low, collapse = ", "),
            "\nStatistical results may be unreliable.")
  }

  # — PSM count validation —
  if (!is.null(psm_counts)) {
    if (!is.numeric(psm_counts)) stop("psm_counts must be numeric")
    if (!is.null(names(psm_counts))) {
      pct_match <- mean(rownames(intensity_matrix) %in% names(psm_counts))
      if (pct_match < 0.5) {
        warning("Only ", round(pct_match * 100), "% of proteins matched in ",
                "psm_counts. Check protein ID format.")
      }
    }
  }

  message("\u2713 Validation passed")
  message("  Proteins  : ", nrow(intensity_matrix))
  message("  Samples   : ", ncol(intensity_matrix))
  message("  Conditions: ",
          paste(levels(metadata[[condition_col]]), collapse = ", "))
  if (!is.null(psm_counts)) {
    message("  PSM counts: ", length(psm_counts), " proteins")
  }

  list(
    intensity_matrix = intensity_matrix,
    metadata         = metadata,
    psm_counts       = psm_counts
  )
}

# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.detect_sep <- function(path, sep) {
  if (!is.null(sep)) return(sep)
  if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
}

#' @keywords internal
.read_metadata <- function(metadata_file) {
  sep  <- .detect_sep(metadata_file, NULL)
  meta <- read.csv(metadata_file, sep = sep,
                   stringsAsFactors = FALSE, check.names = FALSE)

  if (!"condition" %in% colnames(meta)) {
    stop("Metadata file must contain a 'condition' column.")
  }

  # Use 'sample' column as row names if present
  if ("sample" %in% colnames(meta)) {
    rownames(meta) <- meta$sample
  }

  meta$condition <- factor(meta$condition)
  meta
}
