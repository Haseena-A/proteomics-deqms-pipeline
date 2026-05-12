# Data Directory

This directory holds your input data files. Raw mass spectrometry data files
are **not committed to Git** (see `.gitignore`). Only templates are versioned.

## Input Format: PSM-level (recommended)

File: `psm_data.csv` (or `.tsv`)

| Column | Content |
|--------|---------|
| 1 | Any row identifier (ignored) |
| 2 | `gene` — protein/gene symbol |
| 3+ | Raw TMT/LFQ intensities, one column per sample |

- Values should be **raw (non-log-transformed)** intensities
- Missing values: use empty cells or `NA` (not 0)
- Column names for intensity columns **must match** the `sample` column in your metadata file

See `example/psm_template.csv` for the exact format.

## Input Format: Metadata

File: `metadata.csv`

| Column | Content |
|--------|---------|
| `sample` | Sample name — must match intensity column names in PSM data |
| `condition` | Biological condition/group label |

- The first level of `condition` (alphabetically, or as ordered) is treated as the reference
- To control level ordering, edit `load_data.R` line: `factor(conditions, levels = c(...))`

See `example/metadata_template.csv` for the exact format.

## Input Format: Pre-aggregated Protein Matrix (alternative)

File: `protein_matrix.csv`

- Rows = proteins (row names = protein identifiers)
- Columns = samples (column names match metadata `sample` column)
- Values should be **log2-transformed** intensities (set `LOG2_TRANSFORM = TRUE` if not)

Optionally provide `psm_counts.csv` with columns `protein`, `psm_count` for
optimal DEqMS variance correction.

## Large Files

If your data files exceed GitHub's 100 MB limit, consider:
- [Git LFS](https://git-lfs.github.com/)
- Storing in a cloud bucket and downloading in the pipeline
- The built-in `run_example.R` demo avoids this issue entirely
