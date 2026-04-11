## 34_phase8A_build_GSE130874_clinical_PFS_response.R
## Build sample-level PFS table for GSE130874 from supplement

message("=== Phase 8A: Build clinical PFS/response table for GSE130874 ===")

required_pkgs <- c("readxl", "dplyr", "stringr", "readr", "tibble")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) and re-run this script."
    )
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

## ------------------------------------------------------------------
## 1. Paths
## ------------------------------------------------------------------
project_root <- normalizePath("..", winslash = "/", mustWork = TRUE)
meta_dir     <- file.path(project_root, "metadata")

supp_path <- file.path(
  meta_dir,
  "patient outcome and immunophenotype results from 25 dogs with lymphoma.xlsx"
)

if (!file.exists(supp_path)) {
  stop("Supplement file not found at:\n  ", supp_path,
       "\nCopy the .xlsx there and re-run.")
}

message("Project root: ", project_root)
message("Supplement file: ", supp_path)

## ------------------------------------------------------------------
## 2. Read and clean supplement
## ------------------------------------------------------------------
raw <- readxl::read_excel(supp_path, sheet = 1)

# Expect columns: Case, Time to Progression (Days), Censored for PFS
# Drop note rows with NA Case or NA time
clin_clean <- raw %>%
  dplyr::rename(
    Case            = !!names(raw)[1],
    TTP_days        = !!names(raw)[2],
    Censored_for_PFS = !!names(raw)[3]
  ) %>%
  dplyr::filter(!is.na(Case) & !is.na(TTP_days))

message("  Rows with valid Case + time: ", nrow(clin_clean))

## ------------------------------------------------------------------
## 3. Build final clinical table
## ------------------------------------------------------------------
# Case example: "MM-001_S1"
# Expression sample IDs in count matrix: "MM.001_S1"
# => replace '-' with '.' between MM and number

clinical <- clin_clean %>%
  dplyr::mutate(
    sample_id   = stringr::str_replace(Case, "-", "."),
    PFS_days    = as.numeric(TTP_days),
    # According to note: 0 = censored, 1 = progressed
    PFS_status  = as.integer(Censored_for_PFS),
    response_group = NA_character_
  ) %>%
  dplyr::select(sample_id, PFS_days, PFS_status, response_group) %>%
  dplyr::arrange(sample_id)

out_csv <- file.path(meta_dir, "GSE130874_clinical_PFS_response.csv")
readr::write_csv(clinical, out_csv)

message("Saved clinical PFS/response table to:")
message("  ", out_csv)
message("  # dogs (rows): ", nrow(clinical))

## ------------------------------------------------------------------
## 4. Optional sanity check: overlap with module_scores IDs
## ------------------------------------------------------------------
scores_path <- file.path(
  project_root, "results", "tables", "module_scores",
  "module_scores_GSE130874_crossSpeciesBROAD_dog_zscore.tsv"
)

if (file.exists(scores_path)) {
  message("--- Sanity check: overlap with module_scores ---")
  scores <- readr::read_tsv(scores_path, show_col_types = FALSE)
  
  sample_col <- intersect(c("sample_id", "sample", "Sample"), names(scores))
  if (length(sample_col) == 1) {
    ids_scores <- unique(scores[[sample_col]])
    overlap    <- intersect(ids_scores, clinical$sample_id)
    
    message("  Samples in module_scores: ", length(ids_scores))
    message("  Samples in clinical table: ", nrow(clinical))
    message("  Overlap (exact match): ", length(overlap))
    
    if (length(overlap) == 0) {
      warning(
        "No overlap between module_scores sample IDs and clinical sample_id.\n",
        "Check that ID formats match (e.g., 'MM-001_S1' vs 'MM.001_S1')."
      )
    }
  } else {
    message("  Could not auto-detect sample column in module_scores; skipped overlap check.")
  }
} else {
  message("module_scores file for GSE130874 not found; skipping overlap check.")
}

message("=== Phase 8A clinical table build completed ===")
