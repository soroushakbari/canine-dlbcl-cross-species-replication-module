## 38_phase9B_GSE31312_module_scores_Tier2.R
## Compute human DLBCL Tier2 module scores on GSE31312 (GSE56315-based signature)

message("=== Phase 9B: GSE31312 Tier2 module scores (GSE56315 signature) ===")

## ------------------------------------------------------------
## 1. packages
## ------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tibble", "stringr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) و دوباره اسکریپت را اجرا کن."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

## ------------------------------------------------------------
## 2. paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

proc_dir   <- file.path(project_root, "data", "processed")
meta_dir   <- file.path(project_root, "metadata")
sig_dir    <- file.path(project_root, "results", "tables", "signatures")
ms_dir     <- file.path(project_root, "results", "tables", "module_scores")

if (!dir.exists(ms_dir)) {
  dir.create(ms_dir, recursive = TRUE)
  message("Created module_scores directory:\n  ", ms_dir)
}

expr_path   <- file.path(proc_dir, "GSE31312_expr_log2_qcfiltered.tsv")
meta_path   <- file.path(meta_dir, "GSE31312_sample_metadata.csv")

up_sig_path   <- file.path(sig_dir, "human_signature_GSE56315_Tier2_broad_up.txt")
down_sig_path <- file.path(sig_dir, "human_signature_GSE56315_Tier2_broad_down.txt")

if (!file.exists(expr_path)) {
  stop("Expression file for GSE31312 not found at:\n  ", expr_path)
}
if (!file.exists(meta_path)) {
  stop("Sample metadata for GSE31312 not found at:\n  ", meta_path)
}
if (!file.exists(up_sig_path)) {
  stop("Tier2 up-signature file not found at:\n  ", up_sig_path)
}
if (!file.exists(down_sig_path)) {
  stop("Tier2 down-signature file not found at:\n  ", down_sig_path)
}

## ------------------------------------------------------------
## 3. load expression (gene-level) & metadata
## ------------------------------------------------------------
message("\n--- Step 1: Load GSE31312 expression (gene_symbol level) ---")
expr_df <- readr::read_tsv(expr_path, show_col_types = FALSE)

if (!"gene_symbol" %in% names(expr_df)) {
  stop("Expression table must contain a 'gene_symbol' column.")
}

gene_symbols <- expr_df$gene_symbol
expr_mat <- expr_df %>%
  select(-gene_symbol) %>%
  as.matrix()

rownames(expr_mat) <- gene_symbols

message("  Expression matrix: ", nrow(expr_mat), " genes × ",
        ncol(expr_mat), " samples")

## sample metadata
meta_df <- readr::read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

message("  Sample metadata rows: ", nrow(meta_df))

## ------------------------------------------------------------
## 4. load Tier2 signature (human) 
## ------------------------------------------------------------
message("\n--- Step 2: Load Tier2 BROAD human signature (GSE56315) ---")

sig_up <- readr::read_lines(up_sig_path)
sig_down <- readr::read_lines(down_sig_path)

sig_up <- unique(stringr::str_trim(sig_up))
sig_down <- unique(stringr::str_trim(sig_down))

sig_up <- sig_up[sig_up != ""]
sig_down <- sig_down[sig_down != ""]

message("  Tier2 up-signature genes   : ", length(sig_up))
message("  Tier2 down-signature genes : ", length(sig_down))

## ------------------------------------------------------------
## 5. gene overlap & z-score per gene
## ------------------------------------------------------------
message("\n--- Step 3: Compute per-gene z-scores across samples ---")

## ensure numeric matrix
storage.mode(expr_mat) <- "double"

## per-gene mean & sd
gene_means <- rowMeans(expr_mat, na.rm = TRUE)
gene_sds   <- apply(expr_mat, 1L, stats::sd, na.rm = TRUE)

## avoid division by zero
zero_sd <- which(gene_sds == 0 | is.na(gene_sds))
if (length(zero_sd) > 0L) {
  message("  Genes with zero/NA SD across samples (will be set to NA in z-matrix): ",
          length(zero_sd))
  gene_sds[zero_sd] <- NA_real_
}

z_mat <- (expr_mat - gene_means) / gene_sds

## ------------------------------------------------------------
## 6. overlap with Tier2 signature
## ------------------------------------------------------------
message("\n--- Step 4: Overlap with Tier2 signature genes ---")

genes_in_data <- rownames(z_mat)

up_in_data   <- intersect(sig_up,   genes_in_data)
down_in_data <- intersect(sig_down, genes_in_data)

message("  Overlap with expression:")
message("    up   : ", length(up_in_data),   " / ", length(sig_up))
message("    down : ", length(down_in_data), " / ", length(sig_down))

if (length(up_in_data) < 50L || length(down_in_data) < 50L) {
  warning("Overlap between signature and GSE31312 is quite small; ",
          "module scores may be noisy.")
}

## ------------------------------------------------------------
## 7. compute module scores
## ------------------------------------------------------------
message("\n--- Step 5: Compute module scores (z-scored) ---")

if (length(up_in_data) == 0L || length(down_in_data) == 0L) {
  stop("No overlap with signature genes; cannot compute module scores.")
}

up_scores <- colMeans(z_mat[up_in_data, , drop = FALSE], na.rm = TRUE)
down_scores <- colMeans(z_mat[down_in_data, , drop = FALSE], na.rm = TRUE)

score_net <- up_scores - down_scores

scores_df <- tibble::tibble(
  sample_id    = colnames(z_mat),
  score_up_z   = as.numeric(up_scores),
  score_down_z = as.numeric(down_scores),
  score_net_z  = as.numeric(score_net)
)

message("  Computed module scores for ", nrow(scores_df), " samples.")

## ------------------------------------------------------------
## 8. merge with metadata
## ------------------------------------------------------------
message("\n--- Step 6: Merge scores with sample metadata ---")

scores_with_meta <- scores_df %>%
  left_join(meta_df, by = c("sample_id" = "sample_id"))

n_missing_meta <- sum(is.na(scores_with_meta[, setdiff(names(meta_df), "sample_id")]))
if (n_missing_meta > 0L) {
  message("  WARNING: Some metadata fields are NA after merge (this is expected if ",
          "metadata contain many optional columns).")
}

## ------------------------------------------------------------
## 9. save outputs
## ------------------------------------------------------------
out_scores <- file.path(
  ms_dir,
  "module_scores_GSE31312_Tier2_GSE56315_signature_zscore.tsv"
)

out_scores_meta <- file.path(
  ms_dir,
  "module_scores_GSE31312_Tier2_GSE56315_signature_zscore_with_metadata.tsv"
)

readr::write_tsv(scores_df, out_scores)
message("  Saved module scores to:\n  ", out_scores)

readr::write_tsv(scores_with_meta, out_scores_meta)
message("  Saved module scores + metadata to:\n  ", out_scores_meta)

message("\n=== Phase 9B (GSE31312 module scores) completed successfully. ===")
