## 40_phase10A_GSE171272_EVmiRNA_prepare.R
## Download/process GSE171272 exosomal miRNA (human DLBCL) + build group template

message("=== Phase 10A: GSE171272 exosomal miRNA (human DLBCL) ===")

required_pkgs <- c("GEOquery", "Biobase", "dplyr", "tibble", "readr", "stringr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run:\ninstall.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      "))\nthen re-run this script."
    )
  }
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})

## paths
project_root <- normalizePath(file.path(".."), winslash = "/", mustWork = TRUE)
message("Project root: ", project_root)

raw_dir  <- file.path(project_root, "data", "raw", "GEO", "GSE171272")
proc_dir <- file.path(project_root, "data", "processed")
meta_dir <- file.path(project_root, "metadata")

dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

eset_rds <- file.path(proc_dir, "GSE171272_ExpressionSet.rds")

## ------------------------------------------------------------------
## 1. Download / load ExpressionSet
## ------------------------------------------------------------------
if (file.exists(eset_rds)) {
  message("--- Loading cached ExpressionSet ---")
  eset <- readRDS(eset_rds)
} else {
  message("--- Downloading GSE171272 via GEOquery ---")
  gse_list <- GEOquery::getGEO(
    "GSE171272",
    GSEMatrix = TRUE,
    AnnotGPL  = TRUE,
    destdir   = raw_dir
  )
  
  if (length(gse_list) > 1) {
    ns <- vapply(gse_list, function(x) ncol(Biobase::exprs(x)), numeric(1))
    idx <- which.max(ns)
    eset <- gse_list[[idx]]
    message("  Multiple ExpressionSets detected; using index ", idx, " (max samples).")
  } else {
    eset <- gse_list[[1]]
  }
  
  saveRDS(eset, eset_rds)
  message("  Saved ExpressionSet to:\n  ", eset_rds)
}

expr  <- Biobase::exprs(eset)
pheno <- Biobase::pData(eset)
fdat  <- Biobase::fData(eset)

message(
  "  Expression matrix: ", nrow(expr), " probes × ",
  ncol(expr), " samples"
)

rng <- range(expr, na.rm = TRUE)
message(
  "  Expression range (as provided): [",
  sprintf("%.2f", rng[1]), ", ",
  sprintf("%.2f", rng[2]), "]"
)

## ------------------------------------------------------------------
## 2. Derive miRNA ID per probe
## ------------------------------------------------------------------
anno_tbl <- fdat %>%
  as.data.frame() %>%
  rownames_to_column("probe_id")

candidate_id_cols <- c("miRNA_ID", "miRNA", "miRNA_id", "ID", "ID_REF", "Name", "ProbeName")
id_col <- candidate_id_cols[candidate_id_cols %in% colnames(anno_tbl)][1]

if (length(id_col) == 0 || is.na(id_col)) {
  stop(
    "Could not find a miRNA identifier column in featureData.\n",
    "Available columns include:\n  ",
    paste(colnames(anno_tbl), collapse = ", ")
  )
}

message("  Using feature column '", id_col, "' as miRNA ID.")

anno_tbl <- anno_tbl %>%
  mutate(miRNA_id = .data[[id_col]]) %>%
  filter(!is.na(miRNA_id), miRNA_id != "") %>%
  mutate(miRNA_id = stringr::str_squish(miRNA_id))

## ------------------------------------------------------------------
## 3. Collapse probe-level expr to miRNA-level (mean)
## ------------------------------------------------------------------
expr_tbl <- expr %>%
  as.data.frame() %>%
  rownames_to_column("probe_id") %>%
  inner_join(anno_tbl %>% select(probe_id, miRNA_id), by = "probe_id")

message(
  "  Probes with mapped miRNA IDs: ",
  nrow(expr_tbl), " (out of ", nrow(expr), ")"
)

expr_miR <- expr_tbl %>%
  select(-probe_id) %>%
  group_by(miRNA_id) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(miRNA_id)

message("  Collapsed to ", nrow(expr_miR), " unique miRNAs.")

out_expr <- file.path(proc_dir, "GSE171272_EVmiRNA_expr_log2.tsv")
readr::write_tsv(expr_miR, out_expr)
message("  Saved miRNA expression matrix to:\n  ", out_expr)

## ------------------------------------------------------------------
## 4. Save sample metadata + group template
## ------------------------------------------------------------------
pheno_out <- pheno %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

out_meta <- file.path(meta_dir, "GSE171272_sample_metadata.csv")
readr::write_csv(pheno_out, out_meta)
message("  Saved raw GEO sample metadata to:\n  ", out_meta)

template_path <- file.path(meta_dir, "GSE171272_sample_group_template.csv")

if (!file.exists(template_path)) {
  template <- pheno_out %>%
    select(sample_id) %>%
    mutate(
      group = NA_character_
      # لطفاً خودت این ستون را در اکسل پر کن: DLBCL یا Control
    )
  readr::write_csv(template, template_path)
  message(
    "  Wrote group template (sample_id + empty 'group') to:\n  ",
    template_path,
    "\n  EDIT THIS FILE manually so that each sample has group = 'DLBCL' یا 'Control'."
  )
} else {
  message("  Group template already exists at:\n  ", template_path)
}

message("=== Phase 10A completed. Fill the template, then run Phase 10B. ===")
