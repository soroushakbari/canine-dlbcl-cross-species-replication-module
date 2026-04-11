## 37_phase9A_microarray_GSE31312.R
## Process human DLBCL microarray cohort GSE31312 (expression + metadata)

message("=== Phase 9A: Process GSE31312 microarray (human DLBCL) ===")

## ------------------------------------------------------------
## 0. global options (برای دانلود GEO)
## ------------------------------------------------------------
options(timeout = 600)  # 10 دقیقه
options('download.file.method.GEOquery' = 'auto')
options('GEOquery.inmemory.gpl' = FALSE)

## ------------------------------------------------------------
## 1. packages
## ------------------------------------------------------------
required_pkgs <- c(
  "GEOquery", "Biobase", "SummarizedExperiment",
  "readr", "dplyr", "tibble", "stringr"
)

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
  library(GEOquery)
  library(Biobase)
  library(SummarizedExperiment)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

## ------------------------------------------------------------
## 2. paths (اینجا رو فیکس می‌کنیم)
## ------------------------------------------------------------
project_root <- "D:/Research/My Articles/DLBCL drug"
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
message("Project root: ", project_root)

gse_id   <- "GSE31312"
raw_dir  <- file.path(project_root, "data", "raw", "GEO", gse_id)
proc_dir <- file.path(project_root, "data", "processed")
meta_dir <- file.path(project_root, "metadata")

dirs_to_make <- c(raw_dir, proc_dir, meta_dir)
for (dd in dirs_to_make) {
  if (!dir.exists(dd)) {
    dir.create(dd, recursive = TRUE)
    message("Created directory:\n  ", dd)
  }
}

expr_out_path    <- file.path(proc_dir, "GSE31312_expr_log2_qcfiltered.tsv")
se_out_path      <- file.path(proc_dir, "GSE31312_SummarizedExperiment.rds")
sample_out_path  <- file.path(meta_dir, "GSE31312_sample_metadata.csv")
feature_out_path <- file.path(meta_dir, "GSE31312_feature_annotation.csv")

## ------------------------------------------------------------
## 3. download / load GSE31312
## ------------------------------------------------------------
message("\n--- Step 1: Download/load GSE31312 ExpressionSet ---")

gse_list <- GEOquery::getGEO(
  GSEMatrix = TRUE,
  GEO      = gse_id,
  AnnotGPL = TRUE,
  destdir  = raw_dir
)

if (is.list(gse_list) && length(gse_list) > 1L) {
  message("  GSE31312 returned as a list of ", length(gse_list),
          " ExpressionSets; using [[1]].")
  gse <- gse_list[[1]]
} else if (is.list(gse_list)) {
  gse <- gse_list[[1]]
} else {
  gse <- gse_list
}

expr_mat <- Biobase::exprs(gse)
pheno    <- Biobase::pData(gse)
fdat     <- Biobase::fData(gse)

message("  Raw expression matrix: ", nrow(expr_mat), " probes × ",
        ncol(expr_mat), " samples")

## ------------------------------------------------------------
## 4. basic expression sanity check
## ------------------------------------------------------------
expr_range <- range(expr_mat, na.rm = TRUE)
message("  Expression range (before any transform): [",
        sprintf('%.2f', expr_range[1]), ", ",
        sprintf('%.2f', expr_range[2]), "]")

## ------------------------------------------------------------
## 5. derive gene_symbol
## ------------------------------------------------------------
message("\n--- Step 2: Derive gene_symbol from feature data ---")

fdat_df <- fdat %>%
  as.data.frame() %>%
  rownames_to_column("probe_id")

## robust detection (case-insensitive) of a gene-symbol column
symbol_lower_targets <- c(
  "gene symbol", "gene_symbol", "symbol", "genesymbol", "gene.symbol"
)

colnames_lower <- tolower(colnames(fdat_df))
sym_idx <- which(colnames_lower %in% symbol_lower_targets)

if (length(sym_idx) == 0L) {
  stop(
    "No gene symbol-like column found in feature data (case-insensitive search).\n",
    "Available columns (first 20):\n  ",
    paste(head(colnames(fdat_df), 20), collapse = ", ")
  )
}

sym_col <- colnames(fdat_df)[sym_idx[1]]
message("  Using '", sym_col, "' as gene symbol source.")

fdat_df <- fdat_df %>%
  mutate(
    gene_symbol_raw = as.character(.data[[sym_col]]),
    gene_symbol     = gene_symbol_raw
  ) %>%
  mutate(
    ## اگر چند سمبل با ///, //, ;, , جدا شده‌اند، اولی را نگه می‌داریم
    gene_symbol = stringr::str_split(
      gene_symbol,
      pattern = "///|//|;|,",
      simplify = TRUE
    )[, 1],
    gene_symbol = stringr::str_trim(gene_symbol),
    gene_symbol = dplyr::if_else(
      gene_symbol %in% c("", "---", "NA", "NaN"),
      NA_character_,
      gene_symbol
    )
  )

n_non_na <- sum(!is.na(fdat_df$gene_symbol))
message("  Probes with non-missing gene_symbol: ",
        n_non_na, " / ", nrow(fdat_df))

if (n_non_na < 1000L) {
  warning("Number of probes with non-NA gene_symbol is unexpectedly small (",
          n_non_na, "). Annotation may be problematic.")
}

## ------------------------------------------------------------
## 6. collapse probes to genes
## ------------------------------------------------------------
message("\n--- Step 3: Collapse probe-level expression to gene_symbol ---")

expr_df <- expr_mat %>%
  as.data.frame() %>%
  rownames_to_column("probe_id")

expr_annot <- expr_df %>%
  left_join(
    fdat_df %>% select(probe_id, gene_symbol),
    by = "probe_id"
  ) %>%
  filter(!is.na(gene_symbol))

message("  Rows with gene_symbol after join: ", nrow(expr_annot))

if (nrow(expr_annot) == 0L) {
  stop("No rows with valid gene_symbol after join; cannot proceed.")
}

expr_by_gene <- expr_annot %>%
  group_by(gene_symbol) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = mean,
      .names = "{.col}"
    ),
    .groups = "drop"
  )

expr_mat_gene <- expr_by_gene %>%
  as.data.frame()

rownames(expr_mat_gene) <- expr_mat_gene$gene_symbol
expr_mat_gene$gene_symbol <- NULL
expr_mat_gene <- as.matrix(expr_mat_gene)

message("  Collapsed expression matrix: ",
        nrow(expr_mat_gene), " genes × ",
        ncol(expr_mat_gene), " samples")

## ------------------------------------------------------------
## 7. save expression matrix
## ------------------------------------------------------------
message("\n--- Step 4: Save gene-level expression matrix ---")

expr_out_df <- expr_mat_gene %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol")

readr::write_tsv(expr_out_df, expr_out_path)
message("  Saved expression (gene-level, log2) to:\n  ", expr_out_path)

## ------------------------------------------------------------
## 8. save sample metadata
## ------------------------------------------------------------
message("\n--- Step 5: Save sample metadata (raw pheno data) ---")

pheno_df <- pheno %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

readr::write_csv(pheno_df, sample_out_path)
message("  Saved sample metadata to:\n  ", sample_out_path)

## ------------------------------------------------------------
## 9. save feature annotation
## ------------------------------------------------------------
message("\n--- Step 6: Save feature annotation (probe-level) ---")

feature_out_df <- fdat_df %>%
  select(probe_id, gene_symbol_raw, gene_symbol, everything())

readr::write_csv(feature_out_df, feature_out_path)
message("  Saved feature annotation to:\n  ", feature_out_path)

## ------------------------------------------------------------
## 10. build SummarizedExperiment
## ------------------------------------------------------------
message("\n--- Step 7: Build SummarizedExperiment object ---")

row_data <- S4Vectors::DataFrame(
  gene_symbol = rownames(expr_mat_gene),
  row.names   = rownames(expr_mat_gene)
)

col_data <- S4Vectors::DataFrame(
  pheno_df[, !(names(pheno_df) %in% "sample_id"), drop = FALSE],
  row.names = pheno_df$sample_id
)

if (!identical(colnames(expr_mat_gene), rownames(col_data))) {
  warning("Column names of expression matrix do not match rownames of colData; aligning...")
  common <- intersect(colnames(expr_mat_gene), rownames(col_data))
  expr_mat_gene <- expr_mat_gene[, common, drop = FALSE]
  col_data      <- col_data[common, , drop = FALSE]
}

se <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(expr = expr_mat_gene),
  rowData = row_data,
  colData = col_data
)

saveRDS(se, se_out_path)
message("  Saved SummarizedExperiment to:\n  ", se_out_path)

message("\n=== Phase 9A (GSE31312 processing) completed successfully. ===")
