## =====================================================================
## Phase 2D: Canine microarray processing (GSE30881: DLBCL vs healthy LN)
## =====================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
})

## ---------------------------------------------------------------------
## 0. پروژه
## ---------------------------------------------------------------------

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) {
  file.path(project_root, ...)
}

dir.create(path_proj("data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_proj("metadata"), recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1. تابع کمکی: collapse_duplicated_rows
## ---------------------------------------------------------------------

collapse_duplicated_rows <- function(mat) {
  if (!any(duplicated(rownames(mat)))) return(mat)
  
  idx_list <- split(seq_len(nrow(mat)), rownames(mat))
  out <- vapply(
    idx_list,
    function(idx) colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(mat))
  )
  out <- t(out)
  rownames(out) <- names(idx_list)
  out
}

## ---------------------------------------------------------------------
## 2. دانلود و آماده‌سازی GSE30881
## ---------------------------------------------------------------------

message("=== [GSE30881] Downloading / loading from GEO ===")

options('download.file.method.GEOquery' = 'auto')
options('GEOquery.inmemory.gpl' = FALSE)

gse_list <- GEOquery::getGEO("GSE30881", GSEMatrix = TRUE)
es <- if (is.list(gse_list)) gse_list[[1]] else gse_list

expr_raw <- Biobase::exprs(es)
rng <- range(expr_raw, finite = TRUE)
message("[GSE30881] Raw expr range: [", signif(rng[1], 3), ", ", signif(rng[2], 3), "]")

## فرض: اگر max خیلی بالا باشد، log2 نشده → log2(x+1)
if (rng[2] > 50) {
  message("[GSE30881] Expression looks non-log. Applying log2(x + 1).")
  expr_log2 <- log2(expr_raw + 1)
} else {
  message("[GSE30881] Expression appears already log-scale. Using as-is.")
  expr_log2 <- expr_raw
}

## ---------------------------------------------------------------------
## 3. Feature annotation و gene_symbol
## ---------------------------------------------------------------------

fdat <- Biobase::fData(es) %>% as.data.frame()
feature_anno <- fdat %>%
  rownames_to_column("probe_id")

symbol_candidates <- c(
  "gene_symbol", "GENE_SYMBOL", "Gene Symbol", "Gene.symbol",
  "Symbol", "SYMBOL", "geneSymbol"
)
sym_col <- intersect(symbol_candidates, names(feature_anno))
if (length(sym_col) > 0L) {
  sym_col <- sym_col[1L]
  feature_anno$gene_symbol <- as.character(feature_anno[[sym_col]])
  message("[GSE30881] Using '", sym_col, "' as gene_symbol.")
} else {
  feature_anno$gene_symbol <- NA_character_
  warning("[GSE30881] No obvious gene_symbol column in fData. gene_symbol will be NA.")
}

feature_anno_path <- path_proj("metadata", "GSE30881_feature_annotation.csv")
readr::write_csv(feature_anno, feature_anno_path)
message("[GSE30881] Feature annotation written to: ", feature_anno_path)

## ---------------------------------------------------------------------
## 4. ساخت ماتریس بر پایه gene_symbol
## ---------------------------------------------------------------------

gs <- feature_anno$gene_symbol
valid <- !is.na(gs) & gs != ""
expr_sub <- expr_log2[valid, , drop = FALSE]
gs_valid <- gs[valid]
rownames(expr_sub) <- gs_valid

expr_gene <- collapse_duplicated_rows(expr_sub)

message("[GSE30881] Gene-level matrix: ", nrow(expr_gene),
        " genes × ", ncol(expr_gene), " samples.")

## ---------------------------------------------------------------------
## 5. متادیتای نمونه‌ها + group_tn (tumor / normal)
## ---------------------------------------------------------------------

pd <- Biobase::pData(es) %>% as.data.frame() %>%
  rownames_to_column("sample_id")

char_cols <- names(pd)[grepl("characteristics", names(pd), ignore.case = TRUE)]
text_cols <- unique(c("title", "source_name_ch1", char_cols))
text_cols <- intersect(text_cols, names(pd))

if (length(text_cols) == 0L) {
  warning("[GSE30881] No descriptive columns found to infer group_tn. group_tn will be NA.")
  pd$group_tn <- NA_character_
} else {
  text_all <- apply(pd[, text_cols, drop = FALSE], 1L, function(x) {
    paste(x, collapse = " ; ")
  })
  text_low <- tolower(text_all)
  
  is_normal <- stringr::str_detect(
    text_low,
    "normal|healthy|control|non-lymphoma|non lymphoma|non-lymphoid"
  )
  is_tumor <- stringr::str_detect(
    text_low,
    "lymphoma|dlbcl|diffuse large b cell|diffuse large b-cell|b[- ]cell lymphoma|b cell lymphoma|cancer|tumor|tumour|malignant"
  )
  
  group_tn <- ifelse(is_tumor & !is_normal, "tumor",
                     ifelse(is_normal & !is_tumor, "normal", NA_character_))
  pd$group_tn <- group_tn
}

message("[GSE30881] group_tn counts (including NA):")
print(table(pd$group_tn, useNA = "ifany"))

## متادیتا فقط برای sampleهایی که در expr_gene هستند
pd <- pd %>%
  filter(sample_id %in% colnames(expr_gene)) %>%
  distinct(sample_id, .keep_all = TRUE)

meta_path <- path_proj("metadata", "GSE30881_sample_metadata.csv")
readr::write_csv(pd, meta_path)
message("[GSE30881] Sample metadata written to: ", meta_path)

## ---------------------------------------------------------------------
## 6. ذخیره‌ی expr + SE
## ---------------------------------------------------------------------

expr_out <- expr_gene %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol")

expr_path <- path_proj("data", "processed", "GSE30881_expr_log2_qcfiltered.tsv")
readr::write_tsv(expr_out, expr_path)
message("[GSE30881] Expression matrix written to: ", expr_path)

## SummarizedExperiment
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = expr_gene),
  colData = S4Vectors::DataFrame(pd),
  rowData = S4Vectors::DataFrame(
    tibble(gene_symbol = rownames(expr_gene))
  )
)

se_path <- path_proj("data", "processed", "GSE30881_SummarizedExperiment.rds")
saveRDS(se, se_path)
message("[GSE30881] SummarizedExperiment written to: ", se_path)

message("=== [Phase 2D - GSE30881 canine microarray processing] DONE ===")
