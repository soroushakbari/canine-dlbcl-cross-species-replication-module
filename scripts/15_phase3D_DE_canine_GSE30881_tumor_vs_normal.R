## =====================================================================
## Phase 3D: Canine DE analysis (GSE30881: DLBCL vs healthy LN)
## =====================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) {
  file.path(project_root, ...)
}

dir.create(path_proj("results", "tables", "DE"),
           recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1. خواندن expression gene-level
## ---------------------------------------------------------------------

expr_path <- path_proj("data", "processed", "GSE30881_expr_log2_qcfiltered.tsv")
if (!file.exists(expr_path)) {
  stop("Expression file not found: ", expr_path,
       "\nRun 02_phase2D_microarray_canine_GSE30881.R first.")
}

expr_df <- readr::read_tsv(expr_path, show_col_types = FALSE)
if (!"gene_symbol" %in% names(expr_df)) {
  stop("Expression file must have a 'gene_symbol' column as first column.")
}

gene_symbols <- expr_df$gene_symbol %>% as.character()
expr_mat <- expr_df %>%
  select(-gene_symbol) %>%
  as.matrix()
rownames(expr_mat) <- gene_symbols

message("[GSE30881] Expression matrix loaded: ",
        nrow(expr_mat), " genes × ", ncol(expr_mat), " samples.")

## ---------------------------------------------------------------------
## 2. استخراج group_tn از GEO
## ---------------------------------------------------------------------

infer_group_tn_from_geo <- function(gse_id, valid_sample_ids = NULL) {
  message("[", gse_id, "] Inferring group_tn (tumor/normal) from GEO...")
  
  options('download.file.method.GEOquery' = 'auto')
  options('GEOquery.inmemory.gpl' = FALSE)
  
  gse_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  es <- if (is.list(gse_list)) gse_list[[1]] else gse_list
  
  pd <- Biobase::pData(es) %>% as.data.frame()
  pd$sample_id <- rownames(pd)
  
  if (!is.null(valid_sample_ids)) {
    pd <- pd %>% dplyr::filter(sample_id %in% valid_sample_ids)
  }
  
  if (nrow(pd) == 0L) {
    stop("[", gse_id, "] No overlapping samples between GEO and expression matrix.")
  }
  
  char_cols <- names(pd)[grepl("characteristics", names(pd), ignore.case = TRUE)]
  text_cols <- unique(c("title", "source_name_ch1", char_cols))
  text_cols <- intersect(text_cols, names(pd))
  
  if (length(text_cols) == 0L) {
    stop("[", gse_id, "] No descriptive columns found to infer group_tn.")
  }
  
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
  
  out <- tibble(
    sample_id = pd$sample_id,
    group_tn  = group_tn
  )
  
  message("[", gse_id, "] Inferred group_tn counts (including NA):")
  print(table(out$group_tn, useNA = "ifany"))
  
  ## --- Fallback مهم برای این دیتاست ---
  ## اگر فقط 'tumor' و NA داریم، NAها را 'normal' در نظر می‌گیریم
  tab <- table(out$group_tn, useNA = "ifany")
  if ("tumor" %in% names(tab) &&
      !("normal" %in% names(tab)) &&
      any(is.na(out$group_tn))) {
    message("[", gse_id, "] Fallback: only 'tumor' and NA detected. ",
            "Treating NA as 'normal' (healthy LN).")
    out$group_tn[is.na(out$group_tn)] <- "normal"
    message("[", gse_id, "] group_tn counts AFTER fallback:")
    print(table(out$group_tn, useNA = "ifany"))
  }
  ## --------------------------------------
  
  out
}

gt <- infer_group_tn_from_geo("GSE30881", valid_sample_ids = colnames(expr_mat))

## ---------------------------------------------------------------------
## 3. آماده‌سازی برای limma (فقط tumor و normal)
## ---------------------------------------------------------------------

merged <- tibble(sample_id = colnames(expr_mat)) %>%
  left_join(gt, by = "sample_id")

tab_all <- table(merged$group_tn, useNA = "ifany")
message("[GSE30881] group_tn in merged data (including NA):")
print(tab_all)

dn <- merged %>%
  filter(group_tn %in% c("tumor", "normal")) %>%
  mutate(group_tn = as.character(group_tn))

tab_used <- table(dn$group_tn)
message("[GSE30881] group_tn used for DE:")
print(tab_used)

if (!all(c("tumor", "normal") %in% names(tab_used))) {
  stop("[GSE30881] Both 'tumor' and 'normal' are required for DE but not found.")
}

expr_mat_used <- expr_mat[, dn$sample_id, drop = FALSE]

group <- factor(dn$group_tn, levels = c("normal", "tumor"))
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")

message("[GSE30881] Design matrix (first rows):")
print(design[1:min(5, nrow(design)), , drop = FALSE])

## ---------------------------------------------------------------------
## 4. limma: tumor vs normal
## ---------------------------------------------------------------------

fit <- limma::lmFit(expr_mat_used, design)
fit <- limma::eBayes(fit)

tt <- limma::topTable(
  fit,
  coef    = "Tumor_vs_Normal",
  number  = Inf,
  sort.by = "none"
)

tt <- tt %>%
  rownames_to_column("gene_symbol")

## خلاصه‌ی جهت‌ها با threshold معقول
is_sig <- (abs(tt$logFC) >= 1) & (tt$adj.P.Val <= 0.05)
dir_sig <- ifelse(tt$logFC > 0, "up",
                  ifelse(tt$logFC < 0, "down", "nonsig"))
summary_table <- table(direction = dir_sig[is_sig])

message("[GSE30881] Significant genes (|log2FC| >= 1 & adj.P.Val <= 0.05):")
print(summary_table)

## ---------------------------------------------------------------------
## 5. ذخیره‌ی نتایج
## ---------------------------------------------------------------------

de_path <- path_proj("results", "tables", "DE",
                     "DE_GSE30881_tumor_vs_normal.tsv")
readr::write_tsv(tt, de_path)
message("[GSE30881] DE table written to: ", de_path)

message("=== [Phase 3D - GSE30881 DE (tumor vs normal)] DONE ===")
