## ==========================================
## 05_phase3A_DE_human_tumor_normal.R
## Phase 3A - Differential expression (tumor vs normal)
## Human DLBCL: GSE56315, GSE32018
## ==========================================

## ---------- 0) Paths ----------

project_root <- "D:/Research/My Articles/DLBCL drug"

processed_dir <- file.path(project_root, "data", "processed")
metadata_dir  <- file.path(project_root, "metadata")   # فقط اگر خواستی بعداً ازش استفاده کنی
results_dir   <- file.path(project_root, "results", "tables", "DE")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

## ---------- 1) Packages ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkg_bioc <- c("limma", "GEOquery", "Biobase")
pkg_cran <- c("tidyverse")

for (p in pkg_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

for (p in pkg_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(GEOquery)
  library(Biobase)
})

## ---------- 2) Load expression matrix (log2, QC-filtered) ----------

load_expr_geo <- function(gse_id) {
  message("\n==============================")
  message("Loading expression for: ", gse_id)
  message("==============================")
  
  expr_path <- file.path(processed_dir, paste0(gse_id, "_expr_log2_qcfiltered.tsv"))
  if (!file.exists(expr_path)) {
    stop("Expression file not found: ", expr_path)
  }
  
  expr_tbl <- readr::read_tsv(expr_path, show_col_types = FALSE)
  
  if (!"gene_symbol" %in% colnames(expr_tbl)) {
    stop("Expression table for ", gse_id, " has no 'gene_symbol' column.")
  }
  
  expr_mat <- expr_tbl %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%  # اگر ژن تکراریه، اولینش رو نگه می‌داریم
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix()
  
  message("Expression matrix: ",
          nrow(expr_mat), " genes x ",
          ncol(expr_mat), " samples.")
  
  expr_mat
}

## ---------- 3) DE helper (limma) ----------

run_limma_tn <- function(expr_mat, meta_df, gse_id) {
  # فقط نمونه‌هایی که group_tn مشخص (tumor/normal) دارند
  keep <- !is.na(meta_df$group_tn) & meta_df$group_tn %in% c("tumor", "normal")
  
  expr_sub <- expr_mat[, keep, drop = FALSE]
  meta_sub <- meta_df[keep, , drop = FALSE]
  
  n_tumor  <- sum(meta_sub$group_tn == "tumor")
  n_normal <- sum(meta_sub$group_tn == "normal")
  
  message("DE for ", gse_id, ": using ", ncol(expr_sub),
          " samples (", n_tumor, " tumor, ", n_normal, " normal).")
  
  if (length(unique(meta_sub$group_tn)) != 2) {
    stop("For ", gse_id, ": need both tumor and normal, but got only: ",
         paste(unique(meta_sub$group_tn), collapse = ", "))
  }
  
  group_factor <- factor(meta_sub$group_tn, levels = c("normal", "tumor"))
  design <- model.matrix(~ group_factor)
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  
  # coef=2 => tumor vs normal
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("gene_symbol") %>%
    dplyr::rename(
      logFC     = logFC,
      AveExpr   = AveExpr,
      t_stat    = t,
      P.Value   = P.Value,
      adj.P.Val = adj.P.Val,
      B_stat    = B
    )
  
  tt
}

## ---------- 4) Helper: generic GEO-based mapping (sample_ids ↔ metadata) ----------

map_group_from_geo <- function(meta_geo, sample_ids, group_col = "group_tn", gse_id = "GSE") {
  meta_geo$gsm_id <- rownames(meta_geo)
  
  # تمام ستون‌ها (gsm_id + بقیه) رو برای مپ بررسی می‌کنیم
  all_cols <- unique(c("gsm_id", colnames(meta_geo)))
  best_col   <- NA_character_
  best_match <- 0L
  
  for (col in all_cols) {
    vals <- as.character(meta_geo[[col]])
    m <- sum(vals %in% sample_ids)
    if (!is.na(m) && m > best_match) {
      best_match <- m
      best_col   <- col
    }
  }
  
  if (is.na(best_col) || best_match == 0) {
    warning(gse_id, ": could not match sample_ids to any GEO metadata column. ",
            "Returning minimal metadata with NA group_tn.")
    meta_df <- tibble::tibble(
      sample_id = sample_ids,
      group_tn  = NA_character_
    )
    return(meta_df)
  }
  
  message(gse_id, ": best mapping column in GEO metadata: ", best_col,
          " (", best_match, " / ", length(sample_ids), " samples matched).")
  
  meta_map <- tibble::tibble(
    map_val  = as.character(meta_geo[[best_col]]),
    group_tn = meta_geo[[group_col]]
  )
  
  meta_df <- tibble::tibble(sample_id = sample_ids) %>%
    dplyr::left_join(meta_map, by = c("sample_id" = "map_val"))
  
  message(gse_id, ": final group_tn summary after mapping:")
  print(table(meta_df$group_tn, useNA = "ifany"))
  
  meta_df
}

## ---------- 5) GSE56315: build tumor/normal from GEO ----------

build_meta_gse56315 <- function(sample_ids) {
  gse_id <- "GSE56315"
  message("\n[", gse_id, "] Fetching GEO metadata and defining tumor/normal ...")
  
  gse_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  if (length(gse_list) > 1) {
    ns <- sapply(gse_list, function(es) ncol(Biobase::exprs(es)))
    idx <- which.max(ns)
    gset <- gse_list[[idx]]
  } else {
    gset <- gse_list[[1]]
  }
  
  meta_geo <- Biobase::pData(gset) %>%
    as.data.frame()
  
  # همه‌ی ستون‌های کاراکتری را در هم می‌ریزیم
  char_cols <- meta_geo %>%
    dplyr::select(where(is.character))
  
  if (ncol(char_cols) == 0) {
    warning(gse_id, ": GEO metadata has no character columns; cannot infer groups.")
    return(tibble::tibble(sample_id = sample_ids,
                          group_tn  = NA_character_))
  }
  
  text_all <- apply(char_cols, 1, function(x) paste(x, collapse = " | ")) %>%
    tolower()
  
  # طراحی این دیتاست مشخص است: 33 tonsil normal vs 55 DLBCL
  is_normal <- grepl("tonsil|reactive tonsil|normal tonsil|control tonsil",
                     text_all, ignore.case = TRUE)
  
  # هرچیزی که tonsil نیست، DLBCL/lymphoma است
  group_tn_geo <- ifelse(is_normal, "normal", "tumor")
  
  meta_geo$group_tn <- group_tn_geo
  
  message(gse_id, ": GEO-level group_tn table:")
  print(table(meta_geo$group_tn, useNA = "ifany"))
  
  # مپ به sample_ids
  meta_df <- map_group_from_geo(meta_geo, sample_ids, group_col = "group_tn", gse_id = gse_id)
  
  # sanity check: ideally 33 normal, 55 tumor
  tab <- table(meta_df$group_tn, useNA = "ifany")
  message(gse_id, ": mapped group_tn table (should be ~33 normal, 55 tumor):")
  print(tab)
  
  meta_df
}

## ---------- 6) GSE32018: build tumor/normal from GEO ----------

build_meta_gse32018 <- function(sample_ids) {
  gse_id <- "GSE32018"
  message("\n[", gse_id, "] Fetching GEO metadata and defining tumor/normal ...")
  
  gse_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  if (length(gse_list) > 1) {
    ns <- sapply(gse_list, function(es) ncol(Biobase::exprs(es)))
    idx <- which.max(ns)
    gset <- gse_list[[idx]]
  } else {
    gset <- gse_list[[1]]
  }
  
  meta_geo <- Biobase::pData(gset) %>%
    as.data.frame()
  
  char_cols <- meta_geo %>%
    dplyr::select(where(is.character))
  
  if (ncol(char_cols) == 0) {
    warning(gse_id, ": GEO metadata has no character columns; cannot infer groups.")
    return(tibble::tibble(sample_id = sample_ids,
                          group_tn  = NA_character_))
  }
  
  text_all <- apply(char_cols, 1, function(x) paste(x, collapse = " | ")) %>%
    tolower()
  
  # نرمال / واکنشی
  is_normal <- grepl(
    "reactive tonsil|reactive lymph node|reactive lymph nodes|reactive tonsils|normal lymph node|normal tonsil|control",
    text_all,
    ignore.case = TRUE
  )
  
  # تومور (انواع لنفوم و CLL)
  is_tumor <- grepl(
    "lymphoma|dlbcl|diffuse large b|mantle cell|mcl|follicular|malt|mzl|nmzl|cll|non[- ]hodgkin",
    text_all,
    ignore.case = TRUE
  )
  
  group_tn_geo <- ifelse(is_normal & !is_tumor, "normal",
                         ifelse(is_tumor & !is_normal, "tumor", NA))
  
  meta_geo$group_tn <- group_tn_geo
  
  message(gse_id, ": GEO-level group_tn table:")
  print(table(meta_geo$group_tn, useNA = "ifany"))
  
  meta_df <- map_group_from_geo(meta_geo, sample_ids, group_col = "group_tn", gse_id = gse_id)
  
  meta_df
}

## ---------- 7) Main loop over datasets ----------

datasets_tn <- c("GSE56315", "GSE32018")

for (gse_id in datasets_tn) {
  expr_mat <- load_expr_geo(gse_id)
  sample_ids <- colnames(expr_mat)
  
  if (gse_id == "GSE56315") {
    meta_df <- build_meta_gse56315(sample_ids)
  } else if (gse_id == "GSE32018") {
    meta_df <- build_meta_gse32018(sample_ids)
  } else {
    # فعلاً دیتاست دیگری نداریم؛ اگر خواستی بعداً generic logic اضافه کن
    meta_df <- tibble::tibble(sample_id = sample_ids,
                              group_tn  = NA_character_)
  }
  
  # قبل از limma چک می‌کنیم اصلاً دو گروه داریم یا نه
  group_vals <- meta_df$group_tn[!is.na(meta_df$group_tn)]
  unique_groups <- unique(group_vals)
  
  if (length(unique_groups) < 2) {
    warning("For ", gse_id, ": after GEO-based annotation we still do NOT have both tumor and normal.\n",
            "Skipping DE for this dataset.")
    next
  }
  
  # حالا DE را می‌زنیم
  tt <- run_limma_tn(expr_mat, meta_df, gse_id)
  
  out_path <- file.path(results_dir, paste0("DE_", gse_id, "_tumor_vs_normal.tsv"))
  readr::write_tsv(tt, out_path)
  
  message("Saved DE results for ", gse_id, " to:\n  ", out_path)
}

cat("\n[Phase 3A] Human tumor vs normal DE finished.\n")
cat("Results (where available) are under:\n  ", results_dir, "\n\n")
