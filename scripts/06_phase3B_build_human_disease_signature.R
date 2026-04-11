## ==========================================
## 06_phase3B_build_human_disease_signature.R
## Phase 3B - Build human DLBCL disease signature
## From: DE_GSE56315_tumor_vs_normal
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"

de_dir       <- file.path(project_root, "results", "tables", "DE")
sig_dir      <- file.path(project_root, "results", "tables", "signatures")
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

## ---------- 0) Packages ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkg_cran <- c("tidyverse")

for (p in pkg_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---------- 1) پارامترهای امضا (می‌توانی بعداً عوض کنی) ----------

# حد آستانه‌ها برای “core disease signature”
logFC_cut <- 1       # |log2FC| >= 1
padj_cut  <- 0.05    # adj.P.Val <= 0.05

message("Using thresholds: |log2FC| >= ", logFC_cut,
        " and adj.P.Val <= ", padj_cut)

## ---------- 2) خواندن نتایج DE از GSE56315 ----------

de_file <- file.path(de_dir, "DE_GSE56315_tumor_vs_normal.tsv")

if (!file.exists(de_file)) {
  stop("DE file not found: ", de_file,
       "\nRun 05_phase3A_DE_human_tumor_normal.R first.")
}

de <- readr::read_tsv(de_file, show_col_types = FALSE)

required_cols <- c("gene_symbol", "logFC", "adj.P.Val")
missing_cols <- setdiff(required_cols, colnames(de))
if (length(missing_cols) > 0) {
  stop("DE file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

# حذف ژن‌های بدون اسم
de <- de %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "")

message("Loaded DE_GSE56315: ", nrow(de), " genes.")

## ---------- 3) تعیین جهت (up/down) و فیلتر بر اساس آستانه ----------

de_sig <- de %>%
  dplyr::filter(!is.na(adj.P.Val)) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      logFC >=  logFC_cut  & adj.P.Val <= padj_cut ~ "up",
      logFC <= -logFC_cut  & adj.P.Val <= padj_cut ~ "down",
      TRUE ~ "nonsig"
    )
  )

table_dir <- table(de_sig$direction, useNA = "ifany")
message("Direction summary (with thresholds):")
print(table_dir)

# فقط ژن‌های امضای بیماری (up/down) را نگه می‌داریم
sig_core <- de_sig %>%
  dplyr::filter(direction %in% c("up", "down")) %>%
  dplyr::arrange(dplyr::desc(abs(logFC)))

n_up   <- sum(sig_core$direction == "up")
n_down <- sum(sig_core$direction == "down")

message("Core signature size: ", nrow(sig_core),
        " genes (", n_up, " up, ", n_down, " down).")

if (nrow(sig_core) == 0) {
  warning("No genes passed the thresholds for core signature. ",
          "Consider relaxing logFC_cut or padj_cut.")
}

## ---------- 4) ذخیره خروجی‌ها ----------

# 4.1) جدول کامل امضای core
sig_core_path <- file.path(sig_dir, "human_disease_signature_GSE56315_core.tsv")
readr::write_tsv(sig_core, sig_core_path)

# 4.2) فقط لیست ژن‌های up و down (برای CMap / LINCS و غیره)
up_genes <- sig_core %>%
  dplyr::filter(direction == "up") %>%
  dplyr::pull(gene_symbol) %>%
  unique()

down_genes <- sig_core %>%
  dplyr::filter(direction == "down") %>%
  dplyr::pull(gene_symbol) %>%
  unique()

up_path   <- file.path(sig_dir, "human_disease_up_GSE56315_core.txt")
down_path <- file.path(sig_dir, "human_disease_down_GSE56315_core.txt")

readr::write_lines(up_genes,   up_path)
readr::write_lines(down_genes, down_path)

## ---------- 5) خلاصه روی کنسول ----------

cat("\n[Phase 3B] Human DLBCL disease signature built from GSE56315.\n")
cat("Core signature thresholds: |log2FC| >= ", logFC_cut,
    " & adj.P.Val <= ", padj_cut, "\n", sep = "")
cat("Total core genes: ", nrow(sig_core),
    " (", n_up, " up, ", n_down, " down)\n", sep = "")
cat("Core table: ", sig_core_path, "\n", sep = "")
cat("Up genes list: ", up_path, "\n", sep = "")
cat("Down genes list: ", down_path, "\n\n", sep = "")
source("D:/Research/My Articles/DLBCL drug/scripts/06_phase3B_build_human_disease_signature.R")
