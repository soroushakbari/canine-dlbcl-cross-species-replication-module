## ==========================================
## 07_phase3C_refine_human_signature_multi_tier.R
## Phase 3C - Refine human DLBCL disease signature (multi-tier)
## Input: human_disease_signature_GSE56315_core.tsv
## Output: Tier1 (core) و Tier2 (broad) امضاهای up/down
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"

sig_dir_in  <- file.path(project_root, "results", "tables", "signatures")
sig_dir_out <- file.path(project_root, "results", "tables", "signatures")
dir.create(sig_dir_out, recursive = TRUE, showWarnings = FALSE)

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

## ---------- 1) خواندن core signature ----------

core_file <- file.path(sig_dir_in, "human_disease_signature_GSE56315_core.tsv")

if (!file.exists(core_file)) {
  stop("Core signature file not found: ", core_file,
       "\nRun 06_phase3B_build_human_disease_signature.R first.")
}

core <- readr::read_tsv(core_file, show_col_types = FALSE)

required_cols <- c("gene_symbol", "logFC", "adj.P.Val", "direction")
missing_cols <- setdiff(required_cols, colnames(core))
if (length(missing_cols) > 0) {
  stop("Core signature file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

core <- core %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::filter(direction %in% c("up", "down"))

message("Loaded core signature: ", nrow(core),
        " genes (", sum(core$direction == "up"), " up, ",
        sum(core$direction == "down"), " down).")

## ---------- 2) تعریف score برای رتبه‌بندی ----------

core <- core %>%
  dplyr::mutate(
    abs_logFC   = abs(logFC),
    neglog10padj = -log10(adj.P.Val + 1e-300),   # محافظت مقابل صفر
    rank_score  = abs_logFC * neglog10padj
  )

## ---------- 3) پارامترهای Tier ها ----------

# Tier1: signature سخت و فشرده برای CMap/PPI/هسته‌ی cross-species
tier1_logFC_cut <- 1.5
tier1_padj_cut  <- 0.01
tier1_max_up    <- 300   # حداکثر تعداد up
tier1_max_down  <- 300   # حداکثر تعداد down

# Tier2: signature گسترده‌تر برای GSVA/Enrichment
tier2_logFC_cut <- 1.0
tier2_padj_cut  <- 0.05
tier2_max_up    <- 1000
tier2_max_down  <- 1000

## ---------- 4) Helper برای ساخت یک Tier ----------

build_tier <- function(df, logFC_cut, padj_cut, max_up, max_down, tier_name = "TierX") {
  
  message("\nBuilding ", tier_name,
          " with thresholds: |log2FC| >= ", logFC_cut,
          ", adj.P.Val <= ", padj_cut,
          " and caps: up <= ", max_up,
          ", down <= ", max_down)
  
  df_sub <- df %>%
    dplyr::filter(
      abs_logFC >= logFC_cut,
      adj.P.Val <= padj_cut
    )
  
  if (nrow(df_sub) == 0) {
    warning(tier_name, ": no genes pass thresholds.")
    return(list(
      table = df_sub,
      up    = character(0),
      down  = character(0)
    ))
  }
  
  up_tbl <- df_sub %>%
    dplyr::filter(direction == "up") %>%
    dplyr::arrange(dplyr::desc(rank_score))
  
  down_tbl <- df_sub %>%
    dplyr::filter(direction == "down") %>%
    dplyr::arrange(dplyr::desc(rank_score))
  
  if (nrow(up_tbl) > max_up) {
    up_tbl <- up_tbl[seq_len(max_up), , drop = FALSE]
  }
  if (nrow(down_tbl) > max_down) {
    down_tbl <- down_tbl[seq_len(max_down), , drop = FALSE]
  }
  
  tier_tbl <- dplyr::bind_rows(up_tbl, down_tbl) %>%
    dplyr::arrange(direction, dplyr::desc(rank_score))
  
  n_up   <- sum(tier_tbl$direction == "up")
  n_down <- sum(tier_tbl$direction == "down")
  
  message(tier_name, " size: ", nrow(tier_tbl),
          " genes (", n_up, " up, ", n_down, " down).")
  
  list(
    table = tier_tbl,
    up    = unique(up_tbl$gene_symbol),
    down  = unique(down_tbl$gene_symbol)
  )
}

## ---------- 5) ساخت Tier1 و Tier2 ----------

tier1 <- build_tier(
  df         = core,
  logFC_cut  = tier1_logFC_cut,
  padj_cut   = tier1_padj_cut,
  max_up     = tier1_max_up,
  max_down   = tier1_max_down,
  tier_name  = "Tier1_core"
)

tier2 <- build_tier(
  df         = core,
  logFC_cut  = tier2_logFC_cut,
  padj_cut   = tier2_padj_cut,
  max_up     = tier2_max_up,
  max_down   = tier2_max_down,
  tier_name  = "Tier2_broad"
)

## ---------- 6) ذخیره خروجی‌ها ----------

# Tier1
tier1_tbl <- tier1$table
tier1_tbl_path <- file.path(sig_dir_out, "human_signature_GSE56315_Tier1_core.tsv")
readr::write_tsv(tier1_tbl, tier1_tbl_path)

tier1_up_path   <- file.path(sig_dir_out, "human_signature_GSE56315_Tier1_core_up.txt")
tier1_down_path <- file.path(sig_dir_out, "human_signature_GSE56315_Tier1_core_down.txt")

readr::write_lines(tier1$up,   tier1_up_path)
readr::write_lines(tier1$down, tier1_down_path)

# Tier2
tier2_tbl <- tier2$table
tier2_tbl_path <- file.path(sig_dir_out, "human_signature_GSE56315_Tier2_broad.tsv")
readr::write_tsv(tier2_tbl, tier2_tbl_path)

tier2_up_path   <- file.path(sig_dir_out, "human_signature_GSE56315_Tier2_broad_up.txt")
tier2_down_path <- file.path(sig_dir_out, "human_signature_GSE56315_Tier2_broad_down.txt")

readr::write_lines(tier2$up,   tier2_up_path)
readr::write_lines(tier2$down, tier2_down_path)

## ---------- 7) خلاصه ----------

cat("\n[Phase 3C] Refined human DLBCL signature (multi-tier) from core GSE56315.\n")

cat("Tier1_core (strict, for CMap/PPI/cross-species core):\n",
    "  Table: ", tier1_tbl_path, "\n",
    "  Up genes list: ", tier1_up_path, "\n",
    "  Down genes list: ", tier1_down_path, "\n", sep = "")

cat("\nTier2_broad (for GSVA/enrichment):\n",
    "  Table: ", tier2_tbl_path, "\n",
    "  Up genes list: ", tier2_up_path, "\n",
    "  Down genes list: ", tier2_down_path, "\n\n", sep = "")
