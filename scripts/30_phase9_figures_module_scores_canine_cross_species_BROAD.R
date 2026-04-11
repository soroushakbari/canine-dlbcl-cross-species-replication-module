## 30_phase9_figures_module_scores_canine_cross_species_BROAD.R
## هدف:
##  - محاسبه‌ی نمره‌ی ماژول cross-species BROAD (سمت dog) در GSE30881
##  - ساخت شکل Fig. 6A: DLBCL vs normal LN در سگ
##  - استفاده فقط از فایل‌های موجود (expression + metadata + dog up/down)

message("=== Phase 9 (Figures): Canine cross-species BROAD module (GSE30881) – Fig. 6A ===")

required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr", "ggplot2"
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
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

expr_dir <- file.path(project_root, "data", "processed")
meta_dir <- file.path(project_root, "metadata")
sig_dir  <- file.path(project_root, "results", "tables", "signatures")
ms_dir   <- file.path(project_root, "results", "tables", "module_scores")
fig_dir  <- file.path(project_root, "results", "figures")

for (d in c(expr_dir, meta_dir, sig_dir)) {
  if (!dir.exists(d)) {
    stop("Expected directory not found:\n  ", d)
  }
}
if (!dir.exists(ms_dir)) {
  dir.create(ms_dir, recursive = TRUE)
  message("Created module_scores directory:\n  ", ms_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

expr_30881_path <- file.path(expr_dir, "GSE30881_expr_log2_qcfiltered.tsv")
meta_30881_path <- file.path(meta_dir, "GSE30881_sample_metadata.csv")

dog_up_path   <- file.path(sig_dir, "cross_species_module_BROAD_dog_up.txt")
dog_down_path <- file.path(sig_dir, "cross_species_module_BROAD_dog_down.txt")

if (!file.exists(expr_30881_path)) {
  stop("Expression file for GSE30881 not found at:\n  ", expr_30881_path)
}
if (!file.exists(meta_30881_path)) {
  stop("Metadata file for GSE30881 not found at:\n  ", meta_30881_path)
}
if (!file.exists(dog_up_path)) {
  stop("Dog up-signature file not found at:\n  ", dog_up_path)
}
if (!file.exists(dog_down_path)) {
  stop("Dog down-signature file not found at:\n  ", dog_down_path)
}

## --------------------------------------------------------------------
## 2. helper: read signature and compute module scores
## --------------------------------------------------------------------

read_sig_genes <- function(path) {
  x <- readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE)
  v <- unique(na.omit(trimws(as.character(unlist(x)))))
  v[v != ""]
}

# تلاش می‌کنیم ستون ژن را حدس بزنیم (اولین ستون یا اسم‌های رایج)
detect_gene_column <- function(df) {
  cand <- c("gene_symbol", "GeneSymbol", "Gene.symbol",
            "gene", "Gene", "id", "ID")
  present <- cand[cand %in% names(df)]
  if (length(present) > 0) {
    return(present[1])
  } else {
    # fallback: ستون اول
    return(names(df)[1])
  }
}

compute_module_scores_z <- function(expr_df, up_genes, down_genes,
                                    dataset_label = "GSE30881") {
  gene_col <- detect_gene_column(expr_df)
  message("  [", dataset_label, "] Using '", gene_col,
          "' as gene identifier column.")
  
  genes <- expr_df[[gene_col]]
  mat   <- expr_df[, setdiff(names(expr_df), gene_col), drop = FALSE]
  
  # تبدیل به ماتریس عددی
  mat_num <- as.matrix(mat)
  mode(mat_num) <- "numeric"
  
  rownames(mat_num) <- genes
  
  # Z-score روی هر ژن across samples
  z_mat <- t(scale(t(mat_num)))
  # اگر جایی فقط NA شد، scale عددی نمی‌دهد؛ آن‌ها را صفر می‌کنیم
  z_mat[is.na(z_mat)] <- 0
  
  up_present   <- intersect(rownames(z_mat), up_genes)
  down_present <- intersect(rownames(z_mat), down_genes)
  
  message("  [", dataset_label, "] Genes in expression matrix: ", nrow(z_mat))
  message("  [", dataset_label, "] Up genes in signature: ", length(up_genes),
          " ; present: ", length(up_present))
  message("  [", dataset_label, "] Down genes in signature: ", length(down_genes),
          " ; present: ", length(down_present))
  
  if (length(up_present) == 0 || length(down_present) == 0) {
    stop(
      "[", dataset_label, "] No overlap between expression genes and up/down signature.\n",
      "  up_present: ", length(up_present), " ; down_present: ", length(down_present)
    )
  }
  
  score_up   <- colMeans(z_mat[up_present, , drop = FALSE])
  score_down <- colMeans(z_mat[down_present, , drop = FALSE])
  score_net  <- score_up - score_down
  
  tibble(
    sample_id   = colnames(z_mat),
    score_up    = as.numeric(score_up),
    score_down  = as.numeric(score_down),
    score_net   = as.numeric(score_net),
    score_net_z = as.numeric(scale(score_net))
  )
}

## --------------------------------------------------------------------
## 3. load data and compute scores for GSE30881
## --------------------------------------------------------------------
message("--- Loading signatures (dog cross-species BROAD) ---")
dog_up   <- read_sig_genes(dog_up_path)
dog_down <- read_sig_genes(dog_down_path)

message("  Up-signature genes:   ", length(dog_up))
message("  Down-signature genes: ", length(dog_down))

message("\n--- Computing module scores for GSE30881 ---")
expr_30881 <- readr::read_tsv(expr_30881_path, show_col_types = FALSE)
meta_30881 <- readr::read_csv(meta_30881_path, show_col_types = FALSE)

# اطمینان از وجود ستون sample_id در متادیتا
if (!"sample_id" %in% names(meta_30881)) {
  if ("geo_accession" %in% names(meta_30881)) {
    meta_30881 <- meta_30881 %>%
      dplyr::rename(sample_id = geo_accession)
    message("  'sample_id' not found in metadata; using 'geo_accession' as sample_id.")
  } else {
    stop("Metadata for GSE30881 lacks 'sample_id' (or 'geo_accession') column.")
  }
}

scores_30881 <- compute_module_scores_z(
  expr_df  = expr_30881,
  up_genes = dog_up,
  down_genes = dog_down,
  dataset_label = "GSE30881"
)

# join با متادیتا
df_30881 <- scores_30881 %>%
  dplyr::left_join(meta_30881, by = "sample_id")

if (!"group_tn" %in% names(df_30881)) {
  stop("Combined GSE30881 table lacks 'group_tn' column (tumor/normal info).")
}

df_30881 <- df_30881 %>%
  dplyr::mutate(
    group_clean = tolower(as.character(group_tn)),
    group_clean = dplyr::case_when(
      group_clean %in% c("tumor", "dlbcl", "cancer", "case") ~ "DLBCL",
      is.na(group_clean)                                    ~ "Normal LN",
      group_clean %in% c("normal", "control", "healthy")    ~ "Normal LN",
      TRUE                                                  ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group_clean))

message("  GSE30881 group counts:")
message("    - Normal LN: ", sum(df_30881$group_clean == "Normal LN"))
message("    - DLBCL    : ", sum(df_30881$group_clean == "DLBCL"))

if (dplyr::n_distinct(df_30881$group_clean) < 2) {
  stop("After cleaning, fewer than two groups remain in GSE30881.")
}

## ذخیره‌ی جدول نمره‌ها
out_ms_30881 <- file.path(
  ms_dir,
  "module_scores_GSE30881_crossSpecies_BROAD_dog_zscore_with_group.tsv"
)
readr::write_tsv(df_30881, out_ms_30881)
message("Saved GSE30881 module-score table to:")
message("  ", out_ms_30881)

## --------------------------------------------------------------------
## 4. Fig. 6A – violin/boxplot for canine DLBCL vs normal LN
## --------------------------------------------------------------------
message("\n--- Building Fig. 6A: Canine DLBCL vs normal LN (GSE30881) ---")

p6A <- ggplot(df_30881, aes(x = group_clean, y = score_net_z, fill = group_clean)) +
  geom_violin(trim = FALSE, alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, colour = "grey20") +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.5, colour = "black") +
  scale_fill_manual(
    values = c("Normal LN" = "#74a9cf", "DLBCL" = "#fb6a4a"),
    guide  = "none"
  ) +
  labs(
    x = NULL,
    y = "Cross-species BROAD module score (z)",
    title = "GSE30881 (canine) – cross-species module activity in DLBCL vs normal LN"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 11, face = "bold"),
    axis.text.y      = element_text(size = 10),
    plot.title       = element_text(hjust = 0, face = "bold", size = 13),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )

# t-test برای annotation
tt_30881 <- t.test(score_net_z ~ group_clean, data = df_30881)
p_val_30881 <- tt_30881$p.value
p_lab_30881 <- if (p_val_30881 < 1e-4) {
  "p < 1e-4"
} else {
  paste0("p = ", signif(p_val_30881, 3))
}

p6A <- p6A +
  annotate(
    "text",
    x = 1.5,
    y = max(df_30881$score_net_z, na.rm = TRUE) * 1.05,
    label = p_lab_30881,
    size = 3.5
  )

out6A_png <- file.path(fig_dir, "Fig6A_GSE30881_canine_crossSpecies_BROAD_module_tumor_vs_normal.png")
out6A_pdf <- file.path(fig_dir, "Fig6A_GSE30881_canine_crossSpecies_BROAD_module_tumor_vs_normal.pdf")

ggsave(out6A_png, p6A, width = 4.8, height = 5.2, dpi = 400)
ggsave(out6A_pdf, p6A, width = 4.8, height = 5.2)

message("Saved Fig. 6A to:")
message("  ", out6A_png)
message("  ", out6A_pdf)

message("=== Phase 9 (Figures): Canine cross-species BROAD module completed successfully. ===")
