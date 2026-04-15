
## 61_phase9F_MKI67_correlation_module_score.R
## Correlate module score with MKI67 expression in human cohorts

message("=== Phase 9F: MKI67 correlation with module score ===")

required_pkgs <- c("readr", "dplyr", "stringr", "tibble", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      "))"
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

project_root <- normalizePath(
  "D:/Research/My Articles/DLBCL drug",
  winslash = "/",
  mustWork = TRUE
)

data_dir   <- file.path(project_root, "data", "processed")
score_dir  <- file.path(project_root, "results", "tables", "module_scores")
table_dir  <- file.path(project_root, "results", "tables", "survival")
fig_dir    <- file.path(project_root, "results", "figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## --------------------------------------------------------------------
## helpers
## --------------------------------------------------------------------
get_gene_vector <- function(expr_path, gene_symbol = "MKI67") {
  expr <- readr::read_tsv(expr_path, show_col_types = FALSE)
  
  gene_col <- c("gene_symbol", "Gene symbol", "Gene_symbol", "GeneSymbol", "symbol")
  gene_col <- gene_col[gene_col %in% names(expr)][1]
  
  if (is.na(gene_col)) {
    stop("No gene symbol column found in expression file:\n  ", expr_path)
  }
  
  expr2 <- expr %>%
    mutate(gene_use = as.character(.data[[gene_col]])) %>%
    filter(gene_use == gene_symbol)
  
  if (nrow(expr2) == 0) {
    stop("Gene ", gene_symbol, " not found in:\n  ", expr_path)
  }
  
  if (nrow(expr2) > 1) {
    expr2 <- expr2 %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE))
  } else {
    expr2 <- expr2 %>%
      select(where(is.numeric))
  }
  
  vals <- as.numeric(expr2[1, ])
  names(vals) <- colnames(expr2)
  vals
}

run_corr <- function(score_df, score_col, sample_col, mki67_vec, cohort_label) {
  df <- score_df %>%
    mutate(sample_use = as.character(.data[[sample_col]]),
           score_use  = as.numeric(.data[[score_col]])) %>%
    mutate(MKI67_expr = unname(mki67_vec[sample_use])) %>%
    filter(!is.na(score_use), !is.na(MKI67_expr))
  
  if (nrow(df) < 10) {
    stop("Too few matched samples for ", cohort_label)
  }
  
  ct <- cor.test(df$score_use, df$MKI67_expr, method = "pearson")
  
  out <- tibble(
    cohort = cohort_label,
    n = nrow(df),
    score_column = score_col,
    correlation = unname(ct$estimate),
    p_value = ct$p.value
  )
  
  list(df = df, stats = out)
}

## --------------------------------------------------------------------
## GSE56315
## --------------------------------------------------------------------
score_56315_path <- file.path(
  score_dir,
  "module_scores_GSE56315_Tier2_GSE56315_signature_zscore_with_group.tsv"
)
expr_56315_path <- file.path(data_dir, "GSE56315_expr_log2_qcfiltered.tsv")

score_56315 <- readr::read_tsv(score_56315_path, show_col_types = FALSE)
mki67_56315 <- get_gene_vector(expr_56315_path, "MKI67")

res_56315 <- run_corr(
  score_df = score_56315,
  score_col = c("score_net_z", "score_net")[c("score_net_z", "score_net") %in% names(score_56315)][1],
  sample_col = "sample_id",
  mki67_vec = mki67_56315,
  cohort_label = "GSE56315"
)

## --------------------------------------------------------------------
## GSE31312
## --------------------------------------------------------------------
score_31312_path <- file.path(table_dir, "GSE31312_OS_scores_and_clinical_from_pdf.tsv")
expr_31312_path  <- file.path(data_dir, "GSE31312_expr_log2_qcfiltered.tsv")

score_31312 <- readr::read_tsv(score_31312_path, show_col_types = FALSE)
mki67_31312 <- get_gene_vector(expr_31312_path, "MKI67")

sample_col_31312 <- c("sample_id", "case_geo_id")
sample_col_31312 <- sample_col_31312[sample_col_31312 %in% names(score_31312)][1]

res_31312 <- run_corr(
  score_df = score_31312,
  score_col = c("score_net_z", "score_net")[c("score_net_z", "score_net") %in% names(score_31312)][1],
  sample_col = sample_col_31312,
  mki67_vec = mki67_31312,
  cohort_label = "GSE31312"
)

## --------------------------------------------------------------------
## save stats
## --------------------------------------------------------------------
stats_all <- bind_rows(res_56315$stats, res_31312$stats)

out_stats <- file.path(table_dir, "MKI67_module_score_correlations.tsv")
readr::write_tsv(stats_all, out_stats)

message("Correlation summary:")
print(stats_all)

## --------------------------------------------------------------------
## plots
## --------------------------------------------------------------------
make_plot <- function(df, cohort_label, cor_val, p_val) {
  ggplot(df, aes(x = MKI67_expr, y = score_use)) +
    geom_point(size = 1.8, alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    labs(
      x = "MKI67 expression",
      y = "Module score",
      title = cohort_label
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("r = ", sprintf("%.2f", cor_val), "\np = ", format(p_val, digits = 3, scientific = TRUE)),
      hjust = 1.05, vjust = 1.3, size = 4
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

p1 <- make_plot(
  res_56315$df,
  "GSE56315",
  res_56315$stats$correlation,
  res_56315$stats$p_value
)

p2 <- make_plot(
  res_31312$df,
  "GSE31312",
  res_31312$stats$correlation,
  res_31312$stats$p_value
)

ggsave(file.path(fig_dir, "SFig_MKI67_vs_module_GSE56315.png"), p1, width = 5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "SFig_MKI67_vs_module_GSE31312.png"), p2, width = 5, height = 4.5, dpi = 300)

message("Saved to:")
message("  ", out_stats)
message("  ", file.path(fig_dir, "SFig_MKI67_vs_module_GSE56315.png"))
message("  ", file.path(fig_dir, "SFig_MKI67_vs_module_GSE31312.png"))
message("=== Done ===")