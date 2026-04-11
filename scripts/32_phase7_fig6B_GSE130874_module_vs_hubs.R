## 32_phase7_fig6B_GSE130874_module_vs_hubs.R
## هدف: ساخت Fig 6B
##  - نشان دادن اینکه در کوهورت مستقل کانین RNA-seq (GSE130874)،
##    فعالیت ماژول cross-species BROAD (score_net_z) با بیان هاب‌های شبکه
##    (TOP2A, PARP1) هم‌راستا است.
##
## ورودی‌ها:
##  - data/processed/GSE130874_vst_geneSymbol.tsv
##      (خروجی 31_phase2C_update_GSE130874_annotation_and_module_scores.R)
##  - results/tables/module_scores/module_scores_GSE130874_crossSpeciesBROAD_dog_zscore.tsv
##
## خروجی:
##  - results/figures/Fig6B_GSE130874_module_vs_TOP2A_PARP1.png
##  - results/figures/Fig6B_GSE130874_module_vs_TOP2A_PARP1.pdf

message("=== Phase 7 (Figures): Fig 6B – GSE130874 module score vs TOP2A/PARP1 ===")

## ------------------------------------------------------------
## 1. packages
## ------------------------------------------------------------
required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "tidyr", "ggplot2", "purrr"
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
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

## ------------------------------------------------------------
## 2. paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

proc_dir <- file.path(project_root, "data", "processed")
ms_dir   <- file.path(project_root, "results", "tables", "module_scores")
fig_dir  <- file.path(project_root, "results", "figures")

if (!dir.exists(proc_dir)) {
  stop("Processed data dir not found at:\n  ", proc_dir)
}
if (!dir.exists(ms_dir)) {
  stop("Module scores dir not found at:\n  ", ms_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures dir:\n  ", fig_dir)
}

expr_symbol_path <- file.path(proc_dir, "GSE130874_vst_geneSymbol.tsv")
scores_path      <- file.path(
  ms_dir,
  "module_scores_GSE130874_crossSpeciesBROAD_dog_zscore.tsv"
)

if (!file.exists(expr_symbol_path)) {
  stop("GSE130874 geneSymbol VST file not found at:\n  ", expr_symbol_path)
}
if (!file.exists(scores_path)) {
  stop("GSE130874 module score file not found at:\n  ", scores_path)
}

## ------------------------------------------------------------
## 3. load expression + scores
## ------------------------------------------------------------
message("\n--- Step 1: Load expression (geneSymbol) and module scores ---")

expr_gs <- readr::read_tsv(expr_symbol_path, show_col_types = FALSE)

if (!"gene_symbol" %in% names(expr_gs)) {
  stop("Expression file must contain a 'gene_symbol' column.")
}

sample_cols <- setdiff(names(expr_gs), "gene_symbol")

## اطمینان از numeric بودن اکسپرشن
expr_gs <- expr_gs %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

scores_130874 <- readr::read_tsv(scores_path, show_col_types = FALSE)

needed_score_cols <- c("sample_id", "score_net_z")
missing_sc <- setdiff(needed_score_cols, names(scores_130874))
if (length(missing_sc) > 0) {
  stop(
    "Module score file is missing columns: ",
    paste(missing_sc, collapse = ", ")
  )
}

message("  Expression genes: ", nrow(expr_gs))
message("  Samples (expression): ", length(sample_cols))
message("  Samples (scores)    : ", nrow(scores_130874))

## ------------------------------------------------------------
## 4. extract TOP2A / PARP1 expression و z-score
## ------------------------------------------------------------
message("\n--- Step 2: Extract TOP2A and PARP1 expression ---")

hubs_of_interest <- c("TOP2A", "PARP1")

expr_hubs <- expr_gs %>%
  dplyr::filter(gene_symbol %in% hubs_of_interest)

if (nrow(expr_hubs) == 0) {
  stop(
    "Neither TOP2A nor PARP1 found in GSE130874 expression matrix after Dog10K annotation.\n",
    "این غیرعادی است؛ لطفاً فایل GSE130874_vst_geneSymbol.tsv را چک کن."
  )
}

## ممکن است یکی از دو ژن نباشد؛ همان‌هایی که موجودند را نگه می‌داریم
found_hubs <- unique(expr_hubs$gene_symbol)
message("  Hubs found in expression matrix: ",
        paste(found_hubs, collapse = ", "))

expr_hubs_long <- expr_hubs %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(sample_cols),
    names_to  = "sample_id",
    values_to = "expr_raw"
  ) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::mutate(
    expr_z = as.numeric(scale(expr_raw))
  ) %>%
  dplyr::ungroup()

## ------------------------------------------------------------
## 5. merge با module scores
## ------------------------------------------------------------
message("\n--- Step 3: Merge expression with module scores ---")

plot_df <- expr_hubs_long %>%
  dplyr::left_join(scores_130874, by = "sample_id")

## حذف نمونه‌هایی که برایشان score_net_z یا expr_z NA است
n_before <- nrow(plot_df)
plot_df <- plot_df %>%
  dplyr::filter(!is.na(score_net_z), !is.na(expr_z))
n_after <- nrow(plot_df)

message("  Rows in combined table (gene × sample): ", n_after,
        " (from ", n_before, " before filtering NAs)")

if (n_after == 0) {
  stop(
    "No non-NA combination of expression and module scores found.\n",
    "Check if score_net_z and expr values are correctly computed."
  )
}

## ------------------------------------------------------------
## 6. محاسبه r و p برای هر ژن
## ------------------------------------------------------------
message("\n--- Step 4: Compute correlation (expr_z vs score_net_z) for each hub ---")

cor_stats <- plot_df %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    r      = suppressWarnings(
      cor(expr_z, score_net_z, use = "complete.obs", method = "pearson")
    ),
    p      = tryCatch(
      {
        stats::cor.test(expr_z, score_net_z)$p.value
      },
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    r_label = sprintf("r = %.2f", r),
    p_label = dplyr::case_when(
      is.na(p) ~ "p = NA",
      p < 0.001 ~ "p < 0.001",
      TRUE      ~ sprintf("p = %.3f", p)
    ),
    annot   = paste0(r_label, ", ", p_label)
  )

message("Correlation summary:")
print(cor_stats)

## این جدول را برای log و روش می‌نویسیم
cor_out_path <- file.path(
  project_root,
  "results", "tables", "module_scores",
  "GSE130874_module_vs_hubs_correlation_TOP2A_PARP1.tsv"
)
readr::write_tsv(cor_stats, cor_out_path)
message("Saved correlation summary to:")
message("  ", cor_out_path)

## ------------------------------------------------------------
## 7. ساخت شکل (scatter + regression line, facet by gene_symbol)
## ------------------------------------------------------------
message("\n--- Step 5: Build Fig 6B scatter plot ---")

## برای annotate کردن روی هر facette
plot_df_annot <- plot_df %>%
  dplyr::left_join(
    cor_stats %>% dplyr::select(gene_symbol, annot),
    by = "gene_symbol"
  )

p <- ggplot(plot_df_annot,
            aes(x = expr_z, y = score_net_z)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3, color = "grey70") +
  geom_point(
    size  = 2.6,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "lm",
    se     = FALSE,
    size   = 0.8
  ) +
  facet_wrap(~ gene_symbol, nrow = 1, scales = "free_x") +
  labs(
    x = "Hub expression (z-score)",
    y = "Cross-species BROAD module score (z)",
    title    = "GSE130874: module activity vs hub expression",
    subtitle = "Canine DLBCL RNA-seq cohort (Dog10K annotation)"
  ) +
  geom_text(
    data = cor_stats,
    aes(
      x     = -Inf,
      y     = Inf,
      label = annot
    ),
    hjust = -0.05,
    vjust = 1.2,
    size  = 3.2,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

out_png <- file.path(
  fig_dir,
  "Fig6B_GSE130874_module_vs_TOP2A_PARP1.png"
)
out_pdf <- file.path(
  fig_dir,
  "Fig6B_GSE130874_module_vs_TOP2A_PARP1.pdf"
)

ggsave(out_png, p, width = 8.0, height = 4.0, dpi = 400)
ggsave(out_pdf, p, width = 8.0, height = 4.0)

message("Saved Fig 6B to:")
message("  ", out_png)
message("  ", out_pdf)

message("=== Fig 6B generation completed successfully. ===")
