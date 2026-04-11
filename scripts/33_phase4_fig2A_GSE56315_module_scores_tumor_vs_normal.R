## 33_phase4_fig2A_GSE56315_module_scores_tumor_vs_normal.R
## هدف: Fig 2A
##  - نمایش توزیع نمره ماژول human DLBCL (Tier2 BROAD) در GSE56315
##    برای تومور vs نرمال و گزارش p-value و اندازه اثر.

message("=== Phase 4 (Figures): Fig 2A – GSE56315 module scores (tumor vs normal) ===")

## ------------------------------------------------------------
## 1. Packages
## ------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "ggplot2")

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

## ------------------------------------------------------------
## 2. Paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

ms_dir <- file.path(project_root, "results", "tables", "module_scores")
fig_dir <- file.path(project_root, "results", "figures")

if (!dir.exists(ms_dir)) {
  stop("Module scores directory not found at:\n  ", ms_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

ms_path <- file.path(
  ms_dir,
  "module_scores_GSE56315_Tier2_GSE56315_signature_zscore_with_group.tsv"
)

if (!file.exists(ms_path)) {
  stop("Module score file for GSE56315 not found at:\n  ", ms_path)
}

## ------------------------------------------------------------
## 3. Load table & pick columns
## ------------------------------------------------------------
message("\n--- Step 1: Load module scores table ---")

ms <- readr::read_tsv(ms_path, show_col_types = FALSE)

message("  Columns available: ", paste(names(ms), collapse = ", "))

## تشخیص ستون score: اول score_net، اگر نبود score_net_z
score_col <- NULL
if ("score_net" %in% names(ms)) {
  score_col <- "score_net"
} else if ("score_net_z" %in% names(ms)) {
  score_col <- "score_net_z"
} else {
  stop("Could not find 'score_net' or 'score_net_z' in module score file.")
}

## تشخیص ستون گروه
group_col <- NULL
if ("group_tn" %in% names(ms)) {
  group_col <- "group_tn"
} else if ("group" %in% names(ms)) {
  group_col <- "group"
} else {
  stop("Could not find a group column (expected 'group_tn' or 'group').")
}

df <- ms %>%
  dplyr::select(sample_id = dplyr::any_of("sample_id"),
                !!group_col, !!score_col) %>%
  dplyr::rename(
    group_raw = !!group_col,
    score     = !!score_col
  ) %>%
  dplyr::mutate(
    group_clean = stringr::str_to_lower(stringr::str_trim(group_raw))
  ) %>%
  dplyr::filter(group_clean %in% c("tumor", "normal"))

if (nrow(df) == 0) {
  stop("After filtering, no rows with group = 'tumor' or 'normal' remained.\n",
       "Check the group coding in the module score file.")
}

## ترتیب گروه‌ها: نرمال چپ، تومور راست
df <- df %>%
  dplyr::mutate(
    group = factor(
      dplyr::case_when(
        group_clean == "normal" ~ "Normal",
        group_clean == "tumor"  ~ "DLBCL",
        TRUE                    ~ NA_character_
      ),
      levels = c("Normal", "DLBCL")
    )
  ) %>%
  dplyr::filter(!is.na(group))

message("  Sample counts by group:")
print(table(df$group))

## ------------------------------------------------------------
## 4. T-test و اندازه اثر
## ------------------------------------------------------------
message("\n--- Step 2: T-test and effect size ---")

normal_scores <- df %>% dplyr::filter(group == "Normal") %>% dplyr::pull(score)
tumor_scores  <- df %>% dplyr::filter(group == "DLBCL") %>% dplyr::pull(score)

tt <- stats::t.test(tumor_scores, normal_scores, var.equal = FALSE)

p_val <- tt$p.value
mean_normal <- mean(normal_scores)
mean_tumor  <- mean(tumor_scores)
sd_normal   <- sd(normal_scores)
sd_tumor    <- sd(tumor_scores)

## pooled SD برای Cohen's d
n_n <- length(normal_scores)
n_t <- length(tumor_scores)
sd_pooled <- sqrt(((n_n - 1) * sd_normal^2 + (n_t - 1) * sd_tumor^2) /
                    (n_n + n_t - 2))
cohen_d <- (mean_tumor - mean_normal) / sd_pooled

p_label <- dplyr::case_when(
  is.na(p_val)         ~ "p = NA",
  p_val < 1e-20        ~ "p < 1e-20",
  p_val < 0.001        ~ sprintf("p = %.1e", p_val),
  TRUE                 ~ sprintf("p = %.3f", p_val)
)

summary_tbl <- tibble::tibble(
  group       = c("Normal", "DLBCL"),
  n           = c(n_n, n_t),
  mean_score  = c(mean_normal, mean_tumor),
  sd_score    = c(sd_normal, sd_tumor)
) %>%
  tidyr::pivot_wider(
    names_from  = group,
    values_from = c(n, mean_score, sd_score)
  ) %>%
  dplyr::mutate(
    p_value = p_val,
    cohen_d = cohen_d
  )

summary_out <- file.path(
  ms_dir,
  "GSE56315_module_scores_tumor_vs_normal_summary.tsv"
)
readr::write_tsv(summary_tbl, summary_out)

message("  T-test p-value: ", signif(p_val, 3), " (", p_label, ")")
message("  Cohen's d (DLBCL vs Normal): ", round(cohen_d, 2))
message("Saved summary stats to:")
message("  ", summary_out)

## ------------------------------------------------------------
## 5. شکل Fig 2A
## ------------------------------------------------------------
message("\n--- Step 3: Build Fig 2A violin/box/jitter plot ---")

y_range <- range(df$score, na.rm = TRUE)
y_pad   <- diff(y_range) * 0.15
annot_y <- max(y_range) + y_pad * 0.4

p <- ggplot(df, aes(x = group, y = score)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.6
  ) +
  geom_boxplot(
    width         = 0.15,
    outlier.shape = NA,
    alpha         = 0.8
  ) +
  geom_jitter(
    width = 0.08,
    size  = 1.7,
    alpha = 0.75
  ) +
  geom_segment(
    data = NULL,
    aes(
      x    = 1,
      xend = 2,
      y    = annot_y,
      yend = annot_y
    ),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = NULL,
    aes(
      x     = 1.5,
      y     = annot_y,
      label = p_label
    ),
    vjust = -0.4,
    size  = 3.5,
    inherit.aes = FALSE
  ) +
  labs(
    x = NULL,
    y = "Cross-species BROAD module score (GSE56315 Tier2)",
    title    = "Human DLBCL module activity in GSE56315",
    subtitle = "Module scores derived from the GSE56315 Tier2 BROAD signature (tumour vs normal tonsil)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10),
    axis.text.x     = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  )

out_png <- file.path(fig_dir, "Fig2A_GSE56315_module_scores_tumour_vs_normal.png")
out_pdf <- file.path(fig_dir, "Fig2A_GSE56315_module_scores_tumour_vs_normal.pdf")

ggsave(out_png, p, width = 4.8, height = 4.6, dpi = 400)
ggsave(out_pdf, p, width = 4.8, height = 4.6)

message("Saved Fig 2A to:")
message("  ", out_png)
message("  ", out_pdf)

message("=== Fig 2A generation completed successfully. ===")
