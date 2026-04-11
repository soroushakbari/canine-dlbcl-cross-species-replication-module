## 35_phase8B_CMap_drug_level_robustness.R
## Drug-level robustness analysis of CMap hits using empirical tail + random drug sets
## Uses existing per-drug summary table and core network-hit drugs.

message("=== Phase 8B: CMap drug-level robustness (empirical tails + random drug sets) ===")

## ------------------------------------------------------------
## 1. packages
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
## 2. paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

drug_dir <- file.path(project_root, "results", "tables", "Drug")
fig_dir  <- file.path(project_root, "results", "figures")

if (!dir.exists(drug_dir)) {
  stop("Drug results directory not found at:\n  ", drug_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

rank_path <- file.path(
  drug_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_ranked.tsv"
)

core_path <- file.path(
  drug_dir,
  "CMap_core_network_hits_for_manuscript.tsv"
)

if (!file.exists(rank_path)) {
  stop("Per-drug summary table not found at:\n  ", rank_path,
       "\nاین فایل در Phase 7 ساخته می‌شود (20_phase7_...).")
}
if (!file.exists(core_path)) {
  stop("Core network-hit drug table not found at:\n  ", core_path,
       "\nاین فایل در Phase 6B ساخته می‌شود (24_phase6B_...).")
}

## ------------------------------------------------------------
## 3. load per-drug summary
## ------------------------------------------------------------
message("\n--- Step 1: Load per-drug summary table ---")
drug_all <- readr::read_tsv(rank_path, show_col_types = FALSE)

## حدس زدن ستون نام دارو
name_candidates <- c("drug_name", "pert_name", "pert_iname", "name")
name_col <- intersect(name_candidates, names(drug_all))
if (length(name_col) == 0) {
  stop(
    "Could not detect a drug name column in the ranked table.\n",
    "Expected one of: ", paste(name_candidates, collapse = ", ")
  )
}
name_col <- name_col[1]

## ستون score
score_candidates <- c("min_score", "min_norm_cs", "min_connectivity", "min_cs")
score_col <- intersect(score_candidates, names(drug_all))
if (length(score_col) == 0) {
  stop(
    "Could not detect a 'minimum connectivity score' column in the ranked table.\n",
    "Expected one of: ", paste(score_candidates, collapse = ", ")
  )
}
score_col <- score_col[1]

message("  Detected drug name column: ", name_col)
message("  Detected min-score column: ", score_col)

drug_all2 <- drug_all %>%
  dplyr::select(
    drug_name = .data[[name_col]],
    min_score = .data[[score_col]]
  ) %>%
  dplyr::mutate(
    drug_name = stringr::str_squish(as.character(drug_name)),
    min_score = as.numeric(min_score)
  ) %>%
  dplyr::filter(!is.na(drug_name), !is.na(min_score))

n_all <- nrow(drug_all2)
message("  Number of small-molecule drugs with valid min_score: ", n_all)

## ------------------------------------------------------------
## 4. load core network-hit drugs
## ------------------------------------------------------------
message("\n--- Step 2: Load core network-hit drugs ---")
core_tbl_raw <- readr::read_tsv(core_path, show_col_types = FALSE)

core_name_candidates <- c("drug_name", "pert_name", "pert_iname", "name")
core_name_col <- intersect(core_name_candidates, names(core_tbl_raw))
if (length(core_name_col) == 0) {
  stop(
    "Could not detect a drug name column in the core-drug table.\n",
    "Expected one of: ", paste(core_name_candidates, collapse = ", ")
  )
}
core_name_col <- core_name_col[1]

core_tbl <- core_tbl_raw %>%
  dplyr::mutate(
    drug_name = stringr::str_squish(as.character(.data[[core_name_col]]))
  ) %>%
  dplyr::select(dplyr::any_of(c("drug_name", "MOA_category", "moa"))) %>%
  dplyr::distinct()

message("  Core network-hit drugs (from file):")
print(core_tbl)

## ------------------------------------------------------------
## 5. merge: add min_score + empirical rank / percentile
## ------------------------------------------------------------
message("\n--- Step 3: Empirical ranks and percentiles ---")

## تمام داروها: rank و percentile
drug_all3 <- drug_all2 %>%
  dplyr::arrange(min_score) %>%  # بیشتر منفی → بالاتر
  dplyr::mutate(
    rank_min       = dplyr::row_number(),          # 1 = strongest negative
    percentile_min = rank_min / n()
  )

## join با core
core_merged <- core_tbl %>%
  dplyr::left_join(drug_all3, by = "drug_name")

if (any(is.na(core_merged$min_score))) {
  warning("Some core drugs were not found in the ranked table. Check name matching.")
}

message("  Core drugs with empirical position among all drugs:")
print(core_merged)

core_out_path <- file.path(
  drug_dir,
  "CMap_core_network_hits_with_empirical_tail.tsv"
)
readr::write_tsv(core_merged, core_out_path)
message("  Saved core drug empirical table to:\n  ", core_out_path)

## ------------------------------------------------------------
## 6. random-set robustness: best min_score in random k-drug sets
## ------------------------------------------------------------
message("\n--- Step 4: Random k-drug set analysis ---")

## فقط داروهایی را نگه داریم که min_score غیر NA دارند
scores_all <- drug_all3$min_score
names(scores_all) <- drug_all3$drug_name

## core drugs که min_score دارند
core_valid <- core_merged %>%
  dplyr::filter(!is.na(min_score))

k_core <- nrow(core_valid)
if (k_core < 2) {
  stop("Fewer than 2 core drugs with valid min_score; random-set analysis is not meaningful.")
}

message("  Number of core drugs with valid scores: ", k_core)

## observed best (strongest negative) min_score در بین core
obs_best <- min(core_valid$min_score, na.rm = TRUE)
best_drug <- core_valid$drug_name[which.min(core_valid$min_score)]
message(sprintf(
  "  Observed best min_score among core drugs: %.3f (drug: %s)",
  obs_best, best_drug
))

## random sets
set.seed(12345)
n_iter <- 10000L

message("  Running ", n_iter, " random k-drug sets (this may take a few seconds) ...")

best_random <- numeric(n_iter)
n_pool <- length(scores_all)

for (i in seq_len(n_iter)) {
  idx <- sample.int(n_pool, size = k_core, replace = FALSE)
  best_random[i] <- min(scores_all[idx])
}

## empirical p: احتمال این‌که یک set تصادفی k-تایی
## min_score ≤ min_score core داشته باشد (چون هرچه منفی‌تر بهتر است)
emp_p <- mean(best_random <= obs_best)

message(sprintf(
  "  Empirical p (random k-drug set has min_score <= observed): %.4f",
  emp_p
))

## خلاصه و ذخیره
robust_summary <- tibble::tibble(
  n_all_drugs         = n_all,
  n_core_drugs        = k_core,
  observed_best_min   = obs_best,
  observed_best_drug  = best_drug,
  n_iterations        = n_iter,
  empirical_p_best    = emp_p
)

robust_summary_path <- file.path(
  drug_dir,
  "CMap_random_kset_bestMinScore_summary.tsv"
)
readr::write_tsv(robust_summary, robust_summary_path)
message("  Saved robustness summary to:\n  ", robust_summary_path)

## همچنین کل توزیع best_random را (برای histogram) ذخیره کنیم
robust_dist_path <- file.path(
  drug_dir,
  "CMap_random_kset_bestMinScore_distribution.tsv"
)
readr::write_tsv(
  tibble::tibble(best_random_min_score = best_random),
  robust_dist_path
)
message("  Saved random-set distribution to:\n  ", robust_dist_path)

## ------------------------------------------------------------
## 7. plot: histogram of random best min_scores + observed line
## ------------------------------------------------------------
message("\n--- Step 5: Build supplementary robustness figure ---")

plot_df <- tibble::tibble(best_random_min_score = best_random)

p_hist <- ggplot(plot_df, aes(x = best_random_min_score)) +
  geom_histogram(
    bins = 40,
    fill = "#deebf7",
    color = "#08519c",
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = obs_best,
    color = "#cb181d",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = obs_best,
    y = Inf,
    label = sprintf("Observed core\nbest = %.2f", obs_best),
    vjust = 1.6,
    hjust = 0,
    size = 3.3,
    color = "#cb181d"
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "Best (most negative) min_score in random k-drug set",
    y = "Count",
    title = "Robustness of core CMap hits vs random drug sets",
    subtitle = paste0(
      "k = ", k_core, " core drugs; empirical p = ",
      sprintf("%.3g", emp_p)
    )
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10)
  )

fig_png <- file.path(
  fig_dir,
  "SFig_CMap_random_kset_bestMinScore.png"
)
fig_pdf <- file.path(
  fig_dir,
  "SFig_CMap_random_kset_bestMinScore.pdf"
)

ggsave(fig_png, p_hist, width = 6.5, height = 4.8, dpi = 400)
ggsave(fig_pdf, p_hist, width = 6.5, height = 4.8)

message("  Saved robustness figure to:\n  ", fig_png, "\n  ", fig_pdf)

message("\n=== Phase 8B completed successfully. ===")
