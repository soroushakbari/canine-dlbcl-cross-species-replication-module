## 36_phase8B_CMap_core_drug_percentiles.R
## Empirical rank and percentile of core CMap drugs within the full library

message("=== Phase 8B: Empirical ranks/percentiles of core CMap drugs ===")

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
## 3. load global ranked table
## ------------------------------------------------------------
message("\n--- Step 1: Load global per-drug ranked table ---")
drug_all <- readr::read_tsv(rank_path, show_col_types = FALSE)

name_candidates  <- c("drug_name", "pert_name", "pert_iname", "name")
score_candidates <- c("min_score", "min_norm_cs", "min_connectivity", "min_cs")

name_col  <- intersect(name_candidates,  names(drug_all))
score_col <- intersect(score_candidates, names(drug_all))

if (length(name_col) == 0) {
  stop(
    "Could not detect a drug name column in the ranked table.\n",
    "Expected one of: ", paste(name_candidates, collapse = ", ")
  )
}
if (length(score_col) == 0) {
  stop(
    "Could not detect a 'minimum connectivity score' column in the ranked table.\n",
    "Expected one of: ", paste(score_candidates, collapse = ", ")
  )
}

name_col  <- name_col[1]
score_col <- score_col[1]

message("  Detected drug name column : ", name_col)
message("  Detected min-score column : ", score_col)

drug_all2 <- drug_all %>%
  dplyr::transmute(
    drug_name = stringr::str_squish(as.character(.data[[name_col]])),
    min_score = as.numeric(.data[[score_col]])
  ) %>%
  dplyr::filter(!is.na(drug_name), !is.na(min_score))

n_all <- nrow(drug_all2)
message("  Number of small-molecule drugs with valid min_score: ", n_all)

drug_all3 <- drug_all2 %>%
  dplyr::arrange(min_score) %>%           # بیشتر منفی → رتبه بهتر
  dplyr::mutate(
    rank_min       = dplyr::row_number(), # 1 = strongest negative
    percentile_min = rank_min / n()       # 0–1
  )

## ------------------------------------------------------------
## 4. load core drugs and merge
## ------------------------------------------------------------
message("\n--- Step 2: Load core network-hit drugs ---")
core_raw <- readr::read_tsv(core_path, show_col_types = FALSE)

core_name_candidates <- c("drug_name", "pert_name", "pert_iname", "name")
core_name_col <- intersect(core_name_candidates, names(core_raw))
if (length(core_name_col) == 0) {
  stop(
    "Could not detect a drug name column in the core-drug table.\n",
    "Expected one of: ", paste(core_name_candidates, collapse = ", ")
  )
}
core_name_col <- core_name_col[1]

core_tbl <- core_raw %>%
  dplyr::mutate(
    drug_name = stringr::str_squish(as.character(.data[[core_name_col]]))
  ) %>%
  dplyr::select(dplyr::any_of(c("drug_name", "MOA_category", "moa"))) %>%
  dplyr::distinct()

message("  Core drugs (from file):")
print(core_tbl)

## merge
core_res <- core_tbl %>%
  dplyr::left_join(drug_all3, by = "drug_name") %>%
  dplyr::arrange(rank_min)

if (any(is.na(core_res$min_score))) {
  warning("Some core drugs were not found in the ranked table; check spelling.")
}

## ------------------------------------------------------------
## 5. save table + console summary
## ------------------------------------------------------------
out_path <- file.path(
  drug_dir,
  "CMap_core_network_hits_with_empirical_percentiles.tsv"
)
readr::write_tsv(core_res, out_path)
message("\nSaved core drug empirical rank/percentile table to:\n  ", out_path)

message("\n--- Step 3: Console summary (rank & percentile) ---")
summary_tbl <- core_res %>%
  dplyr::mutate(
    top_percent = percentile_min * 100,
    top_label   = sprintf("top %.2f%%", top_percent)
  ) %>%
  dplyr::select(drug_name, MOA_category, min_score,
                rank_min, percentile_min, top_label)

print(summary_tbl)

## ------------------------------------------------------------
## 6. plot: core drugs within global distribution (optional nice figure)
## ------------------------------------------------------------
message("\n--- Step 4: Build supplementary percentile figure ---")

p_core <- ggplot(core_res,
                 aes(x = reorder(drug_name, min_score),
                     y = min_score)) +
  geom_segment(aes(xend = drug_name, y = 0, yend = min_score),
               colour = "grey70", linewidth = 0.6) +
  geom_point(size = 3) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    x = NULL,
    y = "Minimum connectivity score (more negative = stronger reversal)",
    title = "Core CMap drugs within the global connectivity distribution"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  )

fig_png <- file.path(
  fig_dir,
  "SFig_CMap_core_drug_minScore_ranks.png"
)
fig_pdf <- file.path(
  fig_dir,
  "SFig_CMap_core_drug_minScore_ranks.pdf"
)

ggsave(fig_png, p_core, width = 6.0, height = 3.8, dpi = 400)
ggsave(fig_pdf, p_core, width = 6.0, height = 3.8)

message("Saved core-drug rank figure to:\n  ", fig_png, "\n  ", fig_pdf)

message("\n=== Phase 8B (core ranks/percentiles) completed successfully. ===")
