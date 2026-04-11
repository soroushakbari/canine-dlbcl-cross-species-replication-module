## 40_phase9C_GSE31312_OS_from_pdf.R
## Use manual clinical sheet (from PDF) to do OS analysis for GSE31312

message("=== Phase 9C (alt): GSE31312 OS analysis using PDF-derived clinical sheet ===")

## ------------------------------------------------------------
## 1. packages
## ------------------------------------------------------------
required_pkgs <- c(
  "readr", "dplyr", "tibble", "tidyr", "stringr",
  "survival", "survminer"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('"%s"', required_pkgs), collapse = ", "),
      ")) و دوباره اسکریپت را اجرا کن."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(survival)
  library(survminer)
})

## ------------------------------------------------------------
## 2. paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

meta_dir <- file.path(project_root, "metadata")
ms_dir   <- file.path(project_root, "results", "tables", "module_scores")
fig_dir  <- file.path(project_root, "results", "figures")
surv_dir <- file.path(project_root, "results", "tables", "survival")

if (!dir.exists(fig_dir))  dir.create(fig_dir,  recursive = TRUE)
if (!dir.exists(surv_dir)) dir.create(surv_dir, recursive = TRUE)

meta_path_geo   <- file.path(meta_dir, "GSE31312_sample_metadata.csv")
scores_path     <- file.path(
  ms_dir,
  "module_scores_GSE31312_Tier2_GSE56315_signature_zscore.tsv"
)
clin_pdf_path   <- file.path(meta_dir, "GSE31312_clinical_from_pdf.csv")

if (!file.exists(meta_path_geo)) {
  stop("Metadata file GSE31312_sample_metadata.csv not found at:\n  ", meta_path_geo)
}
if (!file.exists(scores_path)) {
  stop("Module scores file not found at:\n  ", scores_path)
}
if (!file.exists(clin_pdf_path)) {
  stop("Clinical PDF-derived CSV not found at:\n  ", clin_pdf_path,
       "\nلطفاً اکسل را به صورت CSV با همین نام ذخیره کن.")
}

## ------------------------------------------------------------
## 3. load module scores & GEO metadata & clinical PDF sheet
## ------------------------------------------------------------
message("\n--- Step 1: Load module scores ---")
scores_df <- readr::read_tsv(scores_path, show_col_types = FALSE)

if (!all(c("sample_id", "score_net_z") %in% names(scores_df))) {
  stop("Scores table must contain 'sample_id' و 'score_net_z'.")
}
message("  Samples with module scores: ", nrow(scores_df))

message("\n--- Step 2: Load GEO sample metadata ---")
meta_geo <- readr::read_csv(meta_path_geo, show_col_types = FALSE)
if (!"sample_id" %in% names(meta_geo)) {
  stop("GEO metadata must contain a 'sample_id' column.")
}
message("  Rows in GEO metadata: ", nrow(meta_geo))

message("\n--- Step 3: Load PDF-derived clinical sheet ---")
clin_pdf <- readr::read_csv(clin_pdf_path, show_col_types = FALSE)
message("  Rows in clinical sheet: ", nrow(clin_pdf))

## ------------------------------------------------------------
## 4. clean clinical sheet: extract N0001-style ID, OS time & event
## ------------------------------------------------------------
message("\n--- Step 4: Clean clinical sheet and extract OS info ---")

if (!"GEO Depository #" %in% names(clin_pdf)) {
  stop("Clinical CSV must contain a column named 'GEO Depository #'.")
}

clin_clean <- clin_pdf %>%
  mutate(
    geo_dep_raw = `GEO Depository #`,
    ## extract N0001-style ID
    case_geo_id = stringr::str_extract(as.character(geo_dep_raw), "N[0-9]+"),
    ## assume PFS and OS are months (numeric-like)
    PFS_months = suppressWarnings(as.numeric(PFS)),
    OS_months  = suppressWarnings(as.numeric(OS)),
    PFScensor_num = suppressWarnings(as.numeric(PFScensor)),
    OScensor_num  = suppressWarnings(as.numeric(OScensor)),
    ## event: 1 = event, 0 = censored
    pfs_event = dplyr::case_when(
      PFScensor_num == 0 ~ 1,
      PFScensor_num == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    os_event = dplyr::case_when(
      OScensor_num == 0 ~ 1,
      OScensor_num == 1 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(case_geo_id))

message("  Clinical rows with valid case_geo_id: ", nrow(clin_clean))

## ------------------------------------------------------------
## 5. map case_geo_id (N0001) to GSMs using GEO metadata
## ------------------------------------------------------------
message("\n--- Step 5: Map case_geo_id to GEO sample_id (GSM) ---")

## اول همه‌ی ستون‌های غیر از sample_id را به character تبدیل می‌کنیم
meta_geo_for_map <- meta_geo %>%
  dplyr::mutate(
    dplyr::across(
      -sample_id,
      ~ as.character(.)
    )
  )

meta_long <- meta_geo_for_map %>%
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "meta_field",
    values_to = "meta_value"
  ) %>%
  dplyr::filter(!is.na(meta_value), meta_value != "") %>%
  dplyr::mutate(
    case_geo_id = stringr::str_extract(meta_value, "N[0-9]+")
  ) %>%
  dplyr::filter(!is.na(case_geo_id)) %>%
  dplyr::distinct(sample_id, case_geo_id)

message("  Sample–case_geo_id map rows: ", nrow(meta_long))

clin_mapped <- clin_clean %>%
  inner_join(meta_long, by = "case_geo_id")

message("  Clinical rows mapped to GSM sample_id: ", nrow(clin_mapped))

if (nrow(clin_mapped) == 0L) {
  stop("No rows in clinical sheet could be mapped to GEO samples via case_geo_id.")
}

clin_mapped_small <- clin_mapped %>%
  select(
    sample_id,
    case_geo_id,
    PFS_months,
    OS_months,
    pfs_event,
    os_event
  )

## ------------------------------------------------------------
## 6. merge module scores with clinical info
## ------------------------------------------------------------
message("\n--- Step 6: Merge module scores with clinical OS data ---")

surv_df <- scores_df %>%
  inner_join(clin_mapped_small, by = "sample_id")

surv_os <- surv_df %>%
  filter(!is.na(OS_months), !is.na(os_event))

n_total  <- nrow(surv_os)
n_events <- sum(surv_os$os_event == 1, na.rm = TRUE)

message("  Samples with complete OS + score: ", n_total)
message("  Number of OS events: ", n_events)

if (n_total < 50L || n_events < 20L) {
  warning("OS dataset is relatively small/low events (n = ", n_total,
          ", events = ", n_events, "); power may be limited.")
}

out_surv_table <- file.path(
  surv_dir,
  "GSE31312_OS_scores_and_clinical_from_pdf.tsv"
)
readr::write_tsv(surv_os, out_surv_table)
message("  Saved merged OS table to:\n  ", out_surv_table)

if (n_total == 0L) {
  stop("No samples with complete OS data after merge; cannot run Cox/KM.")
}

## ------------------------------------------------------------
## 7. Cox model (continuous module score)
## ------------------------------------------------------------
message("\n--- Step 7: Cox model with continuous module score ---")

cox_fit <- survival::coxph(
  Surv(OS_months, os_event) ~ score_net_z,
  data = surv_os
)

cox_sum <- summary(cox_fit)

hr  <- cox_sum$coef["score_net_z", "exp(coef)"]
lcl <- cox_sum$conf.int["score_net_z", "lower .95"]
ucl <- cox_sum$conf.int["score_net_z", "upper .95"]
p   <- cox_sum$coef["score_net_z", "Pr(>|z|)"]

cox_out <- tibble::tibble(
  term          = "score_net_z",
  HR            = hr,
  CI_lower_95   = lcl,
  CI_upper_95   = ucl,
  p_value       = p,
  n_total       = n_total,
  n_events      = n_events
)

out_cox <- file.path(
  surv_dir,
  "GSE31312_OS_cox_Tier2_score_net_from_pdf.tsv"
)
readr::write_tsv(cox_out, out_cox)

message(
  sprintf(
    "  Cox HR per 1 SD increase in score_net_z: %.2f (95%% CI %.2f–%.2f), p = %.3g",
    hr, lcl, ucl, p
  )
)
message("  Saved Cox summary to:\n  ", out_cox)

## ------------------------------------------------------------
## 8. KM high vs low (median split)
## ------------------------------------------------------------
message("\n--- Step 8: Kaplan–Meier (high vs low module score) ---")

median_score <- stats::median(surv_os$score_net_z, na.rm = TRUE)

surv_os <- surv_os %>%
  mutate(
    score_group = if_else(score_net_z >= median_score, "High", "Low")
  )

km_fit <- survival::survfit(
  Surv(OS_months, os_event) ~ score_group,
  data = surv_os
)

survdiff_res <- survival::survdiff(
  Surv(OS_months, os_event) ~ score_group,
  data = surv_os
)

chisq <- survdiff_res$chisq
p_lr  <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)

message("  Log-rank p (High vs Low): ", signif(p_lr, 3))

km_out_table <- file.path(
  surv_dir,
  "GSE31312_OS_KM_high_vs_low_module_scores_from_pdf.tsv"
)
readr::write_tsv(surv_os, km_out_table)
message("  Saved KM sample table to:\n  ", km_out_table)

## ------------------------------------------------------------
## 9. KM figure
## ------------------------------------------------------------
message("\n--- Step 9: Build KM figure ---")

km_plot <- survminer::ggsurvplot(
  km_fit,
  data           = surv_os,
  risk.table     = TRUE,
  pval           = TRUE,
  pval.method    = TRUE,
  conf.int       = FALSE,
  legend.title   = "Module score",
  legend.labs    = c("High", "Low"),
  xlab           = "Overall survival (months)",
  ylab           = "Survival probability",
  ggtheme        = ggplot2::theme_bw(base_size = 12)
)

out_png <- file.path(
  fig_dir,
  "Fig6D_GSE31312_OS_KM_high_vs_low_module_from_pdf.png"
)
out_pdf <- file.path(
  fig_dir,
  "Fig6D_GSE31312_OS_KM_high_vs_low_module_from_pdf.pdf"
)

ggplot2::ggsave(out_png, km_plot$plot, width = 6.5, height = 5, dpi = 400)
ggplot2::ggsave(out_pdf, km_plot$plot, width = 6.5, height = 5)

message("  Saved KM figure to:")
message("    ", out_png)
message("    ", out_pdf)

message("\n=== Phase 9C (OS from PDF sheet) completed successfully. ===")
