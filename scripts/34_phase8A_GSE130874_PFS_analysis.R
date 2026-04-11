## 34_phase8A_GSE130874_PFS_analysis.R
## PFS / response analysis for canine DLBCL RNA-seq cohort (GSE130874)
## using cross-species BROAD module scores

message("=== Phase 8A: GSE130874 PFS / response analysis (cross-species BROAD) ===")

## ------------------------------------------------------------
## 1. packages
## ------------------------------------------------------------
required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "survival", "survminer", "ggplot2"
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
  library(survival)
  library(survminer)
  library(ggplot2)
})

## ------------------------------------------------------------
## 2. paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

meta_dir  <- file.path(project_root, "metadata")
score_dir <- file.path(project_root, "results", "tables", "module_scores")
surv_dir  <- file.path(project_root, "results", "tables", "survival")
fig_dir   <- file.path(project_root, "results", "figures")

if (!dir.exists(surv_dir)) {
  dir.create(surv_dir, recursive = TRUE)
  message("Created survival directory:\n  ", surv_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

scores_path <- file.path(
  score_dir,
  "module_scores_GSE130874_crossSpeciesBROAD_dog_zscore.tsv"
)

clin_path <- file.path(
  meta_dir,
  "GSE130874_clinical_PFS_response.csv"
)

if (!file.exists(scores_path)) {
  stop("Module scores file not found at:\n  ", scores_path)
}
if (!file.exists(clin_path)) {
  stop(
    "Clinical metadata file for GSE130874 not found at:\n  ", clin_path, "\n",
    "Please create 'GSE130874_clinical_PFS_response.csv' in the metadata folder\n",
    "with at least the columns: sample_id, PFS_days, PFS_status (0/1),\n",
    "and optionally response_group."
  )
}

## ------------------------------------------------------------
## 3. load module scores
## ------------------------------------------------------------
message("\n--- Step 1: Load module scores ---")

scores <- readr::read_tsv(scores_path, show_col_types = FALSE)

## تلاش برای پیدا کردن ستون نمونه
sample_col_candidates <- c("sample_id", "sample", "Sample", "gsm")
sample_col <- intersect(sample_col_candidates, names(scores))

if (length(sample_col) == 0) {
  stop(
    "Could not detect a sample ID column in the module scores file.\n",
    "Expected one of: ", paste(sample_col_candidates, collapse = ", ")
  )
}
sample_col <- sample_col[1]

if (!"score_net_z" %in% names(scores)) {
  stop("Module scores table must contain 'score_net_z' column.")
}

scores2 <- scores %>%
  dplyr::select(
    sample_id = .data[[sample_col]],
    score_net_z
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id) %>% stringr::str_squish()
  )

message("  Samples with module scores: ", nrow(scores2))

## ------------------------------------------------------------
## 4. load clinical metadata
## ------------------------------------------------------------
message("\n--- Step 2: Load clinical PFS / response metadata ---")

clin <- readr::read_csv(clin_path, show_col_types = FALSE)

needed_cols <- c("sample_id", "PFS_days", "PFS_status")
missing <- setdiff(needed_cols, names(clin))
if (length(missing) > 0) {
  stop(
    "Clinical file is missing required columns: ",
    paste(missing, collapse = ", "), "\n",
    "Columns present: ", paste(names(clin), collapse = ", ")
  )
}

clin2 <- clin %>%
  dplyr::mutate(
    sample_id  = as.character(sample_id) %>% stringr::str_squish(),
    PFS_days   = as.numeric(PFS_days),
    PFS_status = as.integer(PFS_status)
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    !is.na(PFS_days),
    !is.na(PFS_status)
  )

message("  Clinical records with non-missing PFS: ", nrow(clin2))

## optional response_group
has_response <- "response_group" %in% names(clin2)

if (has_response) {
  clin2 <- clin2 %>%
    dplyr::mutate(
      response_group = as.factor(response_group)
    )
  message("  Distinct response groups: ",
          paste(levels(clin2$response_group), collapse = ", "))
} else {
  message("  No 'response_group' column detected; response analysis will be skipped.")
}

## ------------------------------------------------------------
## 5. merge scores + clinical
## ------------------------------------------------------------
message("\n--- Step 3: Merge scores with clinical data ---")

dat <- scores2 %>%
  dplyr::inner_join(clin2, by = "sample_id")

n_merged <- nrow(dat)
message("  Samples with both module score and PFS: ", n_merged)

if (n_merged < 10) {
  warning("Fewer than 10 samples with complete PFS data; results will be very unstable.")
}

## PFS به ماه (برای خوانایی؛ واحد absolute مهم نیست)
dat <- dat %>%
  dplyr::mutate(
    PFS_months = PFS_days / 30.4375
  )

n_events <- sum(dat$PFS_status == 1, na.rm = TRUE)
message("  Number of PFS events: ", n_events)

if (n_events < 5) {
  warning("Fewer than 5 PFS events; Cox model may be unreliable.")
}

## ------------------------------------------------------------
## 6. Cox model: score_net_z continuous
## ------------------------------------------------------------
message("\n--- Step 4: Cox model (continuous module score) ---")

cox_fit <- survival::coxph(
  survival::Surv(PFS_months, PFS_status) ~ score_net_z,
  data = dat
)

cox_summ <- summary(cox_fit)

HR  <- cox_summ$coef[1, "exp(coef)"]
LCL <- cox_summ$conf.int[1, "lower .95"]
UCL <- cox_summ$conf.int[1, "upper .95"]
pval <- cox_summ$coef[1, "Pr(>|z|)"]

message(sprintf(
  "  Cox HR per 1 SD increase in score_net_z: %.2f (95%% CI %.2f–%.2f), p = %.3g",
  HR, LCL, UCL, pval
))

cox_out <- tibble::tibble(
  cohort          = "GSE130874",
  n_samples       = n_merged,
  n_events        = n_events,
  HR_score_net    = HR,
  CI95_lower      = LCL,
  CI95_upper      = UCL,
  p_value         = pval,
  time_unit       = "months",
  score_variable  = "score_net_z"
)

cox_out_path <- file.path(
  surv_dir,
  "GSE130874_PFS_cox_crossSpeciesBROAD_score_net.tsv"
)
readr::write_tsv(cox_out, cox_out_path)
message("  Saved Cox summary to:\n  ", cox_out_path)

## ------------------------------------------------------------
## 7. KM: high vs low by median module score
## ------------------------------------------------------------
message("\n--- Step 5: Kaplan–Meier (high vs low module score) ---")

median_score <- stats::median(dat$score_net_z, na.rm = TRUE)

dat_km <- dat %>%
  dplyr::mutate(
    score_group = dplyr::if_else(
      score_net_z >= median_score,
      "High module score",
      "Low module score"
    )
  )

## ensure factor order
dat_km <- dat_km %>%
  dplyr::mutate(
    score_group = factor(
      score_group,
      levels = c("Low module score", "High module score")
    )
  )

km_fit <- survival::survfit(
  survival::Surv(PFS_months, PFS_status) ~ score_group,
  data = dat_km
)

## log-rank p
lr_test <- survival::survdiff(
  survival::Surv(PFS_months, PFS_status) ~ score_group,
  data = dat_km
)
lr_p <- 1 - stats::pchisq(lr_test$chisq, df = length(lr_test$n) - 1)

message(sprintf("  Log-rank p (High vs Low): %.3g", lr_p))

km_meta_out <- dat_km %>%
  dplyr::select(
    sample_id, score_net_z, score_group,
    PFS_days, PFS_months, PFS_status
  )

km_meta_path <- file.path(
  surv_dir,
  "GSE130874_PFS_KM_high_vs_low_module_scores.tsv"
)
readr::write_tsv(km_meta_out, km_meta_path)
message("  Saved KM sample-level table to:\n  ", km_meta_path)

## KM figure
km_plot <- survminer::ggsurvplot(
  km_fit,
  data              = dat_km,
  risk.table        = TRUE,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = FALSE,
  legend.title      = "",
  legend.labs       = c("Low module score", "High module score"),
  xlab              = "Progression-free survival (months)",
  ylab              = "PFS probability",
  palette           = c("#6baed6", "#08519c"),
  risk.table.height = 0.22,
  ggtheme           = theme_bw(base_size = 11)
)

fig_km_png <- file.path(
  fig_dir,
  "Fig6C_GSE130874_PFS_KM_high_vs_low_module.png"
)
fig_km_pdf <- file.path(
  fig_dir,
  "Fig6C_GSE130874_PFS_KM_high_vs_low_module.pdf"
)

ggplot2::ggsave(
  filename = fig_km_png,
  plot     = km_plot$plot,
  width    = 6.5,
  height   = 5.5,
  dpi      = 400
)
ggplot2::ggsave(
  filename = fig_km_pdf,
  plot     = km_plot$plot,
  width    = 6.5,
  height   = 5.5
)

message("  Saved KM figure to:\n  ", fig_km_png, "\n  ", fig_km_pdf)

## ------------------------------------------------------------
## 8. Optional: module score vs response group
## ------------------------------------------------------------
if (has_response) {
  message("\n--- Step 6: Module scores vs response_group ---")
  
  ## only keep groups with ≥3 samples
  grp_counts <- dat_km %>%
    dplyr::group_by(response_group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  valid_groups <- grp_counts %>%
    dplyr::filter(n >= 3) %>%
    dplyr::pull(response_group)
  
  if (length(valid_groups) < 2) {
    warning(
      "response_group has fewer than 2 groups with at least 3 samples; ",
      "skipping response-based plots."
    )
  } else {
    dat_resp <- dat_km %>%
      dplyr::filter(response_group %in% valid_groups) %>%
      dplyr::mutate(response_group = droplevels(response_group))
    
    ## simple Wilcoxon (only meaningful for 2 groups)
    if (length(levels(dat_resp$response_group)) == 2) {
      w_test <- stats::wilcox.test(
        score_net_z ~ response_group,
        data = dat_resp
      )
      w_p <- w_test$p.value
      message(sprintf(
        "  Wilcoxon p (module score vs response_group): %.3g",
        w_p
      ))
    } else {
      w_p <- NA_real_
      message(
        "  >2 response groups; skipping Wilcoxon, only plotting distribution."
      )
    }
    
    resp_plot <- ggplot(dat_resp,
                        aes(x = response_group, y = score_net_z)) +
      geom_violin(trim = FALSE, fill = "#deebf7", color = "#08519c") +
      geom_boxplot(width = 0.25, outlier.shape = NA,
                   fill = "white", color = "#08519c") +
      geom_jitter(width = 0.12, alpha = 0.7, size = 2) +
      theme_bw(base_size = 11) +
      labs(
        x = "Response group",
        y = "Cross-species BROAD module score (z)",
        title = "GSE130874: module activity vs clinical response"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 20, hjust = 1)
      )
    
    fig_resp_png <- file.path(
      fig_dir,
      "Fig6D_GSE130874_module_score_vs_response_group.png"
    )
    fig_resp_pdf <- file.path(
      fig_dir,
      "Fig6D_GSE130874_module_score_vs_response_group.pdf"
    )
    
    ggplot2::ggsave(
      filename = fig_resp_png,
      plot     = resp_plot,
      width    = 5.5,
      height   = 4.8,
      dpi      = 400
    )
    ggplot2::ggsave(
      filename = fig_resp_pdf,
      plot     = resp_plot,
      width    = 5.5,
      height   = 4.8
    )
    
    message("  Saved response figure to:\n  ",
            fig_resp_png, "\n  ", fig_resp_pdf)
    
    ## summary table
    resp_summary <- dat_resp %>%
      dplyr::group_by(response_group) %>%
      dplyr::summarise(
        n           = dplyr::n(),
        mean_score  = mean(score_net_z, na.rm = TRUE),
        median_score = stats::median(score_net_z, na.rm = TRUE),
        .groups     = "drop"
      ) %>%
      dplyr::mutate(
        wilcox_p = w_p
      )
    
    resp_out_path <- file.path(
      surv_dir,
      "GSE130874_module_score_vs_response_group_summary.tsv"
    )
    readr::write_tsv(resp_summary, resp_out_path)
    message("  Saved response summary to:\n  ", resp_out_path)
  }
}

message("\n=== Phase 8A completed successfully (اگر فایل کلینیکال درست خورده باشد) ===")
