## ==========================================
## 12_phase4E_TCGA_survival_analysis.R
## Phase 4E - Survival analysis in TCGA-DLBC using Tier2 module score
## Inputs:
##   results/tables/module_scores/module_scores_TCGA_DLBC_Tier2_GSE56315_signature_zscore.tsv
##   metadata/TCGA_DLBC_sample_metadata.csv
## Outputs:
##   results/tables/survival/survival_summary_TCGA_DLBC_Tier2_GSE56315.txt
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"

module_dir <- file.path(project_root, "results", "tables", "module_scores")
metadata_dir <- file.path(project_root, "metadata")
out_dir    <- file.path(project_root, "results", "tables", "survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------- 0) Packages ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkg_bioc <- c("survival")
pkg_cran <- c("tidyverse")

for (p in pkg_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

for (p in pkg_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

## ---------- 1) Load scores ----------

scores_file <- file.path(
  module_dir,
  "module_scores_TCGA_DLBC_Tier2_GSE56315_signature_zscore.tsv"
)

if (!file.exists(scores_file)) {
  stop("Scores file not found for TCGA-DLBC:\n  ", scores_file,
       "\nRun 11_phase4D_TCGA_module_scores_Tier2.R first.")
}

scores <- readr::read_tsv(scores_file, show_col_types = FALSE)

if (!all(c("sample_id", "score_net") %in% colnames(scores))) {
  stop("Scores file must contain 'sample_id' and 'score_net'.")
}

message("Loaded ", nrow(scores), " TCGA-DLBC samples with module scores.")

## ---------- 2) Load metadata & match sample_id ----------

meta_file <- file.path(metadata_dir, "TCGA_DLBC_sample_metadata.csv")

if (!file.exists(meta_file)) {
  stop("TCGA sample metadata file not found:\n  ", meta_file)
}

meta <- readr::read_csv(meta_file, show_col_types = FALSE)

cn_meta <- colnames(meta)

if (!"sample_id" %in% cn_meta) {
  message("No 'sample_id' column in metadata. Trying to detect best matching column...")
  matches <- sapply(cn_meta, function(col) {
    length(intersect(as.character(meta[[col]]), scores$sample_id))
  })
  print(matches)
  
  best_idx <- which.max(matches)
  if (length(best_idx) == 0 || matches[best_idx] == 0) {
    stop("Could not find any metadata column matching scores$sample_id.\n",
         "You may need to add a 'sample_id' column to metadata manually.")
  }
  
  best_col <- cn_meta[best_idx]
  message("Using metadata column ", best_col, " as sample_id.")
  meta <- meta %>%
    dplyr::rename(sample_id = !!best_col)
} else {
  message("Metadata already contains 'sample_id' column.")
}

## ---------- 3) Merge scores + metadata ----------

dat <- scores %>%
  dplyr::left_join(meta, by = "sample_id")

message("Merged scores+metadata: ", nrow(dat), " rows, ",
        ncol(dat), " columns.")

## ---------- 4) Build TCGA-style time & event ----------

cn <- colnames(dat)

# کاندید برای days_to_death و days_to_last_followup
death_cols <- grep("^days_to_death$", cn, ignore.case = TRUE, value = TRUE)
follow_cols <- grep("days_to_last_follow[_]*up", cn, ignore.case = TRUE, value = TRUE)

time_days <- rep(NA_real_, nrow(dat))

if (length(death_cols) > 0 || length(follow_cols) > 0) {
  death_vec <- rep(NA_real_, nrow(dat))
  follow_vec <- rep(NA_real_, nrow(dat))
  
  if (length(death_cols) > 0) {
    death_vec <- suppressWarnings(as.numeric(dat[[death_cols[1]]]))
  }
  if (length(follow_cols) > 0) {
    follow_vec <- suppressWarnings(as.numeric(dat[[follow_cols[1]]]))
  }
  
  time_days <- ifelse(!is.na(death_vec), death_vec, follow_vec)
  message("Built time_days from TCGA columns: ",
          paste(c(death_cols, follow_cols), collapse = ", "))
} else {
  stop("Could not find any 'days_to_death' or 'days_to_last_followup' columns in TCGA metadata.")
}

# event از vital_status یا status
event_cols <- grep("vital_status|status", cn, ignore.case = TRUE, value = TRUE)

if (length(event_cols) == 0) {
  stop("Could not find any 'vital_status' or 'status' columns in TCGA metadata.")
}

event_col <- event_cols[1]
event_raw <- dat[[event_col]]

map_event_to01 <- function(x) {
  v <- as.character(x)
  v_trim <- trimws(tolower(v))
  res <- rep(NA_real_, length(v_trim))
  
  # حالت Alive/Dead
  res[v_trim %in% c("alive", "censored", "no", "0", "n")]  <- 0
  res[v_trim %in% c("dead", "deceased", "yes", "1", "y")] <- 1
  
  # اگر numeric هم باشد
  suppressWarnings({
    vn <- as.numeric(v_trim)
  })
  if (any(!is.na(vn))) {
    u <- sort(unique(vn[!is.na(vn)]))
    if (all(u %in% c(0, 1))) {
      res <- vn
    } else if (all(u %in% c(1, 2))) {
      res <- ifelse(vn == 2, 1,
                    ifelse(vn == 1, 0, NA_real_))
    }
  }
  
  res
}

event_num <- map_event_to01(event_raw)

## ---------- 5) Summaries ----------

time_months <- time_days / 30.44

message("\nSummary of time_months:")
print(summary(time_months))
message("Summary of event_num (0=censored, 1=event, NA=unknown):")
print(table(event_num, useNA = "ifany"))

## ---------- 6) Clean data for survival analysis ----------

df_surv <- dat %>%
  dplyr::mutate(
    time  = time_months,
    event = event_num
  ) %>%
  dplyr::filter(
    !is.na(time),
    !is.na(event),
    time > 0
  )

message("\nAfter cleaning, ", nrow(df_surv), " samples have valid time & event.")

if (nrow(df_surv) < 20) {
  stop("Too few samples with valid survival data (n<20). ",
       "Check TCGA metadata or time/event mapping.")
}

## ---------- 7) Cox model with score_net (continuous) ----------

df_surv <- df_surv %>%
  dplyr::mutate(
    score_net_z = as.numeric(scale(score_net))
  )

cox_fit <- survival::coxph(
  survival::Surv(time, event) ~ score_net_z,
  data = df_surv
)

cox_summary <- summary(cox_fit)

hr      <- cox_summary$coefficients[,"exp(coef)"][1]
hr_low  <- cox_summary$conf.int[,"lower .95"][1]
hr_high <- cox_summary$conf.int[,"upper .95"][1]
p_val   <- cox_summary$wald["pvalue"]

message("\nCox model (score_net_z):")
message("  HR = ", round(hr, 3),
        " (95% CI ", round(hr_low, 3), "–", round(hr_high, 3), ")")
message("  Wald p-value = ", signif(p_val, 3))

## ---------- 8) KM high vs low (median split on score_net_z) ----------

median_cut <- median(df_surv$score_net_z, na.rm = TRUE)

df_surv <- df_surv %>%
  dplyr::mutate(
    group_score = ifelse(score_net_z >= median_cut, "High", "Low")
  )

message("\nGroup sizes (High vs Low):")
print(table(df_surv$group_score))

surv_fit_groups <- survival::survdiff(
  survival::Surv(time, event) ~ group_score,
  data = df_surv
)

chisq_val <- surv_fit_groups$chisq
df_chisq  <- length(surv_fit_groups$n) - 1
p_logrank <- 1 - pchisq(chisq_val, df = df_chisq)

message("Log-rank test (High vs Low score_net_z):")
message("  chisq = ", round(chisq_val, 3),
        ", df = ", df_chisq,
        ", p = ", signif(p_logrank, 3))

## ---------- 9) Save text summary ----------

out_txt <- file.path(out_dir, "survival_summary_TCGA_DLBC_Tier2_GSE56315.txt")

con <- file(out_txt, open = "wt", encoding = "UTF-8")

writeLines(c(
  "Survival analysis for TCGA-DLBC using Human_DLBCL_Tier2_GSE56315 signature",
  "",
  "Time (months) summary:",
  paste(capture.output(print(summary(time_months))), collapse = "\n"),
  "",
  paste0("Event column used: ", event_col),
  "Event summary (0=censored, 1=event):",
  paste(capture.output(print(table(event_num, useNA = "ifany"))), collapse = "\n"),
  "",
  paste0("N with valid survival data: ", nrow(df_surv)),
  "",
  "Cox model: Surv(time, event) ~ score_net_z",
  paste0("  HR (score_net_z) = ", round(hr, 3),
         " (95% CI ", round(hr_low, 3), "–", round(hr_high, 3), ")"),
  paste0("  Wald p-value = ", signif(p_val, 4)),
  "",
  "KM (median split on score_net_z):",
  paste(capture.output(print(table(df_surv$group_score))), collapse = "\n"),
  paste0("  Log-rank chisq = ", round(chisq_val, 3),
         ", df = ", df_chisq,
         ", p = ", signif(p_logrank, 4))
), con = con)

close(con)

message("\n[Phase 4E] Survival analysis for TCGA-DLBC completed.")
message("Summary written to:\n  ", out_txt, "\n")
