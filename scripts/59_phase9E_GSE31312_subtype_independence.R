## 59_phase9E_GSE31312_subtype_independence.R
## هدف:
##  - بررسی اینکه signal ماژول در GSE31312 مستقل از subtype (GCB vs ABC) هست یا نه
##  - 1) مقایسه module score بین GCB و ABC
##  - 2) Cox model با adjust برای subtype
##  - 3) ذخیره جدول‌ها و یک شکل supplementary ساده

message("=== Phase 9E: GSE31312 subtype independence analysis ===")

required_pkgs <- c(
  "readr", "dplyr", "stringr", "tibble",
  "survival", "ggplot2", "readxl", "purrr"
)

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
  library(survival)
  library(ggplot2)
  library(readxl)
  library(purrr)
})

project_root <- normalizePath(
  "D:/Research/My Articles/DLBCL drug",
  winslash = "/",
  mustWork = TRUE
)

surv_dir   <- file.path(project_root, "results", "tables", "survival")
fig_dir    <- file.path(project_root, "results", "figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

score_path <- file.path(surv_dir, "GSE31312_OS_scores_and_clinical_from_pdf.tsv")
if (!file.exists(score_path)) {
  stop("Score/OS merged file not found:\n  ", score_path)
}

message("Using score/OS file:")
message("  ", score_path)

score_df <- readr::read_tsv(score_path, show_col_types = FALSE)

message("Columns in score file:")
message("  ", paste(names(score_df), collapse = ", "))

## -------------------------------------------------------------------
## 1. helper functions
## -------------------------------------------------------------------
pick_col <- function(df, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    if (required) {
      stop("Could not detect ", label, ". Available columns:\n  ",
           paste(names(df), collapse = ", "))
    } else {
      return(NA_character_)
    }
  }
  hit[1]
}

clean_names_simple <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
}

extract_case_id <- function(x) {
  x <- as.character(x)
  out <- str_extract(x, "N\\d{4,}")
  out
}

read_candidate_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("tsv", "txt")) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
    return(df)
  }
  
  if (ext == "csv") {
    df <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
    return(df)
  }
  
  if (ext %in% c("xlsx", "xls")) {
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    if (length(sheets) == 0) return(NULL)
    
    ## اول شیتی را انتخاب کن که احتمال clinical sheet بودنش بیشتر است
    pref <- sheets[str_detect(tolower(sheets), "clinical|pdf|sheet|gse31312|31312")]
    sheets_order <- c(pref, setdiff(sheets, pref))
    
    for (sh in sheets_order) {
      df <- tryCatch(readxl::read_excel(path, sheet = sh), error = function(e) NULL)
      if (!is.null(df) && ncol(df) > 0 && nrow(df) > 0) {
        attr(df, "sheet_used") <- sh
        return(as_tibble(df))
      }
    }
    return(NULL)
  }
  
  NULL
}

find_clinical_sheet_with_subtype <- function(project_root) {
  cand_files <- list.files(
    project_root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(xlsx|xls|csv|tsv|txt)$"
  )
  
  ## فقط فایل‌های محتمل
  cand_files <- cand_files[
    str_detect(
      tolower(cand_files),
      "31312|clinical|pdf|depository|supp|metadata"
    )
  ]
  
  if (length(cand_files) == 0) return(NULL)
  
  message("Searching candidate files for subtype information ...")
  
  for (f in cand_files) {
    df <- read_candidate_file(f)
    if (is.null(df)) next
    
    nm_old <- names(df)
    nm_new <- clean_names_simple(nm_old)
    names(df) <- nm_new
    
    subtype_col <- pick_col(
      df,
      c("gep", "gene_expression_profiling_subgroup", "cell_of_origin", "subtype"),
      "subtype column",
      required = FALSE
    )
    
    case_col <- pick_col(
      df,
      c("geo_depository", "geo_depository_number", "case_geo_id", "geodepository", "geo_depository_"),
      "case id column",
      required = FALSE
    )
    
    if (!is.na(subtype_col) && !is.na(case_col)) {
      message("Found clinical subtype file:")
      message("  File : ", f)
      if (!is.null(attr(df, "sheet_used"))) {
        message("  Sheet: ", attr(df, "sheet_used"))
      }
      return(list(df = df, file = f, subtype_col = subtype_col, case_col = case_col))
    }
  }
  
  NULL
}

## -------------------------------------------------------------------
## 2. detect score/survival columns
## -------------------------------------------------------------------
score_col <- pick_col(
  score_df,
  c("score_net_z", "score_net", "module_score_z", "signature_zscore", "score"),
  "module score column"
)

os_time_col <- pick_col(
  score_df,
  c("OS_months", "OS", "os", "overall_survival_months"),
  "OS time column"
)

event_col <- pick_col(
  score_df,
  c("os_event", "OS_event", "event_os", "status_os", "event", "OScensor"),
  "OS event column"
)

case_id_col <- pick_col(
  score_df,
  c("case_geo_id", "case_id"),
  "case id column in score file"
)

## -------------------------------------------------------------------
## 3. find subtype file and merge
## -------------------------------------------------------------------
subtype_info <- find_clinical_sheet_with_subtype(project_root)

if (is.null(subtype_info)) {
  stop(
    "Could not automatically find a clinical sheet with subtype information.\n",
    "Please place the GSE31312 PDF-derived clinical sheet (with GEP column and GEO Depository #) inside the project folder,\n",
    "then re-run this script."
  )
}

clin_df <- subtype_info$df
clin_subtype_col <- subtype_info$subtype_col
clin_case_col    <- subtype_info$case_col

clin2 <- clin_df %>%
  mutate(
    case_geo_id = extract_case_id(.data[[clin_case_col]]),
    subtype_raw = as.character(.data[[clin_subtype_col]])
  ) %>%
  filter(!is.na(case_geo_id), !is.na(subtype_raw), subtype_raw != "") %>%
  distinct(case_geo_id, .keep_all = TRUE) %>%
  mutate(
    subtype = case_when(
      str_detect(str_to_upper(subtype_raw), "ABC") ~ "ABC",
      str_detect(str_to_upper(subtype_raw), "NON.?GCB") ~ "ABC",
      str_detect(str_to_upper(subtype_raw), "GCB") ~ "GCB",
      str_detect(str_to_upper(subtype_raw), "GERMINAL") ~ "GCB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(subtype)) %>%
  select(case_geo_id, subtype, subtype_raw)

df_sub <- score_df %>%
  mutate(
    case_geo_id = as.character(.data[[case_id_col]]),
    score_use   = as.numeric(.data[[score_col]]),
    os_time     = as.numeric(.data[[os_time_col]]),
    os_event    = as.numeric(.data[[event_col]])
  ) %>%
  left_join(clin2, by = "case_geo_id") %>%
  filter(
    !is.na(score_use),
    !is.na(os_time),
    !is.na(os_event),
    !is.na(subtype)
  ) %>%
  mutate(
    subtype = factor(subtype, levels = c("GCB", "ABC"))
  )

message("N with score + subtype + OS:")
message("  ", nrow(df_sub))
message("Subtype counts:")
print(table(df_sub$subtype))

## -------------------------------------------------------------------
## 4. subtype score comparison
## -------------------------------------------------------------------
subtype_summary <- df_sub %>%
  group_by(subtype) %>%
  summarise(
    n = n(),
    mean_score = mean(score_use, na.rm = TRUE),
    median_score = median(score_use, na.rm = TRUE),
    sd_score = sd(score_use, na.rm = TRUE),
    .groups = "drop"
  )

tt_sub <- t.test(score_use ~ subtype, data = df_sub)
wil_sub <- wilcox.test(score_use ~ subtype, data = df_sub, exact = FALSE)

subtype_test <- tibble(
  score_column = score_col,
  subtype_column = clin_subtype_col,
  n_total = nrow(df_sub),
  n_GCB = sum(df_sub$subtype == "GCB"),
  n_ABC = sum(df_sub$subtype == "ABC"),
  ttest_p = tt_sub$p.value,
  wilcox_p = wil_sub$p.value,
  mean_GCB = mean(df_sub$score_use[df_sub$subtype == "GCB"]),
  mean_ABC = mean(df_sub$score_use[df_sub$subtype == "ABC"])
)

## -------------------------------------------------------------------
## 5. subtype-adjusted Cox model
## -------------------------------------------------------------------
cox_fit <- coxph(Surv(os_time, os_event) ~ score_use + subtype, data = df_sub)
cox_sum <- summary(cox_fit)

coef_df <- as.data.frame(cox_sum$coefficients)
ci_df   <- as.data.frame(cox_sum$conf.int)

coef_df$term <- rownames(coef_df)
ci_df$term   <- rownames(ci_df)

cox_out <- coef_df %>%
  left_join(ci_df, by = "term") %>%
  transmute(
    term = term,
    coef = coef,
    HR = `exp(coef).y`,
    CI_lower = `lower .95`,
    CI_upper = `upper .95`,
    p_value = `Pr(>|z|)`
  )

## -------------------------------------------------------------------
## 6. supplementary figure
## -------------------------------------------------------------------
p_subtype <- ggplot(df_sub, aes(x = subtype, y = score_use, fill = subtype)) +
  geom_violin(trim = FALSE, alpha = 0.7, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.6) +
  labs(
    x = "DLBCL subtype",
    y = "Tier2 module score"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(df_sub$score_use, na.rm = TRUE),
    label = paste0("Wilcoxon p = ", format(wil_sub$p.value, digits = 3, scientific = TRUE)),
    vjust = -0.5,
    size = 4
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

fig_path_png <- file.path(fig_dir, "SFig_GSE31312_Tier2_module_by_subtype.png")
fig_path_pdf <- file.path(fig_dir, "SFig_GSE31312_Tier2_module_by_subtype.pdf")

ggsave(fig_path_png, p_subtype, width = 4.8, height = 4.5, dpi = 300)
ggsave(fig_path_pdf, p_subtype, width = 4.8, height = 4.5)

## -------------------------------------------------------------------
## 7. write outputs
## -------------------------------------------------------------------
out_summary <- file.path(surv_dir, "GSE31312_subtype_module_score_summary.tsv")
out_test    <- file.path(surv_dir, "GSE31312_subtype_module_score_tests.tsv")
out_cox     <- file.path(surv_dir, "GSE31312_subtype_adjusted_cox.tsv")

readr::write_tsv(subtype_summary, out_summary)
readr::write_tsv(subtype_test, out_test)
readr::write_tsv(cox_out, out_cox)

message("\nSubtype summary:")
print(subtype_summary)

message("\nSubtype comparison tests:")
print(subtype_test)

message("\nSubtype-adjusted Cox:")
print(cox_out)

message("\nSaved files:")
message("  ", out_summary)
message("  ", out_test)
message("  ", out_cox)
message("  ", fig_path_png)
message("  ", fig_path_pdf)

message("=== Done ===")