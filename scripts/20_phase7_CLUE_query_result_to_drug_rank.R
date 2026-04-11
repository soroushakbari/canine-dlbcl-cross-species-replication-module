## 20_phase7_CLUE_query_result_to_drug_rank.R
## Parse CLUE query_result.gct as plain TSV and rank drugs by negative connectivity

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) file.path(project_root, ...)

drug_results_dir <- path_proj("results", "tables", "Drug")
dir.create(drug_results_dir, recursive = TRUE, showWarnings = FALSE)

## --------------------------------------------------------------
## 0) مسیر فایل query_result.gct
## --------------------------------------------------------------

clue_raw_root <- file.path(drug_results_dir, "CLUE_raw")

## اگر job جدید ساختی فقط این رشته رو عوض کن
clue_job_id  <- "my_analysis.sig_queryl1k_tool.694fb6107ccc020013c1aa46"

clue_job_dir <- file.path(clue_raw_root, clue_job_id)
arfs_tag_dir <- file.path(clue_job_dir, "arfs", "TAG")
query_gct_path <- file.path(arfs_tag_dir, "query_result.gct")

if (!file.exists(query_gct_path)) {
  stop("[CLUE] query_result.gct not found at: ", query_gct_path)
}

message("[CLUE] Using query_result.gct: ", query_gct_path)

## --------------------------------------------------------------
## 1) تشخیص rough فرمت GCT و خواندن دیتا
## --------------------------------------------------------------

## خط اول را برای چک کردن نسخه بخوانیم
first_line <- readLines(query_gct_path, n = 1)
message("[CLUE] First line of GCT: ", first_line)

## برای GCT 1.2/1.3 استاندارد، دو خط اول metadata هستند → skip = 2
## اگر بعداً دیدیم مشکل دارد، می‌شود تغییر داد، ولی الان همین را فرض می‌گیریم.
gct_df <- readr::read_tsv(
  query_gct_path,
  skip = 2,
  show_col_types = FALSE
)

message("[CLUE] Parsed GCT as TSV: ", nrow(gct_df), " rows × ",
        ncol(gct_df), " columns.")

message("[CLUE] Column names (first 20): ",
        paste(head(colnames(gct_df), 20), collapse = ", "))

## --------------------------------------------------------------
## 2) پیدا کردن ستون‌های کلیدی: دارو، نوع perturbation، score
## --------------------------------------------------------------

## candidate columns for score
score_candidates <- c(
  "norm_cs",
  "raw_cs",
  "cs",
  "wtcs"
)

score_cols_present <- intersect(score_candidates, colnames(gct_df))
if (length(score_cols_present) == 0L) {
  stop("[CLUE] Could not find any score column among: ",
       paste(score_candidates, collapse = ", "),
       "\nAvailable columns: ",
       paste(colnames(gct_df), collapse = ", "))
}

score_col <- score_cols_present[1L]
message("[CLUE] Using '", score_col, "' as score column.")

## candidate columns for drug name
name_candidates <- c("pert_iname", "pert_id", "pert_desc")
name_cols_present <- intersect(name_candidates, colnames(gct_df))
if (length(name_cols_present) == 0L) {
  stop("[CLUE] Could not find any drug-name column among: ",
       paste(name_candidates, collapse = ", "),
       "\nAvailable columns: ",
       paste(colnames(gct_df), collapse = ", "))
}
name_col <- name_cols_present[1L]
message("[CLUE] Using '", name_col, "' as drug-name column.")

## pert_type برای جدا کردن trt_cp
pert_type_col <- if ("pert_type" %in% colnames(gct_df)) "pert_type" else NA_character_

## cell line
cell_candidates <- c("cell_iname", "cell_id")
cell_cols_present <- intersect(cell_candidates, colnames(gct_df))
cell_col <- if (length(cell_cols_present) > 0L) cell_cols_present[1L] else NA_character_

## moa / target اگر موجود باشند
moa_col    <- if ("moa" %in% colnames(gct_df)) "moa" else NA_character_
target_col <- if ("target" %in% colnames(gct_df)) "target" else NA_character_

## fdr_q_nlog10 اگر هست
fdr_col <- if ("fdr_q_nlog10" %in% colnames(gct_df)) "fdr_q_nlog10" else NA_character_

## --------------------------------------------------------------
## 3) آماده‌سازی دیتافریم برای aggregation
## --------------------------------------------------------------

df <- gct_df %>%
  mutate(
    score    = suppressWarnings(as.numeric(.data[[score_col]])),
    pert_name = .data[[name_col]]
  )

if (!is.na(pert_type_col)) {
  df <- df %>%
    mutate(pert_type = .data[[pert_type_col]])
} else {
  df <- df %>% mutate(pert_type = NA_character_)
}

if (!is.na(cell_col)) {
  df <- df %>%
    mutate(cell_use = .data[[cell_col]])
} else {
  df <- df %>% mutate(cell_use = NA_character_)
}

if (!is.na(moa_col)) {
  df <- df %>%
    mutate(moa_use = .data[[moa_col]])
} else {
  df <- df %>% mutate(moa_use = NA_character_)
}

if (!is.na(target_col)) {
  df <- df %>%
    mutate(target_use = .data[[target_col]])
} else {
  df <- df %>% mutate(target_use = NA_character_)
}

if (!is.na(fdr_col)) {
  df <- df %>%
    mutate(fdr_q_nlog10 = suppressWarnings(as.numeric(.data[[fdr_col]])))
} else {
  df <- df %>% mutate(fdr_q_nlog10 = NA_real_)
}

## حذف ردیف‌هایی که score ندارند
df <- df %>%
  filter(!is.na(score))

message("[CLUE] Rows with non-NA scores: ", nrow(df))

## --------------------------------------------------------------
## 4) فیلتر به small molecules (trt_cp)، اگر ستونش هست
## --------------------------------------------------------------

if (!all(is.na(df$pert_type))) {
  message("[CLUE] pert_type distribution:")
  print(table(df$pert_type, useNA = "ifany"))
  df_drug <- df %>%
    filter(pert_type == "trt_cp")
  message("[CLUE] Keeping ", nrow(df_drug), " rows with pert_type == 'trt_cp'.")
} else {
  warning("[CLUE] 'pert_type' not found; keeping all rows.")
  df_drug <- df
}

if (nrow(df_drug) == 0L) {
  stop("[CLUE] No rows remain after filtering to trt_cp.")
}

## --------------------------------------------------------------
## 5) aggregation در سطح دارو
## --------------------------------------------------------------

drug_rank <- df_drug %>%
  group_by(pert_name) %>%
  summarise(
    n_signatures = n(),
    min_score    = min(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    mean_score   = mean(score, na.rm = TRUE),
    best_fdr_q_nlog10 = {
      if (all(is.na(fdr_q_nlog10))) {
        NA_real_
      } else {
        ## امتیاز min_score مربوط به کدوم signature بوده؟
        idx <- which.min(score)
        fdr_q_nlog10[idx]
      }
    },
    best_cell = {
      idx <- which.min(score)
      paste(unique(cell_use[idx]), collapse = ";")
    },
    moa = paste(
      unique(moa_use[!is.na(moa_use) & moa_use != ""]),
      collapse = ";"
    ),
    targets = paste(
      unique(target_use[!is.na(target_use) & target_use != ""]),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  arrange(min_score) %>%
  mutate(
    direction = ifelse(min_score < 0, "negative", "positive_or_zero")
  )

## --------------------------------------------------------------
## 6) ذخیره‌ی نتایج
## --------------------------------------------------------------

drug_rank_path <- file.path(
  drug_results_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_ranked.tsv"
)
readr::write_tsv(drug_rank, drug_rank_path)
message("[CLUE] Drug-level ranking written to: ", drug_rank_path)

top_neg <- drug_rank %>%
  filter(direction == "negative") %>%
  arrange(min_score) %>%
  head(50)

top_neg_path <- file.path(
  drug_results_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_top50_negative.tsv"
)
readr::write_tsv(top_neg, top_neg_path)
message("[CLUE] Top 50 negative-connectivity drugs written to: ", top_neg_path)

message("=== [Phase 7 - Parsed CLUE query_result.gct and ranked drugs] DONE ===")
