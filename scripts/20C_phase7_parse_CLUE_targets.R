## 20C_phase7_parse_CLUE_targets.R
## هدف:
##  - پارس کردن ستون‌های target از فایل‌های CLUE شناسایی شده
##  - ساخت یک فایل جدید برای استفاده در بعدی‌های Drug–Gene–Network analysis

message("=== Phase 7C: Parse CLUE target columns ===")

required_pkgs <- c("readr", "dplyr", "stringr", "tibble")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(\"", pkg, "\") and then re-run this script."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)

message("Project root: ", project_root)

# فایل‌های ورودی
file_query_result <- file.path(
  project_root, "results", "tables", "Drug", "CLUE_raw",
  "my_analysis.sig_queryl1k_tool.694fb6107ccc020013c1aa46", "arfs/TAG/query_result.gct"
)

file_ncs <- file.path(
  project_root, "results", "tables", "Drug", "CLUE_raw",
  "my_analysis.sig_queryl1k_tool.694fb6107ccc020013c1aa46", "ncs.gct"
)

# فایل خروجی
output_drug_targets <- file.path(
  project_root, "metadata", "drug_targets_from_CLUE.tsv"
)

if (!file.exists(file_query_result) || !file.exists(file_ncs)) {
  stop("One or both of the CLUE files not found. Check paths.")
}

# پارس کردن
parse_targets <- function(file) {
  message("Parsing file: ", file)
  
  # خواندن داده‌ها
  df <- readr::read_tsv(file, show_col_types = FALSE, skip = 2)
  
  # پیدا کردن ستون target
  target_cols <- names(df)[stringr::str_detect(names(df), "target")]
  
  if (length(target_cols) == 0) {
    return(NULL)
  }
  
  # ایجاد لیست دارو و target ژن‌ها
  targets <- df %>%
    dplyr::select(drug_name = `pert_name`, dplyr::all_of(target_cols)) %>%
    dplyr::mutate(
      targets = dplyr::coalesce(.data[[target_cols[1]]], .data[[target_cols[2]]]), 
      targets = stringr::str_trim(targets)
    ) %>%
    dplyr::filter(!is.na(targets) & targets != "")
  
  return(targets)
}

# پارس کردن هر فایل
targets_query_result <- parse_targets(file_query_result)
targets_ncs <- parse_targets(file_ncs)

# ترکیب داده‌ها
all_targets <- dplyr::bind_rows(targets_query_result, targets_ncs)

# ذخیره کردن نتیجه
readr::write_tsv(all_targets, output_drug_targets)

message("Saved target data to:\n  ", output_drug_targets)
message("=== Phase 7C completed. ===")
