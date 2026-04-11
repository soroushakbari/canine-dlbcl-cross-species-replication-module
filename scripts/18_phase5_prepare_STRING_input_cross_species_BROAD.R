## =====================================================================
## Phase 5 (fallback): Prepare STRING input for human cross-species BROAD module
##   - Uses cross_species_module_BROAD_table.tsv from phase 17
##   - Produces a clean gene list with direction/logFC/adj.P.Val
##     to upload manually to string-db.org (Multiple proteins)
## =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) file.path(project_root, ...)

sig_dir <- path_proj("results", "tables", "signatures")
net_dir <- path_proj("results", "tables", "network")
dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1. خواندن جدول cross-species BROAD
## ---------------------------------------------------------------------

broad_table_path <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
if (!file.exists(broad_table_path)) {
  stop("BROAD table not found: ", broad_table_path,
       "\nRun 17_phase4D_define_cross_species_module_broad.R first.")
}

cross_broad <- readr::read_tsv(broad_table_path, show_col_types = FALSE)

## انتظار داریم این ستون‌ها وجود داشته باشند (از فاز 17):
needed_cols <- c(
  "human_symbol", "dog_symbol",
  "human_logFC", "human_adj.P",
  "human_dir_broad", "human_is_sig_broad",
  "dog_logFC", "dog_adj.P",
  "dog_dir_broad", "dog_is_sig_broad"
)

missing <- setdiff(needed_cols, names(cross_broad))
if (length(missing) > 0L) {
  stop("Missing expected columns in BROAD table: ",
       paste(missing, collapse = ", "))
}

## ---------------------------------------------------------------------
## 2. ساخت جدول human-centric برای STRING
## ---------------------------------------------------------------------

string_input <- cross_broad %>%
  dplyr::select(
    human_symbol,
    dog_symbol,
    human_logFC,
    human_adj.P,
    human_dir_broad,
    human_is_sig_broad,
    dog_logFC,
    dog_adj.P,
    dog_dir_broad,
    dog_is_sig_broad
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(dplyr::desc(abs(human_logFC)))


## جهت را به صورت ساده up/down برمی‌داریم (nonsig را حذف می‌کنیم)
string_input <- string_input %>%
  filter(human_dir_broad %in% c("up", "down"))

message("[STRING-input] Human BROAD cross-species rows: ", nrow(string_input))

## فایل نهایی برای STRING (فقط human_symbol + جهت + logFC)
string_input_min <- string_input %>%
  transmute(
    gene_symbol = human_symbol,
    direction   = human_dir_broad,
    human_logFC = human_logFC,
    human_adj.P = human_adj.P
  ) %>%
  distinct()

out_path_full <- file.path(
  net_dir,
  "STRING_input_human_cross_species_BROAD_full.tsv"
)
out_path_min <- file.path(
  net_dir,
  "STRING_input_human_cross_species_BROAD_minimal.tsv"
)

readr::write_tsv(string_input,     out_path_full)
readr::write_tsv(string_input_min, out_path_min)

message("[STRING-input] Full table written to: ", out_path_full)
message("[STRING-input] Minimal table (for web upload) written to: ", out_path_min)
message("=== [Phase 5 fallback - Prepared STRING input for human BROAD module] DONE ===")
