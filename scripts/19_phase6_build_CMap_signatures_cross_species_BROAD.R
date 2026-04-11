## =====================================================================
## Phase 6: Build CMap/LINCS disease signatures from cross-species BROAD module
##   - Input: cross_species_module_BROAD_table.tsv (phase 17)
##   - Output: plain-text up/down gene lists for CLUE.io (CMap)
##     e.g., CMap_cross_species_BROAD_up_top100.txt, down_top100.txt
## =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) file.path(project_root, ...)

sig_dir <- path_proj("results", "tables", "signatures")
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1. خواندن جدول cross-species BROAD
## ---------------------------------------------------------------------

broad_table_path <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
if (!file.exists(broad_table_path)) {
  stop("BROAD table not found: ", broad_table_path,
       "\nRun 17_phase4D_define_cross_species_module_broad.R first.")
}

cross_broad <- readr::read_tsv(broad_table_path, show_col_types = FALSE)

needed_cols <- c(
  "human_symbol",
  "human_logFC",
  "human_adj.P",
  "human_dir_broad",
  "human_is_sig_broad"
)

missing <- setdiff(needed_cols, names(cross_broad))
if (length(missing) > 0L) {
  stop("Missing columns in BROAD table: ", paste(missing, collapse = ", "))
}

## فقط ژن‌های up/down واقعی (نه nonsig)
up_tbl <- cross_broad %>%
  filter(human_dir_broad == "up", human_is_sig_broad) %>%
  arrange(desc(human_logFC))

down_tbl <- cross_broad %>%
  filter(human_dir_broad == "down", human_is_sig_broad) %>%
  arrange(human_logFC)  # down: logFC منفی‌تر در بالا

message("[CMap] Candidate up genes:   ", nrow(up_tbl))
message("[CMap] Candidate down genes: ", nrow(down_tbl))

if (nrow(up_tbl) == 0L | nrow(down_tbl) == 0L) {
  stop("[CMap] No significant up/down genes in BROAD table – something is wrong.")
}

## ---------------------------------------------------------------------
## 2. تابع برای نوشتن top-N امضاها
## ---------------------------------------------------------------------

write_topN_signature <- function(N) {
  n_up   <- min(N, nrow(up_tbl))
  n_down <- min(N, nrow(down_tbl))
  
  up_genes   <- up_tbl$human_symbol[1:n_up]   %>% unique()
  down_genes <- down_tbl$human_symbol[1:n_down] %>% unique()
  
  up_path <- file.path(sig_dir,
                       paste0("CMap_cross_species_BROAD_up_top", N, ".txt"))
  down_path <- file.path(sig_dir,
                         paste0("CMap_cross_species_BROAD_down_top", N, ".txt"))
  
  readr::write_lines(up_genes, up_path)
  readr::write_lines(down_genes, down_path)
  
  message("[CMap] N = ", N,
          " → up: ", length(up_genes), " genes written to: ", up_path)
  message("[CMap] N = ", N,
          " → down: ", length(down_genes), " genes written to: ", down_path)
}

## ---------------------------------------------------------------------
## 3. ساخت نسخه‌های top50/top100/top150
## ---------------------------------------------------------------------

Ns <- c(50, 100, 150)

for (N in Ns) {
  write_topN_signature(N)
}

## همچنین کل لیست up/down کامل BROAD (برای GSVA یا استفاده‌های دیگر)
full_up_path   <- file.path(sig_dir, "CMap_cross_species_BROAD_up_FULL.txt")
full_down_path <- file.path(sig_dir, "CMap_cross_species_BROAD_down_FULL.txt")

readr::write_lines(up_tbl$human_symbol %>% unique(),   full_up_path)
readr::write_lines(down_tbl$human_symbol %>% unique(), full_down_path)

message("[CMap] FULL up list   → ", full_up_path)
message("[CMap] FULL down list → ", full_down_path)

## خلاصه‌ی متنی
summary_path <- file.path(sig_dir, "CMap_cross_species_BROAD_signature_summary.txt")
sink(summary_path)
cat("CMap/LINCS disease signatures – cross-species BROAD module\n\n")
cat("Total candidate up genes:   ", nrow(up_tbl), "\n")
cat("Total candidate down genes: ", nrow(down_tbl), "\n\n")
cat("Top-N signatures written for N in {50, 100, 150}.\n")
sink()

message("[CMap] Summary written to: ", summary_path)
message("=== [Phase 6 - CMap signatures from cross-species BROAD] DONE ===")
