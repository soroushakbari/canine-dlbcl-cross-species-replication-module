## 56_build_Table3_core_CMap_drugs_forManuscript.R
## Build clean Table 3 from existing core CMap summary

message("=== Build Table 3 (clean): core CMap drugs & network summary ===")

required_pkgs <- c("readr", "dplyr", "stringr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) and re-run this script."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

## --------------------------------------------------------------------
## paths
## --------------------------------------------------------------------
project_root <- normalizePath(
  file.path(".."),
  winslash = "/",
  mustWork = TRUE
)
message("Project root: ", project_root)

tab_dir <- file.path(project_root, "results", "tables")

in_path  <- file.path(tab_dir, "Table3_core_CMap_drugs_network_summary.tsv")
out_path <- file.path(tab_dir, "Table3_core_CMap_drugs_network_summary_forManuscript.tsv")

if (!file.exists(in_path)) {
  stop("Input Table3 not found at:\n  ", in_path)
}

## --------------------------------------------------------------------
## read & clean
## --------------------------------------------------------------------
t3_raw <- readr::read_tsv(in_path, show_col_types = FALSE)

needed <- c(
  "drug_name", "MOA_category", "moa",
  "n_targets_in_module", "n_targets_in_hubs", "n_targets_in_hub_bottleneck",
  "genes_in_hubs", "genes_in_hub_bottleneck",
  "min_score.x", "rank_min", "percentile_min", "best_cell"
)

missing <- setdiff(needed, names(t3_raw))
if (length(missing) > 0) {
  stop(
    "Table3 is missing required columns: ",
    paste(missing, collapse = ", ")
  )
}

t3_clean <- t3_raw %>%
  ## حذف -666 از moa
  mutate(
    moa_clean = sapply(moa, function(x) {
      if (is.na(x)) return(NA_character_)
      parts <- strsplit(x, ";")[[1]]
      parts <- trimws(parts)
      parts <- parts[parts != "-666" & parts != ""]
      if (length(parts) == 0L) return(NA_character_)
      paste(unique(parts), collapse = "; ")
    }),
    ## اگر hub_bottleneck خالی بود، از hubs استفاده کن
    core_hub_targets = dplyr::coalesce(genes_in_hub_bottleneck, genes_in_hubs)
  ) %>%
  ## مرتب‌سازی بر اساس rank_min
  arrange(rank_min) %>%
  ## ساخت جدول نهایی برای مقاله
  transmute(
    Rank                     = rank_min,
    Drug                     = drug_name,
    `Mechanism class`        = MOA_category,
    `Detailed mechanism`     = moa_clean,
    `Module targets (n)`     = n_targets_in_module,
    `Core hub targets (n)`   = n_targets_in_hubs,
    `Core hub genes`         = core_hub_targets,
    `Minimum connectivity score` = round(min_score.x, 2),
    `Library percentile`     = sprintf("%.6f", percentile_min),
    `Best CMap cell line`    = best_cell
  )

readr::write_tsv(t3_clean, out_path)

message("Saved clean Table 3 to:")
message("  ", out_path)
message("=== Table 3 (clean) build completed successfully. ===")
