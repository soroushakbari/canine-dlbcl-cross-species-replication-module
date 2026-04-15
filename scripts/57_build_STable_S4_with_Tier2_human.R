## 57_build_STable_S4_with_Tier2_human.R
## Build revised Supplementary Table S4 with 4 sheets:
##  - BROAD_full_gene_list
##  - STRICT_core_full_gene_list
##  - Tier1_core_full_gene_list
##  - Tier2_human_GSE56315

message("=== Build revised STable S4 with Tier2 human sheet ===")

required_pkgs <- c("readr", "dplyr", "openxlsx", "stringr")
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
  library(openxlsx)
  library(stringr)
})

project_root <- normalizePath(
  "D:/Research/My Articles/DLBCL drug",
  winslash = "/",
  mustWork = TRUE
)

sig_dir <- file.path(project_root, "results", "tables", "signatures")
out_dir <- file.path(project_root, "results", "tables")

path_broad  <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
path_strict <- file.path(sig_dir, "cross_species_module_STRICT_core_table.tsv")
path_tier1  <- file.path(sig_dir, "cross_species_module_Tier1_core_table.tsv")
path_tier2  <- file.path(sig_dir, "human_signature_GSE56315_Tier2_broad.tsv")

for (p in c(path_broad, path_strict, path_tier1, path_tier2)) {
  if (!file.exists(p)) {
    stop("Required file not found:\n  ", p)
  }
}

message("Reading source tables ...")
broad  <- readr::read_tsv(path_broad,  show_col_types = FALSE)
strict <- readr::read_tsv(path_strict, show_col_types = FALSE)
tier1  <- readr::read_tsv(path_tier1,  show_col_types = FALSE)
tier2  <- readr::read_tsv(path_tier2,  show_col_types = FALSE)

## ---------- helper to detect columns ----------
get_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

## ---------- clean Tier2 sheet for readability ----------
gene_col <- get_col(tier2, c("gene_symbol", "human_symbol", "symbol", "gene"))
fc_col   <- get_col(tier2, c("logFC", "human_logFC", "log2FC"))
fdr_col  <- get_col(tier2, c("adj.P.Val", "padj", "FDR", "adj_p", "human_adj.P"))
dir_col  <- get_col(tier2, c("direction", "human_dir", "dir"))
sig_col  <- get_col(tier2, c("is_sig", "human_is_sig", "significant"))

tier2_clean <- tier2

if (!is.na(gene_col)) tier2_clean <- tier2_clean %>% mutate(gene_symbol = .data[[gene_col]])
if (!is.na(fc_col))   tier2_clean <- tier2_clean %>% mutate(log2FC = as.numeric(.data[[fc_col]]))
if (!is.na(fdr_col))  tier2_clean <- tier2_clean %>% mutate(adj.P.Val = as.numeric(.data[[fdr_col]]))
if (!is.na(sig_col))  tier2_clean <- tier2_clean %>% mutate(is_significant = .data[[sig_col]])
if (!is.na(dir_col))  tier2_clean <- tier2_clean %>% mutate(direction = .data[[dir_col]])

front_cols <- c("gene_symbol", "log2FC", "adj.P.Val", "is_significant", "direction")
front_cols <- front_cols[front_cols %in% names(tier2_clean)]

tier2_clean <- tier2_clean %>%
  select(any_of(front_cols), everything())

## ---------- build workbook ----------
out_xlsx <- file.path(out_dir, "STable_S4_full_gene_lists_with_Tier2.xlsx")

wb <- openxlsx::createWorkbook()

sheet_names <- c(
  "BROAD_full_gene_list",
  "STRICT_core_full_gene_list",
  "Tier1_core_full_gene_list",
  "Tier2_human_GSE56315"
)

sheet_tables <- list(
  broad,
  strict,
  tier1,
  tier2_clean
)

for (i in seq_along(sheet_names)) {
  openxlsx::addWorksheet(wb, sheet_names[i])
  openxlsx::writeData(wb, sheet = sheet_names[i], x = sheet_tables[[i]], withFilter = TRUE)
  openxlsx::freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(sheet_tables[[i]]), widths = "auto")
}

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Saved revised STable S4 to:")
message("  ", out_xlsx)
message("=== Done ===")