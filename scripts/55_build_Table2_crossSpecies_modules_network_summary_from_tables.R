## 55_build_Table2_crossSpecies_modules_network_summary_from_tables.R
## Rebuild Table 2 directly from module tables + PPI nodes
## بدون وابستگی به gene-list های txt (read_gene_list و این داستان‌ها)

message("=== Build Table 2: cross-species modules & PPI summary (FROM TABLES) ===")

required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "purrr")

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
  library(tibble)
  library(stringr)
  library(purrr)
})

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath("D:/Research/My Articles/DLBCL drug",
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

sig_dir <- file.path(project_root, "results", "tables", "signatures")
net_dir <- file.path(project_root, "results", "tables", "network")
out_dir <- file.path(project_root, "results", "tables")

if (!dir.exists(sig_dir)) {
  stop("Signatures directory not found at:\n  ", sig_dir)
}
if (!dir.exists(net_dir)) {
  stop("Network directory not found at:\n  ", net_dir)
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## --------------------------------------------------------------------
## 2. read module tables (BROAD, STRICT, Tier1_core)
## --------------------------------------------------------------------
path_broad  <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
path_strict <- file.path(sig_dir, "cross_species_module_STRICT_core_table.tsv")
path_tier1  <- file.path(sig_dir, "cross_species_module_Tier1_core_table.tsv")

for (p in c(path_broad, path_strict, path_tier1)) {
  if (!file.exists(p)) {
    stop("Module table not found at:\n  ", p)
  }
}

message("Reading module tables ...")

mod_broad <- readr::read_tsv(path_broad, show_col_types = FALSE)
mod_strict <- readr::read_tsv(path_strict, show_col_types = FALSE)
mod_tier1 <- readr::read_tsv(path_tier1, show_col_types = FALSE)

## ستون‌های لازم در هر جدول
needed_cols <- c("human_symbol", "dog_symbol", "human_logFC", "dog_logFC")
missing_broad  <- setdiff(needed_cols, names(mod_broad))
missing_strict <- setdiff(needed_cols, names(mod_strict))
missing_tier1  <- setdiff(needed_cols, names(mod_tier1))

if (length(missing_broad) > 0) {
  stop("BROAD table missing columns: ", paste(missing_broad, collapse = ", "))
}
if (length(missing_strict) > 0) {
  stop("STRICT table missing columns: ", paste(missing_strict, collapse = ", "))
}
if (length(missing_tier1) > 0) {
  stop("Tier1_core table missing columns: ", paste(missing_tier1, collapse = ", "))
}

## --------------------------------------------------------------------
## 3. read PPI node table (only exists for BROAD; others will be intersected)
## --------------------------------------------------------------------
ppi_path <- file.path(net_dir, "PPI_cross_species_BROAD_nodes.tsv")
if (!file.exists(ppi_path)) {
  stop("PPI node table not found at:\n  ", ppi_path)
}

message("Reading PPI node table ...")
ppi_nodes <- readr::read_tsv(ppi_path, show_col_types = FALSE)

## حداقل ستون‌های لازم داخل PPI
needed_ppi_cols <- c("gene", "degree", "betweenness", "hub_type")
missing_ppi <- setdiff(needed_ppi_cols, names(ppi_nodes))
if (length(missing_ppi) > 0) {
  stop("PPI node table missing columns: ", paste(missing_ppi, collapse = ", "))
}

## --------------------------------------------------------------------
## 4. helper: summarise one module
## --------------------------------------------------------------------
summarise_module <- function(module_name, tbl, ppi_nodes) {
  # فقط ردیف‌هایی که gene_symbol دارند
  tbl_clean <- tbl %>%
    dplyr::filter(!is.na(human_symbol),
                  !is.na(dog_symbol))
  
  n_pairs <- nrow(tbl_clean)
  
  # تعداد ژن‌های منحصر به فرد
  n_human_genes <- tbl_clean %>%
    dplyr::summarise(n = dplyr::n_distinct(human_symbol)) %>%
    dplyr::pull(n)
  
  n_dog_genes <- tbl_clean %>%
    dplyr::summarise(n = dplyr::n_distinct(dog_symbol)) %>%
    dplyr::pull(n)
  
  # جهت (up/down) بر اساس sign of logFC
  tbl_dir <- tbl_clean %>%
    dplyr::mutate(
      human_dir = dplyr::case_when(
        human_logFC > 0 ~ "up",
        human_logFC < 0 ~ "down",
        TRUE ~ NA_character_
      ),
      dog_dir = dplyr::case_when(
        dog_logFC > 0 ~ "up",
        dog_logFC < 0 ~ "down",
        TRUE ~ NA_character_
      )
    )
  
  n_up_human   <- sum(tbl_dir$human_dir == "up",   na.rm = TRUE)
  n_down_human <- sum(tbl_dir$human_dir == "down", na.rm = TRUE)
  n_up_dog     <- sum(tbl_dir$dog_dir == "up",     na.rm = TRUE)
  n_down_dog   <- sum(tbl_dir$dog_dir == "down",   na.rm = TRUE)
  
  # اتصال به PPI (فقط بر اساس human_symbol در مقابل ستون 'gene' در PPI)
  module_genes <- unique(tbl_clean$human_symbol)
  
  ppi_sub <- ppi_nodes %>%
    dplyr::filter(gene %in% module_genes)
  
  n_in_PPI <- nrow(ppi_sub)
  frac_in_PPI <- ifelse(n_human_genes > 0,
                        n_in_PPI / n_human_genes,
                        NA_real_)
  
  # hub / bottleneck / hub-bottleneck
  n_hub <- ppi_sub %>%
    dplyr::filter(hub_type %in% c("hub", "hub_bottleneck")) %>%
    nrow()
  
  n_bottleneck <- ppi_sub %>%
    dplyr::filter(hub_type %in% c("bottleneck", "hub_bottleneck")) %>%
    nrow()
  
  n_hub_bottleneck <- ppi_sub %>%
    dplyr::filter(hub_type == "hub_bottleneck") %>%
    nrow()
  
  tibble::tibble(
    module_name      = module_name,
    n_pairs          = n_pairs,
    n_human_genes    = n_human_genes,
    n_dog_genes      = n_dog_genes,
    n_up_human       = n_up_human,
    n_down_human     = n_down_human,
    n_up_dog         = n_up_dog,
    n_down_dog       = n_down_dog,
    n_genes_in_PPI   = n_in_PPI,
    frac_genes_in_PPI = frac_in_PPI,
    n_hubs           = n_hub,
    n_bottlenecks    = n_bottleneck,
    n_hub_bottlenecks = n_hub_bottleneck
  )
}

## --------------------------------------------------------------------
## 5. apply to BROAD / STRICT / Tier1_core
## --------------------------------------------------------------------
message("Summarising modules ...")

row_broad  <- summarise_module("BROAD",      mod_broad,  ppi_nodes)
row_strict <- summarise_module("STRICT",     mod_strict, ppi_nodes)
row_tier1  <- summarise_module("Tier1_core", mod_tier1,  ppi_nodes)

table2 <- dplyr::bind_rows(row_broad, row_strict, row_tier1) %>%
  dplyr::mutate(
    frac_genes_in_PPI = round(frac_genes_in_PPI, 3)
  )

out_path <- file.path(out_dir, "Table2_crossSpecies_modules_network_summary.tsv")
readr::write_tsv(table2, out_path)

message("Saved Table 2 to:")
message("  ", out_path)
message("=== Build Table 2 (FROM TABLES) completed successfully. ===")
