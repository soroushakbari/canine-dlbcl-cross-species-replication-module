## 42_phase10C_GSE171272_EVmiRNA_targets_to_module_PPI.R
## Map significant EV-miRNAs (GSE171272) to cross-species BROAD module & PPI hubs using multiMiR

message("=== Phase 10C: EV-miRNA targets -> cross-species module & PPI hubs ===")

required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "multiMiR"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run:\ninstall.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      "))\nthen re-run this script."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(multiMiR)
})

## ------------------------------------------------------------------
## paths
## ------------------------------------------------------------------
project_root <- normalizePath(file.path(".."), winslash = "/", mustWork = TRUE)
message("Project root: ", project_root)

sig_mir_dir <- file.path(project_root, "results", "tables", "miRNA")
sig_mir_path <- file.path(
  sig_mir_dir,
  "GSE171272_EVmiRNA_DE_sig_FDR0.05_logFC1.tsv"
)

sig_dir <- file.path(project_root, "results", "tables", "signatures")
net_dir <- file.path(project_root, "results", "tables", "network")

if (!file.exists(sig_mir_path)) {
  stop("Significant EV-miRNA DE file not found at:\n  ", sig_mir_path,
       "\nRun Phase 10B first.")
}

## cross-species BROAD module table
broad_path <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
if (!file.exists(broad_path)) {
  stop("BROAD table not found at:\n  ", broad_path)
}

## STRICT core table (optional)
strict_core_path <- file.path(sig_dir, "cross_species_module_STRICT_core_table.tsv")
has_strict_core <- file.exists(strict_core_path)

## PPI hubs/nodes
ppi_nodes_path <- file.path(net_dir, "PPI_cross_species_BROAD_nodes.tsv")
ppi_hubs_path  <- file.path(net_dir, "PPI_cross_species_BROAD_hubs.tsv")

if (!file.exists(ppi_nodes_path)) {
  stop("PPI node table not found at:\n  ", ppi_nodes_path)
}
if (!file.exists(ppi_hubs_path)) {
  stop("PPI hubs table not found at:\n  ", ppi_hubs_path)
}

## ------------------------------------------------------------------
## 1. Load sig miRNAs + module/PPI info
## ------------------------------------------------------------------
mir_sig <- readr::read_tsv(sig_mir_path, show_col_types = FALSE)

if (!"miRNA_id" %in% names(mir_sig)) {
  stop("DE table must contain a 'miRNA_id' column.")
}

message("  Significant EV-miRNAs: ", nrow(mir_sig))

mir_list <- mir_sig %>%
  dplyr::pull(miRNA_id) %>%
  unique() %>%
  sort()

message("  Unique miRNAs passed to multiMiR: ", length(mir_list))

## cross-species BROAD table
broad_tbl <- readr::read_tsv(broad_path, show_col_types = FALSE)

if (!"human_symbol" %in% names(broad_tbl)) {
  stop("BROAD table must contain a 'human_symbol' column.")
}

broad_small <- broad_tbl %>%
  dplyr::select(human_symbol, dplyr::everything()) %>%
  dplyr::distinct(human_symbol, .keep_all = TRUE)

broad_genes <- unique(broad_small$human_symbol)

## STRICT core (optional)
if (has_strict_core) {
  strict_core_tbl <- readr::read_tsv(strict_core_path, show_col_types = FALSE)
  
  if (!"human_symbol" %in% names(strict_core_tbl)) {
    warning("STRICT core table lacks 'human_symbol'; will ignore STRICT core flag.")
    has_strict_core <- FALSE
  } else {
    strict_core_genes <- unique(strict_core_tbl$human_symbol)
  }
}

## PPI nodes + hubs
ppi_nodes <- readr::read_tsv(ppi_nodes_path, show_col_types = FALSE)
ppi_hubs  <- readr::read_tsv(ppi_hubs_path,  show_col_types = FALSE)

## gene column detection
node_gene_col <- dplyr::case_when(
  "gene"   %in% names(ppi_nodes) ~ "gene",
  "symbol" %in% names(ppi_nodes) ~ "symbol",
  "node"   %in% names(ppi_nodes) ~ "node",
  TRUE ~ NA_character_
)

hub_gene_col <- dplyr::case_when(
  "gene"   %in% names(ppi_hubs) ~ "gene",
  "symbol" %in% names(ppi_hubs) ~ "symbol",
  "node"   %in% names(ppi_hubs) ~ "node",
  TRUE ~ NA_character_
)

if (is.na(node_gene_col) || is.na(hub_gene_col)) {
  stop("Could not detect gene column in PPI node / hub tables.")
}

ppi_nodes_small <- ppi_nodes %>%
  dplyr::rename(gene = !!node_gene_col) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

ppi_hubs_small <- ppi_hubs %>%
  dplyr::rename(gene = !!hub_gene_col) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

ppi_core_genes <- unique(ppi_hubs_small$gene)

## ------------------------------------------------------------------
## 2. Query multiMiR for validated targets
## ------------------------------------------------------------------
message("--- Querying multiMiR (validated targets only) ---")
mm_obj <- try(
  multiMiR::get_multimir(
    mirna      = mir_list,
    table      = "validated",
    legacy.out = FALSE
  ),
  silent = TRUE
)

if (inherits(mm_obj, "try-error")) {
  stop(
    "multiMiR::get_multimir failed. Possible causes:\n",
    "  - No internet connection\n",
    "  - multiMiR DB/server unavailable\n",
    "  - Firewall/proxy blocking access\n",
    "Original error:\n  ", as.character(mm_obj)
  )
}

mm_df <- as.data.frame(mm_obj@data)

## حداقل‌ها: miRNA + یک ستون target
if (!"mature_mirna_id" %in% names(mm_df)) {
  stop("multiMiR result lacks 'mature_mirna_id' column.")
}

## کشف نام ستون تارگت
target_col <- dplyr::case_when(
  "target_symbol" %in% names(mm_df) ~ "target_symbol",
  "target_gene"   %in% names(mm_df) ~ "target_gene",
  "target"        %in% names(mm_df) ~ "target",
  TRUE ~ NA_character_
)

if (is.na(target_col)) {
  stop(
    "Could not find a target-gene column in multiMiR result.\n",
    "Available columns:\n  ", paste(names(mm_df), collapse = ", ")
  )
}

## یک کپی با ستون استاندارد target_symbol بساز
mm_df2 <- mm_df
mm_df2$target_symbol <- as.character(mm_df2[[target_col]])

mm_clean <- mm_df2 %>%
  dplyr::filter(!is.na(target_symbol), target_symbol != "") %>%
  dplyr::mutate(
    mature_mirna_id = stringr::str_squish(mature_mirna_id),
    target_symbol   = stringr::str_squish(target_symbol)
  )

message("  Validated miRNA–target pairs from multiMiR: ", nrow(mm_clean))

out_all_targets <- file.path(
  sig_mir_dir,
  "GSE171272_EVmiRNA_multiMiR_validated_targets_all.tsv"
)
readr::write_tsv(mm_clean, out_all_targets)
message("  Saved full validated multiMiR result to:\n  ", out_all_targets)

## ------------------------------------------------------------------
## 3. Filter to BROAD module / PPI hubs & annotate
## ------------------------------------------------------------------
core_genes_of_interest <- c("TOP2A", "PARP1")

targets_filtered <- mm_clean %>%
  dplyr::mutate(
    in_BROAD_module = target_symbol %in% broad_genes,
    in_PPI_core     = target_symbol %in% ppi_core_genes,
    hits_TOP2A      = target_symbol == "TOP2A",
    hits_PARP1      = target_symbol == "PARP1"
  ) %>%
  dplyr::filter(in_BROAD_module | in_PPI_core | hits_TOP2A | hits_PARP1)

message(
  "  Validated pairs hitting BROAD module and/or PPI core: ",
  nrow(targets_filtered)
)

targets_annot <- targets_filtered %>%
  dplyr::left_join(
    broad_small %>% dplyr::select(human_symbol, dplyr::everything()),
    by = c("target_symbol" = "human_symbol")
  ) %>%
  dplyr::left_join(
    ppi_nodes_small %>% dplyr::select(gene, dplyr::everything()),
    by = c("target_symbol" = "gene"),
    suffix = c("", "_PPI")
  ) %>%
  dplyr::left_join(
    ppi_hubs_small %>%
      dplyr::select(gene, hub_type = dplyr::matches("hub_type|class|type")),
    by = c("target_symbol" = "gene")
  )

if (has_strict_core) {
  targets_annot <- targets_annot %>%
    dplyr::mutate(in_STRICT_core = target_symbol %in% strict_core_genes)
}

out_targets_module <- file.path(
  sig_mir_dir,
  "GSE171272_EVmiRNA_targets_in_BROAD_module_PPI.tsv"
)
readr::write_tsv(targets_annot, out_targets_module)
message("  Saved module/PPI-annotated targets to:\n  ", out_targets_module)

## ------------------------------------------------------------------
## 4. Summaries per miRNA
## ------------------------------------------------------------------
summary_by_miR <- targets_annot %>%
  dplyr::group_by(mature_mirna_id) %>%
  dplyr::summarise(
    n_targets_total     = dplyr::n_distinct(target_symbol),
    n_targets_BROAD     = dplyr::n_distinct(target_symbol[in_BROAD_module]),
    n_targets_PPI_core  = dplyr::n_distinct(target_symbol[in_PPI_core]),
    n_hits_TOP2A        = sum(hits_TOP2A),
    n_hits_PARP1        = sum(hits_PARP1),
    hits_any_core_hub   = any(in_PPI_core),
    hits_TOP2A_or_PARP1 = any(hits_TOP2A | hits_PARP1),
    .groups             = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_targets_PPI_core), dplyr::desc(n_targets_BROAD))

## add logFC / adj.P.Val
summary_by_miR <- summary_by_miR %>%
  dplyr::left_join(
    mir_sig %>% dplyr::select(miRNA_id, logFC, adj.P.Val),
    by = c("mature_mirna_id" = "miRNA_id")
  )

out_summary_miR <- file.path(
  sig_mir_dir,
  "GSE171272_EVmiRNA_targets_module_summary_by_miRNA.tsv"
)
readr::write_tsv(summary_by_miR, out_summary_miR)
message("  Saved per-miRNA summary to:\n  ", out_summary_miR)

## ------------------------------------------------------------------
## 5. Console snapshot
## ------------------------------------------------------------------
message("\n--- Top EV-miRNAs by PPI-core hits ---")
print(
  summary_by_miR %>%
    head(15)
)

message("=== Phase 10C completed. Use the summary table to pick EV-miRNAs for the network figure (Phase 10D). ===")
