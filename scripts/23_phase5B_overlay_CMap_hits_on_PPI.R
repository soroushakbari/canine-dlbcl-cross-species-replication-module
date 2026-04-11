## =====================================================================
## Phase 5B: Overlay top CMap hits on PPI hubs
##   - Inputs:
##       results/tables/network/PPI_human_cross_species_BROAD_nodes_with_centrality.tsv
##       results/tables/Drug/CMap_queryl1k_cross_species_BROAD_drug_top50_negative_for_manuscript.tsv
##   - Outputs:
##       results/tables/Drug/CMap_PPI_overlay_cross_species_BROAD_all.tsv
##       results/tables/Drug/CMap_PPI_overlay_cross_species_BROAD_prioritized.tsv
## =====================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) file.path(project_root, ...)

net_dir  <- path_proj("results", "tables", "network")
drug_dir <- path_proj("results", "tables", "Drug")

nodes_cent_path <- file.path(
  net_dir,
  "PPI_human_cross_species_BROAD_nodes_with_centrality.tsv"
)

drug_top_path <- file.path(
  drug_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_top50_negative_for_manuscript.tsv"
)

if (!file.exists(nodes_cent_path)) {
  stop("[Overlay] Node+centrality file not found: ", nodes_cent_path)
}
if (!file.exists(drug_top_path)) {
  stop("[Overlay] Drug top50 file not found: ", drug_top_path)
}

nodes_cent <- readr::read_tsv(nodes_cent_path, show_col_types = FALSE)
drugs      <- readr::read_tsv(drug_top_path,   show_col_types = FALSE)

message("[Overlay] Loaded nodes: ", nrow(nodes_cent), " rows.")
message("[Overlay] Loaded drugs: ", nrow(drugs), " rows.")

## ---------------------------------------------------------------------
## 1. Prep node info
## ---------------------------------------------------------------------

if (!"human_symbol" %in% colnames(nodes_cent)) {
  stop("[Overlay] 'human_symbol' not found in nodes table.")
}

nodes_cent <- nodes_cent %>%
  dplyr::mutate(
    gene_symbol = human_symbol,
    gene_upper  = toupper(human_symbol),
    is_hub      = ifelse(is.na(is_hub), FALSE, is_hub)
  )

gene_set      <- unique(nodes_cent$gene_upper)
hub_gene_set  <- unique(nodes_cent$gene_upper[nodes_cent$is_hub %in% TRUE])

## ---------------------------------------------------------------------
## 2. Prep drug table
## ---------------------------------------------------------------------

## انتظار: ستون‌های زیر وجود دارند:
## rank, pert_name, MOA_category, moa, targets, n_signatures,
## min_score, median_score, mean_score, best_fdr_q_nlog10, best_cell

required_cols <- c(
  "rank", "pert_name", "MOA_category", "moa", "targets",
  "n_signatures", "min_score", "median_score", "mean_score",
  "best_fdr_q_nlog10", "best_cell"
)

missing_cols <- setdiff(required_cols, colnames(drugs))
if (length(missing_cols) > 0L) {
  stop("[Overlay] Missing columns in drug table: ",
       paste(missing_cols, collapse = ", "))
}

drugs <- drugs %>%
  dplyr::mutate(
    targets_raw = ifelse(is.na(targets), "", targets)
  )

## helper برای split targets به ژن‌ها
split_targets <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  parts <- unlist(strsplit(x, ";"))
  parts <- stringr::str_trim(parts)
  parts <- parts[parts != ""]
  parts
}

## ---------------------------------------------------------------------
## 3. برای هر دارو، target ها رو به ژن‌های ماژول map کن
## ---------------------------------------------------------------------

overlay_list <- vector("list", length = nrow(drugs))

for (i in seq_len(nrow(drugs))) {
  drow <- drugs[i, ]
  
  tgt_vec <- split_targets(drow$targets_raw)
  tgt_upper <- toupper(tgt_vec)
  
  ## فقط targetهایی که در ماژول هستند
  tgt_in_module <- unique(tgt_upper[tgt_upper %in% gene_set])
  tgt_in_hubs   <- unique(tgt_in_module[tgt_in_module %in% hub_gene_set])
  
  overlay_list[[i]] <- tibble::tibble(
    rank                 = drow$rank,
    pert_name            = drow$pert_name,
    MOA_category         = drow$MOA_category,
    moa                  = drow$moa,
    n_signatures         = drow$n_signatures,
    min_score            = drow$min_score,
    median_score         = drow$median_score,
    mean_score           = drow$mean_score,
    best_fdr_q_nlog10    = drow$best_fdr_q_nlog10,
    best_cell            = drow$best_cell,
    n_targets_listed     = length(tgt_vec),
    n_targets_in_module  = length(tgt_in_module),
    n_targets_in_hubs    = length(tgt_in_hubs),
    targets_in_module    = paste(tgt_in_module, collapse = ";"),
    targets_in_hubs      = paste(tgt_in_hubs,   collapse = ";")
  )
}

overlay_tbl <- dplyr::bind_rows(overlay_list)

## ---------------------------------------------------------------------
## 4. ذخیره‌ی خروجی‌ها
## ---------------------------------------------------------------------

out_all_path <- file.path(
  drug_dir,
  "CMap_PPI_overlay_cross_species_BROAD_all.tsv"
)

readr::write_tsv(overlay_tbl, out_all_path)
message("[Overlay] Full overlay table written to: ", out_all_path)

## یک نسخه‌ی prioritzed: داروهایی که حداقل ۱ target در ماژول دارند
prioritized_tbl <- overlay_tbl %>%
  dplyr::filter(n_targets_in_module > 0) %>%
  dplyr::arrange(median_score, dplyr::desc(n_targets_in_hubs))

out_prio_path <- file.path(
  drug_dir,
  "CMap_PPI_overlay_cross_species_BROAD_prioritized.tsv"
)

readr::write_tsv(prioritized_tbl, out_prio_path)
message("[Overlay] Prioritized overlay (n_targets_in_module > 0) written to: ",
        out_prio_path)

message("=== [Phase 5B - Overlay CMap hits on PPI hubs] DONE ===")
