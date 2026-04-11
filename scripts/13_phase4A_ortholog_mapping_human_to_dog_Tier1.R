## ==========================================
## 13_phase4A_ortholog_mapping_human_to_dog_Tier1.R
## Phase 4A - Map Human Tier1 core DLBCL signature to Dog orthologs
## Strategy:
##   1) Try Ensembl biomaRt (human -> clfamiliaris homologs via getBM)
##   2) If biomaRt fails or returns zero rows, fallback to shared-symbol
##      mapping using canine feature annotations (GSE39365 & GSE130874)
##
## Inputs:
##   results/tables/signatures/human_signature_GSE56315_Tier1_core_up.txt
##   results/tables/signatures/human_signature_GSE56315_Tier1_core_down.txt
##   metadata/GSE39365_feature_annotation.csv   (if exists)
##   metadata/GSE130874_feature_annotation.csv  (if exists)
##
## Outputs:
##   metadata/orthologs_human_to_dog_Tier1_GSE56315.csv
##   results/tables/signatures/dog_signature_Tier1_from_human_up_initial.txt
##   results/tables/signatures/dog_signature_Tier1_from_human_down_initial.txt
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"

sig_dir      <- file.path(project_root, "results", "tables", "signatures")
metadata_dir <- file.path(project_root, "metadata")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

## ---------- 0) Packages ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt", ask = FALSE, update = FALSE)
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
})

## ---------- 1) Load Human Tier1 core (up/down lists) ----------

tier1_up_file   <- file.path(sig_dir, "human_signature_GSE56315_Tier1_core_up.txt")
tier1_down_file <- file.path(sig_dir, "human_signature_GSE56315_Tier1_core_down.txt")

if (!file.exists(tier1_up_file) || !file.exists(tier1_down_file)) {
  stop("Tier1 up/down files not found.\nExpected:\n  ",
       tier1_up_file, "\n  ", tier1_down_file,
       "\nRun 07_phase3C_refine_human_signature_multi_tier.R first.")
}

genes_up_h   <- readr::read_lines(tier1_up_file)
genes_down_h <- readr::read_lines(tier1_down_file)

genes_up_h   <- unique(genes_up_h[genes_up_h != ""])
genes_down_h <- unique(genes_down_h[genes_down_h != ""])

message("Human Tier1 core: ", length(genes_up_h), " up, ",
        length(genes_down_h), " down.")

if (length(genes_up_h) == 0 || length(genes_down_h) == 0) {
  stop("Tier1 core lists are empty; check Tier1 files.")
}

tier1_human_df <- tibble::tibble(
  human_symbol = c(genes_up_h, genes_down_h),
  direction    = c(rep("up",   length(genes_up_h)),
                   rep("down", length(genes_down_h)))
)

all_human_symbols <- unique(tier1_human_df$human_symbol)

## ---------- 2) Helper: try biomaRt ortholog mapping (human -> dog) ----------

get_orthologs_biomart <- function(human_symbols) {
  message("\n[biomaRt] Connecting to Ensembl hsapiens_gene_ensembl ...")
  mart_hs <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl"
  )
  
  message("[biomaRt] Querying clfamiliaris homologs via getBM ...")
  
  # attributes for dog ortholog; prefix = clfamiliaris (dataset clfamiliaris_gene_ensembl)
  attrs <- c("external_gene_name",
             "clfamiliaris_homolog_associated_gene_name",
             "clfamiliaris_homolog_orthology_type")
  
  res <- biomaRt::getBM(
    attributes = attrs,
    filters    = "external_gene_name",
    values     = human_symbols,
    mart       = mart_hs
  )
  
  if (nrow(res) == 0) {
    message("[biomaRt] getBM returned 0 rows.")
    return(NULL)
  }
  
  colnames(res) <- c("human_symbol", "dog_symbol", "orthology_type")
  
  res <- res %>%
    dplyr::filter(!is.na(dog_symbol), dog_symbol != "") %>%
    dplyr::distinct()
  
  if (nrow(res) == 0) {
    message("[biomaRt] After filtering empty dog_symbol, 0 rows remain.")
    return(NULL)
  }
  
  message("[biomaRt] Raw biomaRt orthologs (human->dog): ", nrow(res))
  
  res
}

orth_bm <- NULL

orth_bm <- tryCatch(
  get_orthologs_biomart(all_human_symbols),
  error = function(e) {
    message("[biomaRt] ERROR while querying Ensembl: ", conditionMessage(e))
    message("[biomaRt] Will fall back to shared-symbol mapping using canine annotations.")
    NULL
  }
)

## ---------- 3) Helper: fallback shared-symbol mapping ----------

get_dog_symbols_from_annotations <- function() {
  
  dog_symbols <- character(0)
  
  annot_files <- c(
    file.path(metadata_dir, "GSE39365_feature_annotation.csv"),
    file.path(metadata_dir, "GSE130874_feature_annotation.csv")
  )
  
  for (f in annot_files) {
    if (!file.exists(f)) {
      message("[Fallback] Annotation file not found (skipping): ", f)
      next
    }
    
    message("[Fallback] Reading canine annotation: ", f)
    annot <- readr::read_csv(f, show_col_types = FALSE)
    
    cn <- colnames(annot)
    symbol_cols <- grep("symbol|gene_name", cn, ignore.case = TRUE, value = TRUE)
    
    if (length(symbol_cols) == 0) {
      message("[Fallback] No symbol/gene_name column detected in: ", f)
      next
    }
    
    sym_col <- symbol_cols[1]
    message("[Fallback] Using column '", sym_col,
            "' as gene symbol in ", basename(f))
    
    syms <- annot[[sym_col]]
    syms <- unique(syms[!is.na(syms) & syms != ""])
    
    dog_symbols <- union(dog_symbols, syms)
  }
  
  dog_symbols
}

## ---------- 4) Build final ortholog table (biomaRt OR fallback) ----------

if (!is.null(orth_bm)) {
  # biomaRt خودش orthology_type می‌دهد؛ ما فقط Tier1 را روی آن ساب‌ست می‌کنیم.
  message("\nUsing biomaRt orthologs for Tier1 mapping ...")
  
  orth_df <- orth_bm %>%
    dplyr::inner_join(tier1_human_df, by = "human_symbol") %>%
    dplyr::distinct()
  
  n_human_mapped <- length(unique(orth_df$human_symbol))
  n_dog_genes    <- length(unique(orth_df$dog_symbol))
  
  message("Ortholog pairs that are in Tier1 human signature: ", nrow(orth_df))
  message("Unique human Tier1 genes with a dog ortholog (biomaRt): ", n_human_mapped)
  message("Unique dog genes mapped from Tier1 (biomaRt): ", n_dog_genes)
  
  mapping_source <- "biomaRt_hsapiens_clfamiliaris"
} else {
  message("\n[Fallback] Building orthologs by shared symbols between human Tier1 and canine annotations ...")
  
  dog_symbols <- get_dog_symbols_from_annotations()
  
  if (length(dog_symbols) == 0) {
    stop("[Fallback] Could not obtain any dog gene symbols from annotations. ",
         "Check canine feature annotation files.")
  }
  
  shared_syms <- intersect(all_human_symbols, dog_symbols)
  
  if (length(shared_syms) == 0) {
    stop("[Fallback] No shared gene symbols between human Tier1 and canine annotations.")
  }
  
  orth_df <- tibble::tibble(
    human_symbol = shared_syms
  ) %>%
    dplyr::left_join(tier1_human_df, by = "human_symbol") %>%
    dplyr::mutate(
      dog_symbol    = human_symbol,
      orthology_type = "shared_symbol"
    ) %>%
    dplyr::distinct()
  
  n_human_mapped <- length(unique(orth_df$human_symbol))
  n_dog_genes    <- length(unique(orth_df$dog_symbol))
  
  message("[Fallback] Shared-symbol ortholog pairs (human->dog): ", nrow(orth_df))
  message("[Fallback] Unique human Tier1 genes with shared-symbol match: ", n_human_mapped)
  message("[Fallback] Unique dog genes mapped from Tier1 (shared-symbol): ", n_dog_genes)
  
  mapping_source <- "shared_symbol_via_canine_annotations"
}

## ---------- 5) Save mapping table ----------

orth_tier1 <- orth_df %>%
  dplyr::select(human_symbol, dog_symbol, direction, orthology_type) %>%
  dplyr::arrange(direction, human_symbol, dog_symbol)

out_csv <- file.path(metadata_dir, "orthologs_human_to_dog_Tier1_GSE56315.csv")
readr::write_csv(orth_tier1, out_csv)

message("\n[Phase 4A] Ortholog mapping (Human Tier1 -> Dog) completed.")
message("Mapping source: ", mapping_source)
message("Mapping table saved to:\n  ", out_csv)

## ---------- 6) Build initial Dog up/down lists (projected Tier1) ----------

dog_up <- orth_tier1 %>%
  dplyr::filter(direction == "up") %>%
  dplyr::pull(dog_symbol) %>%
  unique() %>%
  sort()

dog_down <- orth_tier1 %>%
  dplyr::filter(direction == "down") %>%
  dplyr::pull(dog_symbol) %>%
  unique() %>%
  sort()

message("\nDog Tier1-projected genes: ", length(dog_up), " up, ",
        length(dog_down), " down.")

dog_up_file   <- file.path(sig_dir, "dog_signature_Tier1_from_human_up_initial.txt")
dog_down_file <- file.path(sig_dir, "dog_signature_Tier1_from_human_down_initial.txt")

readr::write_lines(dog_up,   dog_up_file)
readr::write_lines(dog_down, dog_down_file)

message("Initial dog up list saved to:\n  ", dog_up_file)
message("Initial dog down list saved to:\n  ", dog_down_file)
