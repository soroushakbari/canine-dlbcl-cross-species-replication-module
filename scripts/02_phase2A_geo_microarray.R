## ==========================================
## 02_phase2A_geo_microarray.R
## Phase 2A - Download & process GEO microarrays
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"

raw_geo_dir   <- file.path(project_root, "data", "raw", "GEO")
processed_dir <- file.path(project_root, "data", "processed")
metadata_dir  <- file.path(project_root, "metadata")

dir.create(raw_geo_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir,  recursive = TRUE, showWarnings = FALSE)

inventory_path <- file.path(metadata_dir, "dataset_inventory_phase1.csv")
if (!file.exists(inventory_path)) {
  stop("Inventory file not found: ", inventory_path,
       "\nRun 01_phase1_dataset_inventory.R first.")
}

# پکیج‌ها
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkg_bioc <- c("GEOquery", "Biobase", "SummarizedExperiment", "S4Vectors")
pkg_cran <- c("tidyverse")

for (p in pkg_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

for (p in pkg_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(Biobase)
  library(SummarizedExperiment)
  library(S4Vectors)
})

dataset_inventory <- readr::read_csv(inventory_path, show_col_types = FALSE)

# ---------- توابع کمکی ----------

guess_symbol_column <- function(fdat) {
  cn <- colnames(fdat)
  candidates_exact <- c("Gene.symbol", "Gene Symbol", "GENE_SYMBOL",
                        "SYMBOL", "Symbol", "gene_symbol", "GeneSymbol")
  hit <- candidates_exact[candidates_exact %in% cn][1]
  if (!is.na(hit)) return(hit)
  
  hits <- grep("symbol", cn, ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(hits[1])
  
  stop("No obvious gene symbol column found in featureData. Check platform annotation manually.")
}

collapse_probes_to_genes <- function(expr_mat, gene_symbols,
                                     fun = c("median", "mean")) {
  fun <- match.arg(fun)
  stopifnot(nrow(expr_mat) == length(gene_symbols))
  
  df <- as.data.frame(expr_mat)
  df$gene_symbol <- gene_symbols
  
  df <- df %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "")
  
  df_sum <- df %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        .fns = ~ if (fun == "median") median(., na.rm = TRUE) else mean(., na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  mat <- as.matrix(df_sum[, -1, drop = FALSE])
  rownames(mat) <- df_sum$gene_symbol
  mat
}

compute_qc_metrics <- function(expr_gene) {
  if (!is.matrix(expr_gene)) {
    expr_gene <- as.matrix(expr_gene)
  }
  
  libsize <- colSums(expr_gene, na.rm = TRUE)
  
  qs <- apply(expr_gene, 2, function(x) {
    stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  })
  
  cor_mat <- stats::cor(expr_gene, method = "spearman", use = "pairwise.complete.obs")
  mean_cor <- apply(cor_mat, 2, function(x) mean(x, na.rm = TRUE))
  
  qc <- tibble::tibble(
    sample_id = colnames(expr_gene),
    libsize   = as.numeric(libsize),
    q25       = as.numeric(qs[1, ]),
    median    = as.numeric(qs[2, ]),
    q75       = as.numeric(qs[3, ]),
    mean_cor  = as.numeric(mean_cor)
  )
  
  med_cor <- stats::median(qc$mean_cor, na.rm = TRUE)
  iqr_cor <- stats::IQR(qc$mean_cor, na.rm = TRUE)
  threshold <- med_cor - 3 * iqr_cor
  qc$qc_outlier <- qc$mean_cor < threshold
  
  qc
}

process_geo_microarray <- function(gse_id,
                                   species = c("human", "dog"),
                                   collapse_fun = "median") {
  species <- match.arg(species)
  
  message("\n==============================")
  message("Processing GEO dataset: ", gse_id, " (species: ", species, ")")
  message("==============================")
  
  gse_list <- GEOquery::getGEO(
    GEO = gse_id,
    GSEMatrix = TRUE,
    getGPL = TRUE,
    destdir = raw_geo_dir
  )
  
  if (length(gse_list) > 1) {
    ns <- sapply(gse_list, function(es) ncol(Biobase::exprs(es)))
    idx <- which.max(ns)
    message("Multiple platforms detected; using platform index: ", idx)
    gset <- gse_list[[idx]]
  } else {
    gset <- gse_list[[1]]
  }
  
  platform_id <- Biobase::annotation(gset)
  
  expr <- Biobase::exprs(gset)
  meta <- Biobase::pData(gset)
  fdat <- Biobase::fData(gset)
  
  message("Original expression matrix: ",
          nrow(expr), " probes x ", ncol(expr), " samples.")
  message("Platform: ", platform_id)
  
  qx <- stats::quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  logc <- (qx[5] - qx[1] > 100) || (qx[5] > 50)
  if (logc) {
    message("Data looks non-log; applying log2(x + 1) transformation.")
    expr <- log2(expr + 1)
  } else {
    message("Data already appears to be on log scale.")
  }
  
  sym_col <- guess_symbol_column(fdat)
  message("Using gene symbol column: ", sym_col)
  
  gene_symbols <- as.character(fdat[[sym_col]])
  expr_gene <- collapse_probes_to_genes(expr, gene_symbols, fun = collapse_fun)
  
  message("Collapsed expression matrix: ",
          nrow(expr_gene), " genes x ", ncol(expr_gene), " samples.")
  
  qc <- compute_qc_metrics(expr_gene)
  n_outliers <- sum(qc$qc_outlier, na.rm = TRUE)
  message("QC: flagged ", n_outliers, " sample(s) as potential outliers.")
  
  keep_samples <- qc$sample_id[!qc$qc_outlier]
  expr_gene_qc <- expr_gene[, keep_samples, drop = FALSE]
  
  meta$sample_id <- rownames(meta)
  meta <- meta %>%
    dplyr::left_join(qc, by = c("sample_id" = "sample_id"))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr_gene),
    colData = S4Vectors::DataFrame(meta)
  )
  
  expr_tbl <- expr_gene %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_symbol")
  
  expr_full_path <- file.path(processed_dir, paste0(gse_id, "_expr_log2.tsv"))
  readr::write_tsv(expr_tbl, expr_full_path)
  
  expr_qc_tbl <- expr_gene_qc %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_symbol")
  
  expr_qc_path <- file.path(processed_dir, paste0(gse_id, "_expr_log2_qcfiltered.tsv"))
  readr::write_tsv(expr_qc_tbl, expr_qc_path)
  
  meta_path <- file.path(metadata_dir, paste0(gse_id, "_sample_metadata.csv"))
  readr::write_csv(meta, meta_path)
  
  fdat_tbl <- fdat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("probe_id")
  
  fdat_path <- file.path(metadata_dir, paste0(gse_id, "_feature_annotation.csv"))
  readr::write_csv(fdat_tbl, fdat_path)
  
  se_path <- file.path(processed_dir, paste0(gse_id, "_SummarizedExperiment.rds"))
  saveRDS(se, se_path)
  
  message("Saved:\n  ", expr_full_path,
          "\n  ", expr_qc_path,
          "\n  ", meta_path,
          "\n  ", fdat_path,
          "\n  ", se_path)
  
  if ("dataset_id" %in% colnames(dataset_inventory)) {
    idx <- which(dataset_inventory$dataset_id == gse_id)
    if (length(idx) == 1) {
      dataset_inventory$platform[idx]  <- platform_id
      dataset_inventory$n_samples[idx] <- ncol(expr_gene)
      readr::write_csv(dataset_inventory, inventory_path)
      message("Updated inventory for ", gse_id, " (platform, n_samples).")
    } else {
      warning("Could not uniquely match dataset_id in inventory for: ", gse_id)
    }
  } else {
    warning("Inventory has no 'dataset_id' column. Skipping inventory update.")
  }
  
  invisible(list(
    expr_gene = expr_gene,
    expr_gene_qc = expr_gene_qc,
    meta = meta,
    qc = qc,
    se = se
  ))
}

# ---------- اجرا برای GEO microarrayها ----------

microarray_ids <- c(
  "GSE10846",
  "GSE56315",
  "GSE32018",
  "GSE30881",
  "GSE43664",
  "GSE39365"
)

species_map <- dataset_inventory %>%
  dplyr::filter(dataset_id %in% microarray_ids) %>%
  dplyr::select(dataset_id, species)

default_species <- tibble::tibble(
  dataset_id = microarray_ids,
  species    = "human"
)

if (nrow(species_map) == 0) {
  species_map <- default_species
} else {
  species_map <- default_species %>%
    dplyr::left_join(species_map, by = "dataset_id", suffix = c("_default", "_inv")) %>%
    dplyr::mutate(
      species_final = dplyr::coalesce(species_inv, species_default)
    ) %>%
    dplyr::select(dataset_id, species = species_final)
}

results_list <- list()

for (i in seq_len(nrow(species_map))) {
  ds <- species_map$dataset_id[i]
  sp <- species_map$species[i]
  
  res <- try(
    process_geo_microarray(gse_id = ds, species = sp),
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    message(">>> ERROR while processing ", ds, ":")
    message(as.character(res))
  } else {
    results_list[[ds]] <- res
    message("Finished processing ", ds, ".")
  }
}

cat("\n[Phase 2A] GEO microarrays processed (where possible).\n")
cat("Check 'data/processed' and 'metadata' under:\n", project_root, "\n\n")
