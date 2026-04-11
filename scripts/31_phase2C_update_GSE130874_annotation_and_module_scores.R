## 31_phase2C_update_GSE130874_annotation_and_module_scores.R
## Re-annotate GSE130874 (canine RNA-seq) and compute cross-species BROAD module scores
## + ذخیره‌ی ماتریس gene_symbol برای استفاده در شکل‌ها (Fig 6B)

message("=== Phase 2C update: Re-annotate GSE130874 and compute cross-species BROAD module scores ===")

## ------------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tibble", "stringr")
bioc_pkgs     <- c("rtracklayer")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run:\n",
      "  install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      "))\n",
      "و بعد دوباره اسکریپت را اجرا کن."
    )
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Bioconductor package '", pkg, "' is required but not installed.\n",
      "Install it via:\n",
      "  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")\n",
      "  BiocManager::install(\"", pkg, "\")\n"
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(rtracklayer)
})

## ------------------------------------------------------------------
## 1. Paths
## ------------------------------------------------------------------
project_root <- normalizePath(file.path(".."), winslash = "/", mustWork = TRUE)
message("Project root: ", project_root)

data_dir  <- file.path(project_root, "data", "processed")
meta_dir  <- file.path(project_root, "metadata")
sig_dir   <- file.path(project_root, "results", "tables", "signatures")
score_dir <- file.path(project_root, "results", "tables", "module_scores")

if (!dir.exists(score_dir)) {
  dir.create(score_dir, recursive = TRUE)
  message("Created module_scores dir:\n  ", score_dir)
}

vst_path <- file.path(data_dir, "GSE130874_vst.tsv")
if (!file.exists(vst_path)) {
  stop("VST matrix for GSE130874 not found at:\n  ", vst_path)
}

dog_up_path   <- file.path(sig_dir, "cross_species_module_BROAD_dog_up.txt")
dog_down_path <- file.path(sig_dir, "cross_species_module_BROAD_dog_down.txt")

if (!file.exists(dog_up_path)) {
  stop("Dog BROAD up-signature file not found at:\n  ", dog_up_path)
}
if (!file.exists(dog_down_path)) {
  stop("Dog BROAD down-signature file not found at:\n  ", dog_down_path)
}

## ------------------------------------------------------------------
## 2. Helper: compute module scores (Z-score per gene)
## ------------------------------------------------------------------
compute_module_scores_z <- function(expr_df,
                                    id_col    = "gene_symbol",
                                    up_genes,
                                    down_genes) {
  if (!id_col %in% names(expr_df)) {
    stop("compute_module_scores_z: id_col '", id_col, "' not found in expr_df.")
  }
  
  expr_df <- expr_df %>%
    dplyr::filter(!is.na(.data[[id_col]]), .data[[id_col]] != "")
  
  gene_ids    <- expr_df[[id_col]]
  sample_cols <- setdiff(names(expr_df), id_col)
  
  mat <- expr_df %>%
    dplyr::select(dplyr::all_of(sample_cols)) %>%
    as.matrix()
  
  rownames(mat) <- gene_ids
  
  up_use   <- intersect(up_genes, rownames(mat))
  down_use <- intersect(down_genes, rownames(mat))
  
  message("  Overlap with dog BROAD signature after annotation:")
  message("    up  : ", length(up_use))
  message("    down: ", length(down_use))
  
  if (length(up_use) < 15 || length(down_use) < 15) {
    stop(
      "Overlap between GSE130874 and dog BROAD signature is too small (",
      length(up_use), " up, ", length(down_use), " down).\n",
      "ساختن module score با این تعداد ژن بیشتر نویز است تا بیولوژی."
    )
  }
  
  zmat <- t(scale(t(mat)))
  zmat[is.nan(zmat)] <- NA_real_
  
  score_up   <- colMeans(zmat[up_use, , drop = FALSE], na.rm = TRUE)
  score_down <- colMeans(zmat[down_use, , drop = FALSE], na.rm = TRUE)
  score_net  <- score_up - score_down
  
  tibble::tibble(
    sample_id    = colnames(zmat),
    score_up_z   = as.numeric(score_up),
    score_down_z = as.numeric(score_down),
    score_net_z  = as.numeric(score_net)
  )
}

## ------------------------------------------------------------------
## 3. Step 1 – Load VST matrix
## ------------------------------------------------------------------
message("\n--- Step 1: Load GSE130874 VST matrix ---")

expr_raw <- readr::read_tsv(vst_path, show_col_types = FALSE)

if (!"gene_id" %in% names(expr_raw)) {
  message("  'gene_id' column not found; assuming first column is gene_id.")
  names(expr_raw)[1] <- "gene_id"
}

expr_raw <- expr_raw %>%
  dplyr::mutate(gene_id = as.character(gene_id))

sample_cols <- setdiff(names(expr_raw), "gene_id")

expr_mat <- expr_raw %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

message("  Using 'gene_id' as gene identifier column.")
message("  Number of rows (genes): ", nrow(expr_mat))
message("  Number of samples    : ", length(sample_cols))

## ------------------------------------------------------------------
## 4. Step 2 – Load dog BROAD signature
## ------------------------------------------------------------------
message("\n--- Step 2: Load dog BROAD signature ---")

dog_up <- readr::read_lines(dog_up_path) %>%
  stringr::str_squish() %>%
  .[. != ""]
dog_down <- readr::read_lines(dog_down_path) %>%
  stringr::str_squish() %>%
  .[. != ""]

message("  Dog BROAD up-signature genes   : ", length(dog_up))
message("  Dog BROAD down-signature genes : ", length(dog_down))

## ------------------------------------------------------------------
## 5. Step 3 – Direct overlap check
## ------------------------------------------------------------------
message("\n--- Step 3: Check direct overlap (is gene_id already a symbol?) ---")

sym_vec <- expr_mat$gene_id

overlap_up_direct   <- intersect(sym_vec, dog_up)
overlap_down_direct <- intersect(sym_vec, dog_down)

message("  Direct overlap with dog signature (using 'gene_id' as symbol):")
message("    up  : ", length(overlap_up_direct))
message("    down: ", length(overlap_down_direct))

use_direct_symbols <- (length(overlap_up_direct) >= 50 &&
                         length(overlap_down_direct) >= 50)

if (use_direct_symbols) {
  message("  --> gene_id already looks like dog gene symbols; no re-annotation needed.")
  annot_df <- expr_mat %>%
    dplyr::mutate(gene_symbol = gene_id) %>%
    dplyr::select(gene_symbol, dplyr::all_of(sample_cols))
} else {
  ## --------------------------------------------------------------
  ## 6. Step 4 – Dog10K GFF3 mapping: ENSCAFG → gene_symbol
  ## --------------------------------------------------------------
  message("  --> Direct overlap is low; will attempt Ensembl Dog10K GFF3-based mapping.")
  
  if (!dir.exists(meta_dir)) {
    dir.create(meta_dir, recursive = TRUE)
  }
  
  gff_local <- file.path(
    meta_dir,
    "Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.110.gff3.gz"
  )
  
  gff_url <- paste0(
    "https://ftp.ensembl.org/pub/release-110/gff3/",
    "canis_lupus_familiarisboxer/",
    "Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.110.gff3.gz"
  )
  
  if (!file.exists(gff_local)) {
    message("  Downloading Dog10K GFF3 from Ensembl (one-time, ~20MB) ...")
    ok <- TRUE
    tryCatch(
      {
        utils::download.file(
          url      = gff_url,
          destfile = gff_local,
          mode     = "wb",
          quiet    = TRUE
        )
      },
      error = function(e) {
        ok <<- FALSE
        message("  ERROR during download: ", conditionMessage(e))
      }
    )
    if (!ok || !file.exists(gff_local)) {
      stop(
        "Failed to download Dog10K GFF3 from Ensembl.\n",
        "بدون این فایل، امکان mapping علمی برای GSE130874 وجود ندارد."
      )
    }
  } else {
    message("  Using existing local Dog10K GFF3:\n    ", gff_local)
  }
  
  message("  Importing GFF3 (this may take some time) ...")
  gff <- rtracklayer::import(gff_local)
  
  meta_all <- as.data.frame(mcols(gff))
  if (!"type" %in% names(meta_all)) {
    stop("GFF3 metadata does not contain a 'type' column; cannot select gene rows.")
  }
  
  gene_idx <- which(meta_all$type == "gene")
  if (length(gene_idx) == 0) {
    stop("No 'gene' rows found in Dog10K GFF3; something is wrong with the file.")
  }
  
  gene_meta <- meta_all[gene_idx, , drop = FALSE]
  colnames_meta <- colnames(gene_meta)
  
  id_col <- if ("gene_id" %in% colnames_meta) {
    "gene_id"
  } else if ("ID" %in% colnames_meta) {
    "ID"
  } else {
    NA_character_
  }
  
  sym_col <- if ("gene_name" %in% colnames_meta) {
    "gene_name"
  } else if ("Name" %in% colnames_meta) {
    "Name"
  } else {
    NA_character_
  }
  
  if (is.na(id_col) || is.na(sym_col)) {
    stop(
      "Dog10K GFF3 does not contain both a usable gene ID column and a gene symbol/name column.\n",
      "بدون این اطلاعات، نمی‌توان ENSCAFG را به symbol نگاشت."
    )
  }
  
  mapping <- tibble::tibble(
    ensembl_id  = as.character(gene_meta[[id_col]]),
    gene_symbol = as.character(gene_meta[[sym_col]])
  ) %>%
    dplyr::filter(
      !is.na(ensembl_id), ensembl_id != "",
      !is.na(gene_symbol), gene_symbol != ""
    ) %>%
    dplyr::distinct()
  
  message("  Mapping rows (Dog10K GFF3, ENSCAFG → symbol): ", nrow(mapping))
  
  if (nrow(mapping) == 0) {
    stop(
      "Dog10K GFF3 mapping yielded 0 rows با هر دو فیلد ID و gene_symbol.\n",
      "این یعنی عملاً راهی برای ENSCAFG → symbol نداریم."
    )
  }
  
  expr_ids <- expr_mat %>%
    dplyr::mutate(
      ensembl_id = stringr::str_replace(gene_id, "\\.\\d+$", "")
    ) %>%
    dplyr::select(ensembl_id) %>%
    dplyr::distinct()
  
  expr_ids_annot <- expr_ids %>%
    dplyr::inner_join(mapping, by = "ensembl_id")
  
  message(
    "  Expression genes with Dog10K annotation (after join): ",
    nrow(expr_ids_annot)
  )
  
  if (nrow(expr_ids_annot) == 0) {
    stop(
      "After Dog10K mapping, 0 expression genes had a gene_symbol.\n",
      "عملاً gene_idهای GSE130874 با Dog10K هم سازگار نیستند."
    )
  }
  
  expr_with_map <- expr_mat %>%
    dplyr::mutate(
      ensembl_id = stringr::str_replace(gene_id, "\\.\\d+$", "")
    ) %>%
    dplyr::inner_join(mapping, by = "ensembl_id")
  
  message(
    "  Rows in expression matrix with gene_symbol after join: ",
    nrow(expr_with_map)
  )
  
  if (nrow(expr_with_map) == 0) {
    stop(
      "Inner join with Dog10K mapping produced 0 rows.\n",
      "نمی‌توان بدون annotation درست، module score ساخت."
    )
  }
  
  annot_df <- expr_with_map %>%
    dplyr::select(gene_symbol, dplyr::all_of(sample_cols)) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  message(
    "  Distinct gene symbols after collapsing duplicates: ",
    nrow(annot_df)
  )
}

## ------------------------------------------------------------------
## 6. NEW: save gene-symbol matrix برای استفاده در شکل‌ها
## ------------------------------------------------------------------
expr_symbol_out <- file.path(data_dir, "GSE130874_vst_geneSymbol.tsv")
readr::write_tsv(annot_df, expr_symbol_out)
message("Saved gene-symbol VST matrix to:")
message("  ", expr_symbol_out)

## ------------------------------------------------------------------
## 7. Step 5 – Compute module scores
## ------------------------------------------------------------------
message("\n--- Step 5: Compute cross-species BROAD module scores for GSE130874 ---")

scores_130874 <- compute_module_scores_z(
  expr_df    = annot_df,
  id_col     = "gene_symbol",
  up_genes   = dog_up,
  down_genes = dog_down
)

message("  Number of samples with scores: ", nrow(scores_130874))
message("  Mean score_net_z: ", round(mean(scores_130874$score_net_z, na.rm = TRUE), 3))

## ------------------------------------------------------------------
## 8. Step 6 – Merge with metadata (اگر وجود دارد) و ذخیره
## ------------------------------------------------------------------
meta_path <- file.path(meta_dir, "GSE130874_sample_metadata.csv")
if (file.exists(meta_path)) {
  message("  Merging with sample metadata ...")
  smeta <- readr::read_csv(meta_path, show_col_types = FALSE)
  
  if (!"sample_id" %in% names(smeta)) {
    candidate <- intersect(names(smeta), scores_130874$sample_id)
    if (length(candidate) == 1) {
      smeta <- smeta %>%
        dplyr::rename(sample_id = !!candidate)
    }
  }
  
  if (!"sample_id" %in% names(smeta)) {
    message("  WARNING: metadata has no 'sample_id' column; saving scores بدون merge.")
    scores_final <- scores_130874
  } else {
    scores_final <- scores_130874 %>%
      dplyr::left_join(smeta, by = "sample_id")
  }
} else {
  message("  Metadata file for GSE130874 not found; saving scores بدون merge.")
  scores_final <- scores_130874
}

out_path <- file.path(
  score_dir,
  "module_scores_GSE130874_crossSpeciesBROAD_dog_zscore.tsv"
)

readr::write_tsv(scores_final, out_path)
message("Saved GSE130874 module scores to:")
message("  ", out_path)

message("=== Phase 2C update completed successfully (اگر اینجا رسید یعنی mapping علمی انجام شده) ===")
