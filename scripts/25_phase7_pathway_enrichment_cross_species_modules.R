## 25_phase7_pathway_enrichment_cross_species_modules.R
## هدف:
##  - Pathway enrichment روی cross-species STRICT و BROAD (سمت انسان)
##  - استفاده از KEGG و Reactome (انسان)
##  - خروجی: جدول‌های TSV آماده برای manuscript و تفسیر بیولوژیک

message("=== Phase 7: Pathway enrichment for cross-species STRICT & BROAD modules (human) ===")

## --------------------------------------------------------------------
## 0. پکیج‌ها
## --------------------------------------------------------------------
required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run:\n  install.packages(c(\"readr\",\"dplyr\",\"tibble\",\"stringr\",\"clusterProfiler\"))\n",
      "و همچنین برای بیولوژی:\n  BiocManager::install(c(\"org.Hs.eg.db\", \"ReactomePA\"))\n",
      "سپس اسکریپت را دوباره اجرا کن."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

## --------------------------------------------------------------------
## 1. مسیرها
## --------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)

message("\nProject root: ", project_root)

sig_dir <- file.path(
  project_root, "results", "tables", "signatures"
)

if (!dir.exists(sig_dir)) {
  stop("Signatures directory not found at:\n  ", sig_dir)
}

## BROAD: همون قبلی
path_broad <- file.path(
  sig_dir,
  "cross_species_module_BROAD_table.tsv"
)

if (!file.exists(path_broad)) {
  stop("cross_species_module_BROAD_table.tsv not found at:\n  ", path_broad)
}

## STRICT: اینجا robust می‌شیم و هر دو اسم ممکن رو چک می‌کنیم
strict_candidates <- c(
  "cross_species_module_STRICT_table.tsv",
  "cross_species_module_STRICT_core_table.tsv"
)

path_strict <- NULL
for (fname in strict_candidates) {
  cand <- file.path(sig_dir, fname)
  if (file.exists(cand)) {
    path_strict <- cand
    break
  }
}

if (is.null(path_strict)) {
  stop(
    "Could not find a STRICT module table.\n",
    "Tried the following filenames under ", sig_dir, ":\n  ",
    paste(strict_candidates, collapse = "\n  ")
  )
}

message("Using STRICT module file:\n  ", path_strict)

out_dir <- file.path(
  project_root, "results", "tables", "network"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## --------------------------------------------------------------------
## 2. تابع کمکی: آماده‌سازی لیست ژن انسانی و map به ENTREZ
## --------------------------------------------------------------------
prepare_gene_list <- function(df, module_name = "BROAD") {
  ## حداقل باید یک ستون حاوی سمبل انسانی داشته باشیم
  symbol_cols <- c("human_symbol", "human_gene", "gene_human", "symbol_human")
  symbol_col <- symbol_cols[symbol_cols %in% names(df)][1]
  
  if (is.na(symbol_col)) {
    stop(
      "Could not find a human gene symbol column in cross-species module (",
      module_name, ").\nChecked: ", paste(symbol_cols, collapse = ", ")
    )
  }
  
  message("  [", module_name, "] Using '", symbol_col, "' as human symbol column.")
  
  genes <- df[[symbol_col]] %>%
    as.character() %>%
    toupper() %>%
    trimws()
  
  genes <- unique(genes[genes != "" & !is.na(genes)])
  
  message("  [", module_name, "] Unique human symbols: ", length(genes))
  
  if (length(genes) == 0) {
    stop("[", module_name, "] No valid human gene symbols.")
  }
  
  ## map به ENTREZ
  bitr_res <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  if (is.null(bitr_res) || nrow(bitr_res) == 0) {
    stop(
      "[", module_name, "] bitr mapping returned no ENTREZIDs.\n",
      "Check that org.Hs.eg.db is installed and symbols are human."
    )
  }
  
  bitr_res <- bitr_res %>% distinct(SYMBOL, ENTREZID)
  
  message("  [", module_name, "] Successfully mapped to ENTREZ IDs: ", nrow(bitr_res))
  
  list(
    symbols = genes,
    mapping = bitr_res
  )
}

## --------------------------------------------------------------------
## 3. تابع کمکی برای enrichment و save
## --------------------------------------------------------------------
run_and_save_enrichment <- function(entrez_ids,
                                    module_name = "BROAD",
                                    suffix      = "all") {
  
  if (length(entrez_ids) < 5) {
    message("  [", module_name, "/", suffix, "] Too few genes (", length(entrez_ids),
            ") for meaningful enrichment. Skipping.")
    return(invisible(NULL))
  }
  
  message("  [", module_name, "/", suffix, "] Running KEGG enrichment ...")
  
  ekegg <- tryCatch(
    {
      clusterProfiler::enrichKEGG(
        gene         = entrez_ids,
        organism     = "hsa",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.25
      )
    },
    error = function(e) {
      message("    KEGG enrichment failed: ", e$message)
      NULL
    }
  )
  
  message("  [", module_name, "/", suffix, "] Running Reactome enrichment ...")
  
  ereact <- tryCatch(
    {
      ReactomePA::enrichPathway(
        gene         = entrez_ids,
        organism     = "human",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.25,
        readable     = TRUE
      )
    },
    error = function(e) {
      message("    Reactome enrichment failed: ", e$message)
      NULL
    }
  )
  
  ## ذخیره نتایج در صورت موجود بودن
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    fname_kegg <- file.path(
      out_dir,
      paste0("pathway_enrichment_cross_species_", module_name,
             "_", suffix, "_KEGG.tsv")
    )
    readr::write_tsv(as.data.frame(ekegg), fname_kegg)
    message("    [", module_name, "/", suffix,
            "] Saved KEGG enrichment to:\n      ", fname_kegg)
  } else {
    message("    [", module_name, "/", suffix, "] No significant KEGG terms.")
  }
  
  if (!is.null(ereact) && nrow(as.data.frame(ereact)) > 0) {
    fname_react <- file.path(
      out_dir,
      paste0("pathway_enrichment_cross_species_", module_name,
             "_", suffix, "_Reactome.tsv")
    )
    readr::write_tsv(as.data.frame(ereact), fname_react)
    message("    [", module_name, "/", suffix,
            "] Saved Reactome enrichment to:\n      ", fname_react)
  } else {
    message("    [", module_name, "/", suffix, "] No significant Reactome terms.")
  }
  
  invisible(list(KEGG = ekegg, Reactome = ereact))
}

## --------------------------------------------------------------------
## 4. خواندن cross-species BROAD و STRICT
## --------------------------------------------------------------------
message("\nReading cross-species BROAD module table ...")
broad_tbl <- readr::read_tsv(path_broad, show_col_types = FALSE)

message("Reading cross-species STRICT module table ...")
strict_tbl <- readr::read_tsv(path_strict, show_col_types = FALSE)

## --------------------------------------------------------------------
## 5. آماده‌سازی لیست ژن‌ها و enrichment
## --------------------------------------------------------------------
## BROAD
message("\n=== BROAD module ===")
broad_prep <- prepare_gene_list(broad_tbl, module_name = "BROAD")
broad_entrez_all <- unique(broad_prep$mapping$ENTREZID)

res_broad_all <- run_and_save_enrichment(
  entrez_ids  = broad_entrez_all,
  module_name = "BROAD",
  suffix      = "all"
)

## STRICT
message("\n=== STRICT module ===")
strict_prep <- prepare_gene_list(strict_tbl, module_name = "STRICT")
strict_entrez_all <- unique(strict_prep$mapping$ENTREZID)

res_strict_all <- run_and_save_enrichment(
  entrez_ids  = strict_entrez_all,
  module_name = "STRICT",
  suffix      = "all"
)

message("\n=== Phase 7 pathway enrichment completed. Check KEGG/Reactome TSVs in: ===")
message("  ", out_dir)
