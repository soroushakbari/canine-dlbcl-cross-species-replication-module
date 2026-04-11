## =====================================================================
## Phase 4D: Cross-species disease modules (strict + broad)
## Uses:
##  - Human DE:  DE_GSE56315_tumor_vs_normal.tsv
##  - Dog   DE:  DE_GSE30881_tumor_vs_normal.tsv
##  - Orthologs: fresh biomaRt mapping for ALL human DE genes
## Outputs:
##  - cross_species_module_strict_* (Tier1-style)
##  - cross_species_module_broad_*  (Tier2-style, برای PPI و CMap)
## =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(biomaRt)
})

## ---------------------------------------------------------------------
## 0. مسیر پروژه
## ---------------------------------------------------------------------

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) {
  file.path(project_root, ...)
}

dir.create(path_proj("metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(path_proj("results", "tables", "signatures"),
           recursive = TRUE, showWarnings = FALSE)

sig_dir <- path_proj("results", "tables", "signatures")

## ---------------------------------------------------------------------
## 1. خواندن DE انسان و سگ
## ---------------------------------------------------------------------

load_de_table <- function(path, label) {
  if (!file.exists(path)) {
    stop("[", label, "] DE file not found: ", path)
  }
  df <- read_tsv(path, show_col_types = FALSE)
  
  if (!"gene_symbol" %in% names(df)) {
    gene_col_candidates <- intersect(
      names(df),
      c("GeneSymbol", "symbol", "SYMBOL", "gene", "ID", "probe_id", "Gene_Symbol")
    )
    if (length(gene_col_candidates) == 0L) {
      stop("[", label, "] No 'gene_symbol' column and no obvious alternative.")
    } else {
      df <- df %>% rename(gene_symbol = !!gene_col_candidates[1L])
    }
  }
  if (!"logFC" %in% names(df)) {
    stop("[", label, "] No 'logFC' column in DE table.")
  }
  if (!("adj.P.Val" %in% names(df) | "adj.P.Val." %in% names(df))) {
    stop("[", label, "] No 'adj.P.Val' column in DE table.")
  }
  if ("adj.P.Val." %in% names(df) & !"adj.P.Val" %in% names(df)) {
    df <- df %>% rename(adj.P.Val = adj.P.Val.)
  }
  
  df <- df %>%
    mutate(
      gene_symbol = as.character(gene_symbol)
    )
  
  message("[", label, "] Loaded DE: ", nrow(df), " rows.")
  df
}

de_h_path <- path_proj("results", "tables", "DE", "DE_GSE56315_tumor_vs_normal.tsv")
de_d_path <- path_proj("results", "tables", "DE", "DE_GSE30881_tumor_vs_normal.tsv")

de_h <- load_de_table(de_h_path, "Human (GSE56315)")
de_d <- load_de_table(de_d_path, "Dog   (GSE30881)")

## ---------------------------------------------------------------------
## 2. biomaRt ortholog mapping برای همه‌ی ژن‌های انسانی (DE table)
## ---------------------------------------------------------------------

mapping_path <- path_proj("metadata", "orthologs_human_to_dog_all_DE_GSE56315.csv")

if (file.exists(mapping_path)) {
  orth_all <- read_csv(mapping_path, show_col_types = FALSE)
  message("[Ortholog] Using existing mapping file: ", mapping_path,
          " (", nrow(orth_all), " rows).")
} else {
  human_genes <- unique(de_h$gene_symbol)
  human_genes <- human_genes[!is.na(human_genes) & human_genes != ""]
  
  
  message("[Ortholog] Querying Ensembl via biomaRt for ", length(human_genes),
          " human genes (might take a while).")
  
  ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  attrs <- c(
    "external_gene_name",
    "clfamiliaris_homolog_associated_gene_name",
    "clfamiliaris_homolog_orthology_type"
  )
  
  bm <- biomaRt::getBM(
    attributes = attrs,
    filters    = "external_gene_name",
    values     = human_genes,
    mart       = ensembl
  )
  
  # اگر به هر دلیل اسم ستون‌ها چیزی غیر از attributes شد، ما خودمون ست می‌کنیم
  if (ncol(bm) < 3) {
    stop("[Ortholog] biomaRt returned fewer than 3 columns; attributes changed?")
  }
  
  # attributes به‌ترتیب: external_gene_name, clfamiliaris_homolog_associated_gene_name, clfamiliaris_homolog_orthology_type
  colnames(bm)[1:3] <- c("human_symbol", "dog_symbol", "orthology_type")
  
  orth_all <- bm %>%
    mutate(
      human_symbol = as.character(human_symbol),
      dog_symbol   = as.character(dog_symbol)
    ) %>%
    filter(!is.na(dog_symbol), dog_symbol != "")
  
  orth_all <- orth_all %>%
    distinct(human_symbol, dog_symbol, orthology_type)
  
  write_csv(orth_all, mapping_path)
  message("[Ortholog] Mapping saved to: ", mapping_path,
          " (", nrow(orth_all), " rows).")
}

message("[Ortholog] Mapping table ready: ", nrow(orth_all), " rows; columns: ",
        paste(names(orth_all), collapse = ", "))

## ---------------------------------------------------------------------
## 3. آماده‌سازی DE با threshold و جهت (برای دو لایه: strict و broad)
## ---------------------------------------------------------------------

## thresholds
## strict: بیشتر شبیه Tier1
thr_h_strict_fc  <- 1.5
thr_h_strict_fdr <- 0.01
thr_d_strict_fc  <- 1.0
thr_d_strict_fdr <- 0.05

## broad: شبیه Tier2 + سگ relaxed
thr_h_broad_fc  <- 1.0
thr_h_broad_fdr <- 0.05
thr_d_broad_fc  <- 0.5
thr_d_broad_fdr <- 0.10

annotate_de <- function(df, label,
                        fc_strict, fdr_strict,
                        fc_broad, fdr_broad) {
  df %>%
    mutate(
      is_sig_strict = (abs(logFC) >= fc_strict & adj.P.Val <= fdr_strict),
      dir_strict = case_when(
        is_sig_strict & logFC > 0  ~ "up",
        is_sig_strict & logFC < 0  ~ "down",
        TRUE                      ~ "nonsig"
      ),
      is_sig_broad = (abs(logFC) >= fc_broad & adj.P.Val <= fdr_broad),
      dir_broad = case_when(
        is_sig_broad & logFC > 0  ~ "up",
        is_sig_broad & logFC < 0  ~ "down",
        TRUE                     ~ "nonsig"
      )
    )
}

de_h_annot <- annotate_de(
  de_h, "Human",
  fc_strict = thr_h_strict_fc,  fdr_strict = thr_h_strict_fdr,
  fc_broad  = thr_h_broad_fc,   fdr_broad  = thr_h_broad_fdr
)

de_d_annot <- annotate_de(
  de_d, "Dog",
  fc_strict = thr_d_strict_fc,  fdr_strict = thr_d_strict_fdr,
  fc_broad  = thr_d_broad_fc,   fdr_broad  = thr_d_broad_fdr
)

message("[Human] Strict sig (|logFC|>=", thr_h_strict_fc,
        " & FDR<=", thr_h_strict_fdr, "):")
print(table(direction = de_h_annot$dir_strict[de_h_annot$is_sig_strict]))

message("[Dog] Strict sig (|logFC|>=", thr_d_strict_fc,
        " & FDR<=", thr_d_strict_fdr, "):")
print(table(direction = de_d_annot$dir_strict[de_d_annot$is_sig_strict]))

message("[Human] Broad sig (|logFC|>=", thr_h_broad_fc,
        " & FDR<=", thr_h_broad_fdr, "):")
print(table(direction = de_h_annot$dir_broad[de_h_annot$is_sig_broad]))

message("[Dog] Broad sig (|logFC|>=", thr_d_broad_fc,
        " & FDR<=", thr_d_broad_fdr, "):")
print(table(direction = de_d_annot$dir_broad[de_d_annot$is_sig_broad]))

## ---------------------------------------------------------------------
## 4. Merge orthologs با DEها
## ---------------------------------------------------------------------

de_h_sel <- de_h_annot %>%
  dplyr::select(
    human_symbol = gene_symbol,
    human_logFC  = logFC,
    human_adj.P  = adj.P.Val,
    human_is_sig_strict = is_sig_strict,
    human_dir_strict    = dir_strict,
    human_is_sig_broad  = is_sig_broad,
    human_dir_broad     = dir_broad
  )


de_d_sel <- de_d_annot %>%
  dplyr::select(
    dog_symbol = gene_symbol,
    dog_logFC  = logFC,
    dog_adj.P  = adj.P.Val,
    dog_is_sig_strict = is_sig_strict,
    dog_dir_strict    = dir_strict,
    dog_is_sig_broad  = is_sig_broad,
    dog_dir_broad     = dir_broad
  )


cross_all <- orth_all %>%
  left_join(de_h_sel, by = "human_symbol") %>%
  left_join(de_d_sel, by = "dog_symbol")

message("[Cross] After merge: ", nrow(cross_all), " ortholog pairs.")

## ---------------------------------------------------------------------
## 5. تعریف cross-species strict core (Tier1-style)
## ---------------------------------------------------------------------

cross_strict <- cross_all %>%
  filter(
    !is.na(human_is_sig_strict),
    !is.na(dog_is_sig_strict),
    human_is_sig_strict,
    dog_is_sig_strict,
    human_dir_strict %in% c("up", "down"),
    dog_dir_strict   %in% c("up", "down"),
    human_dir_strict == dog_dir_strict
  )

n_total <- nrow(cross_all)
n_strict <- nrow(cross_strict)

message("[Cross-STRICT] Core size: ", n_strict, " pairs (from ", n_total, " ortholog pairs).")
tab_strict <- table(direction = cross_strict$human_dir_strict)
message("[Cross-STRICT] Direction breakdown:")
print(tab_strict)

strict_up   <- cross_strict %>% filter(human_dir_strict == "up")
strict_down <- cross_strict %>% filter(human_dir_strict == "down")

strict_human_up   <- strict_up$human_symbol   %>% unique() %>% sort()
strict_human_down <- strict_down$human_symbol %>% unique() %>% sort()
strict_dog_up     <- strict_up$dog_symbol     %>% unique() %>% sort()
strict_dog_down   <- strict_down$dog_symbol   %>% unique() %>% sort()

message("[Cross-STRICT] Human up:   ", length(strict_human_up))
message("[Cross-STRICT] Human down: ", length(strict_human_down))
message("[Cross-STRICT] Dog up:     ", length(strict_dog_up))
message("[Cross-STRICT] Dog down:   ", length(strict_dog_down))

## ذخیره‌ی strict
write_tsv(cross_strict,
          file.path(sig_dir, "cross_species_module_STRICT_core_table.tsv"))
write_lines(strict_human_up,
            file.path(sig_dir, "cross_species_module_STRICT_human_up.txt"))
write_lines(strict_human_down,
            file.path(sig_dir, "cross_species_module_STRICT_human_down.txt"))
write_lines(strict_dog_up,
            file.path(sig_dir, "cross_species_module_STRICT_dog_up.txt"))
write_lines(strict_dog_down,
            file.path(sig_dir, "cross_species_module_STRICT_dog_down.txt"))

## ---------------------------------------------------------------------
## 6. تعریف cross-species broad module (Tier2-style)
## ---------------------------------------------------------------------

cross_broad <- cross_all %>%
  filter(
    !is.na(human_is_sig_broad),
    !is.na(dog_is_sig_broad),
    human_is_sig_broad,
    dog_is_sig_broad,
    human_dir_broad %in% c("up", "down"),
    dog_dir_broad   %in% c("up", "down"),
    human_dir_broad == dog_dir_broad
  )

n_broad <- nrow(cross_broad)

message("[Cross-BROAD] Module size: ", n_broad, " pairs (from ", n_total, " ortholog pairs).")
tab_broad <- table(direction = cross_broad$human_dir_broad)
message("[Cross-BROAD] Direction breakdown:")
print(tab_broad)

broad_up   <- cross_broad %>% filter(human_dir_broad == "up")
broad_down <- cross_broad %>% filter(human_dir_broad == "down")

broad_human_up   <- broad_up$human_symbol   %>% unique() %>% sort()
broad_human_down <- broad_down$human_symbol %>% unique() %>% sort()
broad_dog_up     <- broad_up$dog_symbol     %>% unique() %>% sort()
broad_dog_down   <- broad_down$dog_symbol   %>% unique() %>% sort()

message("[Cross-BROAD] Human up:   ", length(broad_human_up))
message("[Cross-BROAD] Human down: ", length(broad_human_down))
message("[Cross-BROAD] Dog up:     ", length(broad_dog_up))
message("[Cross-BROAD] Dog down:   ", length(broad_dog_down))

## ذخیره‌ی broad
write_tsv(cross_broad,
          file.path(sig_dir, "cross_species_module_BROAD_table.tsv"))
write_lines(broad_human_up,
            file.path(sig_dir, "cross_species_module_BROAD_human_up.txt"))
write_lines(broad_human_down,
            file.path(sig_dir, "cross_species_module_BROAD_human_down.txt"))
write_lines(broad_dog_up,
            file.path(sig_dir, "cross_species_module_BROAD_dog_up.txt"))
write_lines(broad_dog_down,
            file.path(sig_dir, "cross_species_module_BROAD_dog_down.txt"))

## ---------------------------------------------------------------------
## 7. خلاصه‌ی متنی
## ---------------------------------------------------------------------

summary_path <- file.path(sig_dir, "cross_species_module_STRICT_and_BROAD_summary.txt")
sink(summary_path)
cat("Cross-species disease modules (STRICT + BROAD)\n\n")
cat("Total ortholog pairs in mapping: ", n_total, "\n\n")

cat("STRICT module (Tier1-style):\n")
cat("  Thresholds (human): |logFC| >=", thr_h_strict_fc,
    ", FDR <=", thr_h_strict_fdr, "\n")
cat("  Thresholds (dog):   |logFC| >=", thr_d_strict_fc,
    ", FDR <=", thr_d_strict_fdr, "\n")
cat("  Pairs: ", n_strict, "\n")
cat("  Direction breakdown:\n")
print(tab_strict)
cat("  Human up:   ", length(strict_human_up), "\n")
cat("  Human down: ", length(strict_human_down), "\n")
cat("  Dog up:     ", length(strict_dog_up), "\n")
cat("  Dog down:   ", length(strict_dog_down), "\n\n")

cat("BROAD module (Tier2-style):\n")
cat("  Thresholds (human): |logFC| >=", thr_h_broad_fc,
    ", FDR <=", thr_h_broad_fdr, "\n")
cat("  Thresholds (dog):   |logFC| >=", thr_d_broad_fc,
    ", FDR <=", thr_d_broad_fdr, "\n")
cat("  Pairs: ", n_broad, "\n")
cat("  Direction breakdown:\n")
print(tab_broad)
cat("  Human up:   ", length(broad_human_up), "\n")
cat("  Human down: ", length(broad_human_down), "\n")
cat("  Dog up:     ", length(broad_dog_up), "\n")
cat("  Dog down:   ", length(broad_dog_down), "\n")
sink()

message("[Summary] Written to: ", summary_path)
message("=== [Phase 4D - Cross-species modules (STRICT + BROAD)] DONE ===")
