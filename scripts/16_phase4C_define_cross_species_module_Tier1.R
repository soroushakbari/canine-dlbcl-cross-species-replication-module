## =====================================================================
## Phase 4C: Define cross-species disease module (Tier1 core)
## Using:
##  - Human DE:  DE_GSE56315_tumor_vs_normal.tsv
##  - Canine DE: DE_GSE30881_tumor_vs_normal.tsv
##  - Ortholog map: orthologs_human_to_dog_Tier1_GSE56315.csv
## =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})

## ---------------------------------------------------------------------
## 0. Project paths
## ---------------------------------------------------------------------

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) {
  file.path(project_root, ...)
}

dir.create(path_proj("results", "tables", "signatures"),
           recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1. Load DE tables (human + dog)
## ---------------------------------------------------------------------

load_de_table <- function(path, label) {
  if (!file.exists(path)) {
    stop("[", label, "] DE file not found: ", path)
  }
  
  df <- readr::read_tsv(path, show_col_types = FALSE)
  
  # انتظار: ستون gene_symbol، logFC، adj.P.Val
  if (!"gene_symbol" %in% names(df)) {
    # تلاش برای کشف ستون gene
    gene_col_candidates <- intersect(
      names(df),
      c("GeneSymbol", "symbol", "SYMBOL", "gene", "ID", "probe_id", "Gene_Symbol")
    )
    if (length(gene_col_candidates) == 0L) {
      stop("[", label, "] No 'gene_symbol' column and no obvious alternative in DE table.")
    } else {
      gene_col <- gene_col_candidates[1L]
      df <- df %>%
        rename(gene_symbol = !!gene_col)
    }
  }
  
  if (!"logFC" %in% names(df)) {
    stop("[", label, "] DE table lacks 'logFC' column.")
  }
  if (!("adj.P.Val" %in% names(df) | "adj.P.Val." %in% names(df))) {
    stop("[", label, "] DE table lacks 'adj.P.Val' column.")
  }
  
  if ("adj.P.Val." %in% names(df) & !"adj.P.Val" %in% names(df)) {
    df <- df %>% rename(adj.P.Val = adj.P.Val.)
  }
  
  df <- df %>%
    mutate(
      gene_symbol = as.character(gene_symbol),
      is_sig = (abs(logFC) >= 1 & adj.P.Val <= 0.05),
      direction = case_when(
        is_sig & logFC > 0  ~ "up",
        is_sig & logFC < 0  ~ "down",
        TRUE                ~ "nonsig"
      )
    )
  
  message("[", label, "] Loaded DE table: ", nrow(df), " rows.")
  sig_tab <- table(direction = df$direction[df$is_sig])
  message("[", label, "] Significant (|log2FC|>=1 & adj.P.Val<=0.05):")
  print(sig_tab)
  
  df
}

de_h_path <- path_proj("results", "tables", "DE", "DE_GSE56315_tumor_vs_normal.tsv")
de_d_path <- path_proj("results", "tables", "DE", "DE_GSE30881_tumor_vs_normal.tsv")

de_h <- load_de_table(de_h_path, "Human (GSE56315)")
de_d <- load_de_table(de_d_path, "Dog   (GSE30881)")

## ---------------------------------------------------------------------
## 2. Load ortholog mapping (Tier1 human -> dog)
## ---------------------------------------------------------------------

orth_path <- path_proj("metadata", "orthologs_human_to_dog_Tier1_GSE56315.csv")
if (!file.exists(orth_path)) {
  stop("[Ortholog] Mapping file not found: ", orth_path)
}

orth <- readr::read_csv(orth_path, show_col_types = FALSE)
message("[Ortholog] Loaded mapping table: ", nrow(orth), " rows; columns: ",
        paste(names(orth), collapse = ", "))

## در این پروژه، فایل mapping همین حالا ستون‌های human_symbol و dog_symbol را دارد.
if (!all(c("human_symbol", "dog_symbol") %in% names(orth))) {
  stop("[Ortholog] Expected columns 'human_symbol' and 'dog_symbol' not found in mapping file.")
}

orth2 <- orth %>%
  mutate(
    human_symbol = as.character(human_symbol),
    dog_symbol   = as.character(dog_symbol)
  )


## ---------------------------------------------------------------------
## 3. Merge orthologs with DE (human + dog)
## ---------------------------------------------------------------------

de_h_sel <- de_h %>%
  select(
    human_symbol = gene_symbol,
    human_logFC  = logFC,
    human_adj.P  = adj.P.Val,
    human_is_sig = is_sig,
    human_dir    = direction
  )

de_d_sel <- de_d %>%
  select(
    dog_symbol = gene_symbol,
    dog_logFC  = logFC,
    dog_adj.P  = adj.P.Val,
    dog_is_sig = is_sig,
    dog_dir    = direction
  )

cross_raw <- orth2 %>%
  left_join(de_h_sel, by = "human_symbol") %>%
  left_join(de_d_sel, by = "dog_symbol")

message("[Cross] After merge: ", nrow(cross_raw), " ortholog pairs.")

## ---------------------------------------------------------------------
## 4. Define cross-species disease module (Tier1 core)
##    Criteria:
##      - human_is_sig == TRUE
##      - dog_is_sig   == TRUE
##      - human_dir, dog_dir ∈ {up, down}
##      - human_dir == dog_dir  (هم‌جهت: up/up یا down/down)
## ---------------------------------------------------------------------

cross_core <- cross_raw %>%
  filter(
    !is.na(human_is_sig),
    !is.na(dog_is_sig),
    human_is_sig,
    dog_is_sig,
    human_dir %in% c("up", "down"),
    dog_dir   %in% c("up", "down"),
    human_dir == dog_dir
  )

n_total <- nrow(cross_raw)
n_core  <- nrow(cross_core)

message("[Cross] Cross-species core module size: ", n_core,
        " pairs (from ", n_total, " ortholog pairs).")

tab_core <- table(direction = cross_core$human_dir)
message("[Cross] Direction breakdown in core (by human_dir = dog_dir):")
print(tab_core)

## ---------------------------------------------------------------------
## 5. Split into up / down and prepare gene lists
## ---------------------------------------------------------------------

cross_up <- cross_core %>% filter(human_dir == "up")
cross_down <- cross_core %>% filter(human_dir == "down")

## Human gene lists
human_up_genes   <- cross_up$human_symbol   %>% unique() %>% sort()
human_down_genes <- cross_down$human_symbol %>% unique() %>% sort()

## Dog gene lists
dog_up_genes   <- cross_up$dog_symbol   %>% unique() %>% sort()
dog_down_genes <- cross_down$dog_symbol %>% unique() %>% sort()

message("[Cross] Human up genes:   ", length(human_up_genes))
message("[Cross] Human down genes: ", length(human_down_genes))
message("[Cross] Dog up genes:     ", length(dog_up_genes))
message("[Cross] Dog down genes:   ", length(dog_down_genes))

## ---------------------------------------------------------------------
## 6. Write outputs
## ---------------------------------------------------------------------

sig_dir <- path_proj("results", "tables", "signatures")

## جدول کامل core module
core_table_path <- file.path(sig_dir, "cross_species_module_Tier1_core_table.tsv")
readr::write_tsv(cross_core, core_table_path)
message("[Cross] Core table written to: ", core_table_path)

## لیست‌ها – انسان
human_up_path   <- file.path(sig_dir, "cross_species_module_Tier1_core_human_up.txt")
human_down_path <- file.path(sig_dir, "cross_species_module_Tier1_core_human_down.txt")
readr::write_lines(human_up_genes,   human_up_path)
readr::write_lines(human_down_genes, human_down_path)
message("[Cross] Human up list   -> ", human_up_path)
message("[Cross] Human down list -> ", human_down_path)

## لیست‌ها – سگ
dog_up_path   <- file.path(sig_dir, "cross_species_module_Tier1_core_dog_up.txt")
dog_down_path <- file.path(sig_dir, "cross_species_module_Tier1_core_dog_down.txt")
readr::write_lines(dog_up_genes,   dog_up_path)
readr::write_lines(dog_down_genes, dog_down_path)
message("[Cross] Dog up list   -> ", dog_up_path)
message("[Cross] Dog down list -> ", dog_down_path)

## خلاصه‌ی متنی برای رجوع سریع
summary_path <- file.path(sig_dir, "cross_species_module_Tier1_core_summary.txt")
sink(summary_path)
cat("Cross-species disease module (Tier1 core)\n\n")
cat("Total ortholog pairs in mapping: ", n_total, "\n")
cat("Core ortholog pairs (sig human & dog, same direction): ", n_core, "\n\n")
cat("Direction breakdown (human_dir = dog_dir):\n")
print(tab_core)
cat("\nHuman up genes:   ", length(human_up_genes), "\n")
cat("Human down genes: ", length(human_down_genes), "\n")
cat("Dog up genes:     ", length(dog_up_genes), "\n")
cat("Dog down genes:   ", length(dog_down_genes), "\n")
sink()

message("[Cross] Summary written to: ", summary_path)
message("=== [Phase 4C - Cross-species module (Tier1 core)] DONE ===")
