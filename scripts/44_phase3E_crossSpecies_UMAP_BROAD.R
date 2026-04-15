
## 44_phase3E_crossSpecies_UMAP_BROAD.R
## Cross-species UMAP using human–dog BROAD module (GSE56315 & GSE30881)
## Revised for reviewer response:
##  - fixed cohort annotation (human controls + canine controls)
##  - explicit within-cohort gene-wise z-scoring
##  - cleaned duplicate code
##  - improved legend for Figure 1A

message("=== Phase 3E: Cross-species UMAP using BROAD module (GSE56315 + GSE30881) ===")

required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "ggplot2", "uwot"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run:\n",
      "  install.packages(c(",
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
  library(ggplot2)
  library(uwot)
})

## ---------------------------------------------------------------------
## paths
## ---------------------------------------------------------------------
project_root <- normalizePath(
  "D:/Research/My Articles/DLBCL drug",
  winslash = "/",
  mustWork = TRUE
)
message("Project root: ", project_root)

data_dir <- file.path(project_root, "data", "processed")
meta_dir <- file.path(project_root, "metadata")
sig_dir  <- file.path(project_root, "results", "tables", "signatures")
fig_dir  <- file.path(project_root, "results", "figures")
aux_dir  <- file.path(project_root, "results", "tables", "network")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(aux_dir)) dir.create(aux_dir, recursive = TRUE)

path_cs_broad   <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
path_expr_human <- file.path(data_dir, "GSE56315_expr_log2_qcfiltered.tsv")
path_meta_human <- file.path(meta_dir, "GSE56315_sample_metadata.csv")
path_expr_dog   <- file.path(data_dir, "GSE30881_expr_log2_qcfiltered.tsv")
path_meta_dog   <- file.path(meta_dir, "GSE30881_sample_metadata.csv")

for (p in c(path_cs_broad, path_expr_human, path_meta_human,
            path_expr_dog, path_meta_dog)) {
  if (!file.exists(p)) stop("Required file not found:\n  ", p)
}

## ---------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------
load_expr_gene_symbol <- function(path_expr, label) {
  message("Loading expression for ", label, " from:\n  ", path_expr)
  
  df <- readr::read_tsv(path_expr, show_col_types = FALSE)
  
  gene_col <- intersect(
    c("gene_symbol", "Gene symbol", "Gene_symbol", "GeneSymbol", "symbol"),
    names(df)
  )[1]
  
  if (is.na(gene_col)) {
    stop("No gene_symbol-like column found in expression table for ", label)
  }
  
  g <- df[[gene_col]]
  
  if (anyDuplicated(g) > 0) {
    message("  [", label, "] collapsing duplicated gene symbols by mean")
    df <- df %>%
      group_by(.data[[gene_col]]) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
    g <- df[[gene_col]]
  }
  
  mat <- df %>%
    select(-all_of(gene_col)) %>%
    as.matrix()
  
  rownames(mat) <- g
  message("  [", label, "] Genes: ", nrow(mat), " ; Samples: ", ncol(mat))
  mat
}

collapse_text_fields <- function(df, cols) {
  cols_use <- cols[cols %in% names(df)]
  if (length(cols_use) == 0) return(rep("", nrow(df)))
  
  out <- apply(df[, cols_use, drop = FALSE], 1, function(x) {
    paste(as.character(x), collapse = " | ")
  })
  out[is.na(out)] <- ""
  out
}

infer_condition_human <- function(meta_human) {
  txt <- collapse_text_fields(
    meta_human,
    c(
      "title",
      "source_name_ch1",
      "tissue:ch1",
      "cell type:ch1",
      "description",
      "abc/gcb/nc subclass:ch1"
    )
  )
  txt <- stringr::str_to_lower(txt)
  
  dplyr::case_when(
    str_detect(txt, "dlbcl|lymphoma") ~ "Tumor",
    str_detect(txt, "\\babc\\b|\\bgcb\\b") ~ "Tumor",
    str_detect(txt, "\\bnc\\b") ~ "Normal",
    str_detect(txt, "normal|reactive|benign|control|tonsil|hyperplasia|non.?malignant") ~ "Normal",
    TRUE ~ "Unknown"
  )
}

infer_condition_dog <- function(meta_dog) {
  txt <- collapse_text_fields(
    meta_dog,
    c("title", "source_name_ch1", "disease state:ch1", "characteristics_ch1")
  )
  txt <- stringr::str_to_lower(txt)
  
  dplyr::case_when(
    "group_tn" %in% names(meta_dog) & meta_dog$group_tn == "tumor" ~ "Tumor",
    "group_tn" %in% names(meta_dog) & is.na(meta_dog$group_tn) ~ "Normal",
    str_detect(txt, "normal|control|healthy") ~ "Normal",
    str_detect(txt, "dlbcl|lymphoma|tumou?r") ~ "Tumor",
    TRUE ~ "Unknown"
  )
}

zscore_by_gene <- function(mat) {
  z <- t(scale(t(mat), center = TRUE, scale = TRUE))
  z[is.na(z)] <- 0
  z
}

## ---------------------------------------------------------------------
## 1. cross-species BROAD table
## ---------------------------------------------------------------------
message("--- Step 1: Load cross-species BROAD module table ---")
cs_broad <- readr::read_tsv(path_cs_broad, show_col_types = FALSE) %>%
  filter(!is.na(human_symbol), !is.na(dog_symbol))
message("  cross-species BROAD pairs: ", nrow(cs_broad))

## ---------------------------------------------------------------------
## 2. load human and dog expression + metadata
## ---------------------------------------------------------------------
message("--- Step 2: Load GSE56315 (human) ---")
expr_human <- load_expr_gene_symbol(path_expr_human, "GSE56315 (human)")
meta_human <- readr::read_csv(path_meta_human, show_col_types = FALSE)

if (!"sample_id" %in% names(meta_human)) {
  stop("GSE56315 metadata must contain 'sample_id'.")
}
meta_human <- meta_human %>%
  mutate(
    sample_id = as.character(sample_id),
    condition = infer_condition_human(.)
  )

message("Human condition counts:")
print(table(meta_human$condition, useNA = "ifany"))

message("--- Step 3: Load GSE30881 (dog) ---")
expr_dog <- load_expr_gene_symbol(path_expr_dog, "GSE30881 (dog)")
meta_dog <- readr::read_csv(path_meta_dog, show_col_types = FALSE)

if (!"sample_id" %in% names(meta_dog)) {
  stop("GSE30881 metadata must contain 'sample_id'.")
}
meta_dog <- meta_dog %>%
  mutate(
    sample_id = as.character(sample_id),
    condition = infer_condition_dog(.)
  )

message("Dog condition counts:")
print(table(meta_dog$condition, useNA = "ifany"))

## ---------------------------------------------------------------------
## 4. build common human_symbol gene space
## ---------------------------------------------------------------------
message("--- Step 4: Build common human_symbol gene space ---")

human_genes <- rownames(expr_human)
dog_genes   <- rownames(expr_dog)

cs_pairs_common <- cs_broad %>%
  filter(human_symbol %in% human_genes,
         dog_symbol   %in% dog_genes)

if (nrow(cs_pairs_common) < 50) {
  stop("Too few cross-species pairs present in both matrices (", nrow(cs_pairs_common), ").")
}

## human-side matrix on shared human symbols
genes_common_h <- sort(unique(cs_pairs_common$human_symbol))
expr_human_cs  <- expr_human[genes_common_h, , drop = FALSE]

## map dog -> human_symbol space
df_dog <- as.data.frame(expr_dog) %>%
  tibble::rownames_to_column("dog_symbol")

df_dog_mapped <- cs_pairs_common %>%
  select(human_symbol, dog_symbol) %>%
  distinct() %>%
  inner_join(df_dog, by = "dog_symbol")

df_dog_h <- df_dog_mapped %>%
  group_by(human_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

expr_dog_cs <- df_dog_h %>%
  as.data.frame()
rownames(expr_dog_cs) <- expr_dog_cs$human_symbol
expr_dog_cs <- expr_dog_cs %>%
  select(-human_symbol) %>%
  as.matrix()

## final shared genes
genes_final <- intersect(rownames(expr_human_cs), rownames(expr_dog_cs))
expr_human_cs <- expr_human_cs[genes_final, , drop = FALSE]
expr_dog_cs   <- expr_dog_cs[genes_final, , drop = FALSE]

message("  Final shared genes for UMAP: ", length(genes_final))

## ---------------------------------------------------------------------
## 5. gene-wise z-scores within each cohort
## ---------------------------------------------------------------------
message("--- Step 5: Gene-wise z-scores per dataset ---")
expr_human_z <- zscore_by_gene(expr_human_cs)
expr_dog_z   <- zscore_by_gene(expr_dog_cs)

## ---------------------------------------------------------------------
## 6. combine samples and annotation
## ---------------------------------------------------------------------
message("--- Step 6: Combine samples & annotation ---")

samples_h <- colnames(expr_human_z)
samples_d <- colnames(expr_dog_z)

annot_h <- meta_human %>%
  filter(sample_id %in% samples_h) %>%
  mutate(dataset = "GSE56315", species = "Human") %>%
  select(sample_id, dataset, species, condition)

annot_d <- meta_dog %>%
  filter(sample_id %in% samples_d) %>%
  mutate(dataset = "GSE30881", species = "Dog") %>%
  select(sample_id, dataset, species, condition)

annot_all <- bind_rows(annot_h, annot_d)

expr_all_z <- cbind(expr_human_z, expr_dog_z)
expr_all_z <- expr_all_z[, annot_all$sample_id, drop = FALSE]

annot_all <- annot_all %>%
  mutate(
    species_condition = case_when(
      species == "Human" & condition == "Tumor"  ~ "Human DLBCL",
      species == "Human" & condition == "Normal" ~ "Human normal LN",
      species == "Dog"   & condition == "Tumor"  ~ "Canine DLBCL",
      species == "Dog"   & condition == "Normal" ~ "Canine normal LN",
      species == "Dog"   & condition == "Unknown" ~ "Canine DLBCL (unlabelled)",
      TRUE ~ NA_character_
    ),
    species_condition = factor(
      species_condition,
      levels = c(
        "Human normal LN",
        "Human DLBCL",
        "Canine normal LN",
        "Canine DLBCL",
        "Canine DLBCL (unlabelled)"
      )
    )
  )

message("Species-condition counts in plotted data:")
print(table(annot_all$species_condition, useNA = "ifany"))

## ---------------------------------------------------------------------
## 7. PCA -> UMAP
## ---------------------------------------------------------------------
message("--- Step 7: PCA then UMAP ---")
set.seed(12345)

pca <- prcomp(t(expr_all_z), center = TRUE, scale. = FALSE)
n_pcs <- min(20, ncol(pca$x))
pc_scores <- pca$x[, seq_len(n_pcs), drop = FALSE]

umap_res <- uwot::umap(
  pc_scores,
  n_neighbors  = 20,
  min_dist     = 0.3,
  metric       = "euclidean",
  n_components = 2,
  verbose      = TRUE
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sample_id <- annot_all$sample_id

plot_df <- umap_df %>%
  left_join(annot_all, by = "sample_id") %>%
  filter(!is.na(species_condition))

coords_out <- file.path(aux_dir, "UMAP_crossSpecies_BROAD_PCAinput_coordinates.tsv")
readr::write_tsv(plot_df, coords_out)
message("  Saved UMAP coordinates to:\n  ", coords_out)

## ---------------------------------------------------------------------
## 8. plot
## ---------------------------------------------------------------------
message("--- Step 8: Plot UMAP ---")

col_vals <- c(
  "Human normal LN"           = "#1f78b4",
  "Human DLBCL"               = "#a6cee3",
  "Canine normal LN"          = "#33a02c",
  "Canine DLBCL"              = "#b2df8a",
  "Canine DLBCL (unlabelled)" = "#999999"
)

p <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
  geom_point(
    aes(color = species_condition, shape = dataset),
    size = 2.9, alpha = 0.95
  ) +
  scale_color_manual(
    name = NULL,
    values = col_vals,
    drop = FALSE,
    na.translate = FALSE
  ) +
  scale_shape_manual(
    name   = "Dataset",
    values = c(GSE30881 = 17, GSE56315 = 16)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 16, size = 4, alpha = 1)
    ),
    shape = guide_legend(
      override.aes = list(color = "black", size = 4, alpha = 1)
    )
  ) +
  labs(
    title = "Cross-species DLBCL landscape in BROAD module space",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.box       = "vertical",
    plot.title       = element_text(face = "bold", hjust = 0,
                                    margin = margin(b = 6)),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )

out_png <- file.path(fig_dir, "Fig1A_crossSpecies_UMAP_BROAD_PCAinput.png")
out_pdf <- file.path(fig_dir, "Fig1A_crossSpecies_UMAP_BROAD_PCAinput.pdf")

ggsave(out_png, p, width = 7.5, height = 6.2, dpi = 400)
ggsave(out_pdf, p, width = 7.5, height = 6.2)

message("Saved UMAP figure to:")
message("  ", out_png)
message("  ", out_pdf)
message("=== Phase 3E UMAP (PCA->UMAP) completed successfully. ===")
