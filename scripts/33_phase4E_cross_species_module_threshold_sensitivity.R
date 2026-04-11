## 33_phase4E_cross_species_module_threshold_sensitivity.R
## Sensitivity analysis of cross-species module to DE thresholds
##
## ورودی‌ها:
##  - results/tables/DE/DE_GSE56315_tumor_vs_normal.tsv
##  - results/tables/DE/DE_GSE30881_tumor_vs_normal.tsv
##  - metadata/orthologs_human_to_dog_all_DE_GSE56315.csv
##  - results/tables/signatures/cross_species_module_BROAD_table.tsv
##
## خروجی‌ها:
##  - results/tables/signatures/cross_species_threshold_sensitivity_summary.tsv
##  - results/figures/SFig_threshold_sensitivity_crossSpeciesModule.png / .pdf

message("=== Phase 4E: Cross-species module threshold sensitivity analysis ===")

## ------------------------------------------------------------
## 1. packages
## ------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "ggplot2", "tidyr", "purrr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) و دوباره اسکریپت را اجرا کن."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

## ------------------------------------------------------------
## 2. paths
## ------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

de_dir   <- file.path(project_root, "results", "tables", "DE")
meta_dir <- file.path(project_root, "metadata")
sig_dir  <- file.path(project_root, "results", "tables", "signatures")
fig_dir  <- file.path(project_root, "results", "figures")

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures dir:\n  ", fig_dir)
}

de_human_path <- file.path(de_dir,  "DE_GSE56315_tumor_vs_normal.tsv")
de_dog_path   <- file.path(de_dir,  "DE_GSE30881_tumor_vs_normal.tsv")

orth_path     <- file.path(meta_dir, "orthologs_human_to_dog_all_DE_GSE56315.csv")

broad_table_path <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")

for (p in c(de_human_path, de_dog_path, orth_path, broad_table_path)) {
  if (!file.exists(p)) {
    stop("Required input not found at:\n  ", p)
  }
}

## ------------------------------------------------------------
## 3. helper: detect column names
## ------------------------------------------------------------
detect_column <- function(df, patterns) {
  cols <- names(df)
  hit <- cols[
    vapply(
      cols,
      function(x) any(stringr::str_detect(tolower(x), patterns)),
      logical(1)
    )
  ]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

## ------------------------------------------------------------
## 4. load DE tables
## ------------------------------------------------------------
message("\n--- Step 1: Load DE tables ---")

de_h <- readr::read_tsv(de_human_path, show_col_types = FALSE)
de_d <- readr::read_tsv(de_dog_path,   show_col_types = FALSE)

## human: symbol, logFC, FDR
sym_h_col <- detect_column(de_h, c("symbol", "gene", "hgnc"))
log_h_col <- detect_column(de_h, "logfc")
fdr_h_col <- detect_column(de_h, c("adj", "fdr", "padj"))

if (any(is.na(c(sym_h_col, log_h_col, fdr_h_col)))) {
  stop("Could not detect symbol/logFC/FDR columns in DE_GSE56315 table.")
}

## dog: symbol, logFC, FDR
sym_d_col <- detect_column(de_d, c("symbol", "gene"))
log_d_col <- detect_column(de_d, "logfc")
fdr_d_col <- detect_column(de_d, c("adj", "fdr", "padj"))

if (any(is.na(c(sym_d_col, log_d_col, fdr_d_col)))) {
  stop("Could not detect symbol/logFC/FDR columns in DE_GSE30881 table.")
}

message("  Human DE columns: symbol = ", sym_h_col,
        " ; logFC = ", log_h_col,
        " ; FDR = ", fdr_h_col)
message("  Dog   DE columns: symbol = ", sym_d_col,
        " ; logFC = ", log_d_col,
        " ; FDR = ", fdr_d_col)

de_h2 <- de_h %>%
  dplyr::transmute(
    human_symbol = .data[[sym_h_col]] %>% stringr::str_squish(),
    logFC_h      = as.numeric(.data[[log_h_col]]),
    FDR_h        = as.numeric(.data[[fdr_h_col]])
  ) %>%
  dplyr::filter(!is.na(human_symbol), human_symbol != "")

de_d2 <- de_d %>%
  dplyr::transmute(
    dog_symbol = .data[[sym_d_col]] %>% stringr::str_squish(),
    logFC_d    = as.numeric(.data[[log_d_col]]),
    FDR_d      = as.numeric(.data[[fdr_d_col]])
  ) %>%
  dplyr::filter(!is.na(dog_symbol), dog_symbol != "")

## ------------------------------------------------------------
## 5. load ortholog mapping + BROAD ref
## ------------------------------------------------------------
message("\n--- Step 2: Load ortholog mapping and BROAD reference module ---")

orth <- readr::read_csv(orth_path, show_col_types = FALSE)

## سعی می‌کنیم نام ستون‌ها را حدس بزنیم
h_orth_col <- detect_column(orth, c("human", "hgnc", "symbol"))
d_orth_col <- detect_column(orth, c("dog", "canis", "symbol"))

if (any(is.na(c(h_orth_col, d_orth_col)))) {
  stop("Could not detect human/dog symbol columns in ortholog table.")
}

orth2 <- orth %>%
  dplyr::transmute(
    human_symbol = .data[[h_orth_col]] %>% stringr::str_squish(),
    dog_symbol   = .data[[d_orth_col]] %>% stringr::str_squish()
  ) %>%
  dplyr::filter(
    !is.na(human_symbol), human_symbol != "",
    !is.na(dog_symbol),   dog_symbol   != ""
  ) %>%
  dplyr::distinct()

message("  Ortholog mapping rows (unique pairs): ", nrow(orth2))

broad_ref <- readr::read_tsv(broad_table_path, show_col_types = FALSE)

if (!"human_symbol" %in% names(broad_ref)) {
  stop("cross_species_module_BROAD_table.tsv must contain 'human_symbol' column.")
}

ref_human_set <- unique(broad_ref$human_symbol)

message("  Reference BROAD module (human genes): ", length(ref_human_set))

## ------------------------------------------------------------
## 6. function: build cross-species module under thresholds
## ------------------------------------------------------------
build_module <- function(lfc_h_cut, lfc_d_cut,
                         fdr_h_cut = 0.05,
                         fdr_d_cut = 0.10) {
  
  h_sig <- de_h2 %>%
    dplyr::filter(
      !is.na(logFC_h), !is.na(FDR_h),
      abs(logFC_h) >= lfc_h_cut,
      FDR_h <= fdr_h_cut
    )
  
  d_sig <- de_d2 %>%
    dplyr::filter(
      !is.na(logFC_d), !is.na(FDR_d),
      abs(logFC_d) >= lfc_d_cut,
      FDR_d <= fdr_d_cut
    )
  
  joined <- orth2 %>%
    dplyr::inner_join(h_sig, by = "human_symbol") %>%
    dplyr::inner_join(d_sig, by = "dog_symbol") %>%
    dplyr::mutate(
      dir_h = dplyr::if_else(logFC_h >= 0, "up", "down"),
      dir_d = dplyr::if_else(logFC_d >= 0, "up", "down")
    ) %>%
    dplyr::filter(dir_h == dir_d)
  
  n_pairs <- nrow(joined)
  n_up    <- sum(joined$dir_h == "up")
  n_down  <- sum(joined$dir_h == "down")
  
  mod_human <- unique(joined$human_symbol)
  
  jaccard <- if (length(mod_human) == 0) {
    NA_real_
  } else {
    intersect_n <- length(intersect(mod_human, ref_human_set))
    union_n     <- length(union(mod_human, ref_human_set))
    if (union_n == 0) NA_real_ else intersect_n / union_n
  }
  
  tibble::tibble(
    lfc_h_cut = lfc_h_cut,
    lfc_d_cut = lfc_d_cut,
    fdr_h_cut = fdr_h_cut,
    fdr_d_cut = fdr_d_cut,
    n_pairs   = n_pairs,
    n_up      = n_up,
    n_down    = n_down,
    jaccard_with_BROAD = jaccard
  )
}

## ------------------------------------------------------------
## 7. grid of thresholds + run
## ------------------------------------------------------------
message("\n--- Step 3: Run threshold grid ---")

lfc_h_grid <- c(1.0, 1.25, 1.5)
lfc_d_grid <- c(0.5, 0.75, 1.0)

grid <- tidyr::expand_grid(
  lfc_h_cut = lfc_h_grid,
  lfc_d_cut = lfc_d_grid
)

sens_tbl <- purrr::map2_dfr(
  grid$lfc_h_cut,
  grid$lfc_d_cut,
  ~ build_module(lfc_h_cut = .x, lfc_d_cut = .y)
)

sens_out_path <- file.path(
  sig_dir,
  "cross_species_threshold_sensitivity_summary.tsv"
)
readr::write_tsv(sens_tbl, sens_out_path)

message("Saved threshold sensitivity summary to:")
message("  ", sens_out_path)

print(sens_tbl)

## ------------------------------------------------------------
## 8. Heatmap figure (Jaccard vs thresholds)
## ------------------------------------------------------------
message("\n--- Step 4: Build sensitivity heatmap figure ---")

plot_df <- sens_tbl %>%
  dplyr::mutate(
    lfc_h_lab = paste0("|logFC| ≥ ", lfc_h_cut),
    lfc_d_lab = paste0("|logFC| ≥ ", lfc_d_cut)
  )

p <- ggplot(plot_df,
            aes(x = lfc_d_lab, y = lfc_h_lab, fill = jaccard_with_BROAD)) +
  geom_tile(color = "grey80") +
  geom_text(
    aes(label = ifelse(is.na(n_pairs), "0", as.character(n_pairs))),
    size = 3
  ) +
  scale_fill_gradient(
    name = "Jaccard\nvs BROAD",
    low  = "#eff3ff",
    high = "#084594",
    na.value = "white",
    limits = c(0, 1)
  ) +
  labs(
    x = "Dog DE threshold",
    y = "Human DE threshold",
    title = "Robustness of the cross-species DLBCL module to DE cut-offs",
    subtitle = "Tiles coloured by Jaccard overlap with canonical BROAD module; numbers show # ortholog pairs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10)
  )

out_png <- file.path(
  fig_dir,
  "SFig_threshold_sensitivity_crossSpeciesModule.png"
)
out_pdf <- file.path(
  fig_dir,
  "SFig_threshold_sensitivity_crossSpeciesModule.pdf"
)

ggsave(out_png, p, width = 7.0, height = 4.8, dpi = 400)
ggsave(out_pdf, p, width = 7.0, height = 4.8)

message("Saved sensitivity figure to:")
message("  ", out_png)
message("  ", out_pdf)

message("=== Phase 4E threshold sensitivity completed successfully. ===")
