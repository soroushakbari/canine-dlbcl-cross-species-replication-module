## 41_phase10B_GSE171272_EVmiRNA_DE.R
## Differential expression for GSE171272 exosomal miRNA (DLBCL vs Control)

message("=== Phase 10B: GSE171272 EV-miRNA DE (DLBCL vs Control) ===")

required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "limma")

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
  library(limma)
})

## paths
project_root <- normalizePath(file.path(".."), winslash = "/", mustWork = TRUE)
message("Project root: ", project_root)

proc_dir <- file.path(project_root, "data", "processed")
meta_dir <- file.path(project_root, "metadata")
res_dir  <- file.path(project_root, "results", "tables", "miRNA")

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(proc_dir, "GSE171272_EVmiRNA_expr_log2.tsv")
tmpl_path <- file.path(meta_dir, "GSE171272_sample_group_template.csv")

if (!file.exists(expr_path)) {
  stop("Expression matrix not found at:\n  ", expr_path)
}
if (!file.exists(tmpl_path)) {
  stop("Group template not found at:\n  ", tmpl_path,
       "\nRun Phase 10A first and fill the 'group' column.")
}

## ------------------------------------------------------------------
## 1. Load expr + template
## ------------------------------------------------------------------
expr_tbl <- readr::read_tsv(expr_path, show_col_types = FALSE)
meta     <- readr::read_csv(tmpl_path, show_col_types = FALSE)

if (!all(c("sample_id", "group") %in% colnames(meta))) {
  stop("Template must contain 'sample_id' and 'group' columns.")
}

meta_clean <- meta %>%
  filter(!is.na(group), group != "") %>%
  mutate(
    group = stringr::str_squish(group),
    group = factor(group)
  )

if (!all(c("DLBCL", "Control") %in% levels(meta_clean$group))) {
  stop(
    "Group codes must include exactly 'DLBCL' and 'Control'.\n",
    "Detected levels: ", paste(levels(meta_clean$group), collapse = ", ")
  )
}

sample_cols <- setdiff(colnames(expr_tbl), "miRNA_id")
common_samples <- intersect(sample_cols, meta_clean$sample_id)

if (length(common_samples) < 6) {
  stop(
    "Too few samples with non-missing group labels.\n",
    "Common samples between expression and template: ", length(common_samples)
  )
}

message("  Samples used in DE analysis: ", length(common_samples))

expr_sub <- expr_tbl %>%
  select(miRNA_id, all_of(common_samples))

## align meta to expression columns
meta_sub <- meta_clean %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

stopifnot(identical(meta_sub$sample_id, colnames(expr_sub)[-1]))

expr_mat <- as.matrix(expr_sub[ , -1, drop = FALSE])
rownames(expr_mat) <- expr_sub$miRNA_id
colnames(expr_mat) <- meta_sub$sample_id

rng <- range(expr_mat, na.rm = TRUE)
message(
  "  Expression range (used as input to limma): [",
  sprintf("%.2f", rng[1]), ", ",
  sprintf("%.2f", rng[2]), "]"
)

## ------------------------------------------------------------------
## 2. limma: DLBCL vs Control
## ------------------------------------------------------------------
group <- meta_sub$group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

message("  Design matrix columns: ", paste(colnames(design), collapse = ", "))

contrast_mat <- limma::makeContrasts(
  DLBCL_vs_Control = DLBCL - Control,
  levels = design
)

fit <- limma::lmFit(expr_mat, design)
fit2 <- limma::contrasts.fit(fit, contrast_mat)
fit2 <- limma::eBayes(fit2)

tt <- limma::topTable(
  fit2,
  coef     = "DLBCL_vs_Control",
  number   = Inf,
  sort.by  = "P"
) %>%
  rownames_to_column("miRNA_id")

message("  Total miRNAs tested: ", nrow(tt))

## add simple direction labels
tt <- tt %>%
  mutate(
    direction = case_when(
      logFC >  0 ~ "Up_in_DLBCL",
      logFC <  0 ~ "Down_in_DLBCL",
      TRUE       ~ "No_change"
    )
  )

out_all <- file.path(res_dir, "GSE171272_EVmiRNA_DE_full.tsv")
readr::write_tsv(tt, out_all)
message("  Saved full DE table to:\n  ", out_all)

sig <- tt %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

out_sig <- file.path(res_dir, "GSE171272_EVmiRNA_DE_sig_FDR0.05_logFC1.tsv")
readr::write_tsv(sig, out_sig)
message(
  "  Significant miRNAs (FDR < 0.05, |logFC| ≥ 1): ",
  nrow(sig),
  "\n  Saved to:\n  ",
  out_sig
)

message("=== Phase 10B completed. Next: map sig miRNAs onto module / PPI (Phase 10C). ===")
