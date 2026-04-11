## ==========================================
## 03_phase2B_rnaseq_TCGA_DLBC.R
## Phase 2B - Download & process TCGA-DLBC RNA-seq
## Project: DLBCL cross-species module & drug repurposing
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"

raw_tcga_dir  <- file.path(project_root, "data", "raw", "TCGA")
processed_dir <- file.path(project_root, "data", "processed")
metadata_dir  <- file.path(project_root, "metadata")

dir.create(raw_tcga_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir,  recursive = TRUE, showWarnings = FALSE)

inventory_path <- file.path(metadata_dir, "dataset_inventory_phase1.csv")
if (!file.exists(inventory_path)) {
  stop("Inventory file not found: ", inventory_path,
       "\nRun 01_phase1_dataset_inventory.R first.")
}

# ---------- 1) پکیج‌ها ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkg_bioc <- c("TCGAbiolinks", "SummarizedExperiment", "S4Vectors", "DESeq2")
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
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(DESeq2)
})

dataset_inventory <- readr::read_csv(inventory_path, show_col_types = FALSE)

# ---------- 2) توابع کمکی ----------

# حدس ستون سمبل ژن در rowData
guess_gene_symbol_column_tcga <- function(se) {
  rd <- SummarizedExperiment::rowData(se)
  cn <- colnames(rd)
  candidates <- c("external_gene_name", "gene_name", "symbol", "hgnc_symbol")
  hit <- candidates[candidates %in% cn][1]
  if (!is.na(hit)) return(hit)
  
  hits <- grep("symbol|gene", cn, ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(hits[1])
  
  stop("Could not find a gene symbol-like column in rowData(SE).")
}

# QC ساده روی RNA-seq نرمال‌شده (VST)
compute_qc_metrics_rnaseq <- function(expr_mat) {
  if (!is.matrix(expr_mat)) {
    expr_mat <- as.matrix(expr_mat)
  }
  
  libsize <- colSums(expr_mat, na.rm = TRUE)
  
  qs <- apply(expr_mat, 2, function(x) {
    stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  })
  
  cor_mat <- stats::cor(expr_mat, method = "spearman", use = "pairwise.complete.obs")
  mean_cor <- apply(cor_mat, 2, function(x) mean(x, na.rm = TRUE))
  
  qc <- tibble::tibble(
    sample_id = colnames(expr_mat),
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

# ---------- 3) Query + download + prepare TCGA-DLBC (بدون دردسر workflow.type) ----------

message("\n==============================")
message("Phase 2B: TCGA-DLBC RNA-seq")
message("==============================")

message("Querying TCGA-DLBC (Transcriptome Profiling / Gene Expression Quantification) ...")

# این همون چیزی‌ه که الان دیدی جواب می‌ده
query_dlbc <- TCGAbiolinks::GDCquery(
  project = "TCGA-DLBC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)

# اگر هیچ نتیجه‌ای نباشه، stop می‌کنیم
res_df <- TCGAbiolinks::getResults(query_dlbc)
if (nrow(res_df) == 0) {
  stop("GDCquery returned 0 results for TCGA-DLBC (Transcriptome Profiling / Gene Expression Quantification).")
}

message("Number of files found: ", nrow(res_df))

# 3.2) دانلود
TCGAbiolinks::GDCdownload(query_dlbc, directory = raw_tcga_dir)

# 3.3) SummarizedExperiment
se_dlbc <- TCGAbiolinks::GDCprepare(
  query = query_dlbc,
  directory = raw_tcga_dir
)

# ذخیره نسخه خام
se_raw_path <- file.path(processed_dir, "TCGA_DLBC_SummarizedExperiment_raw.rds")
saveRDS(se_dlbc, se_raw_path)

message("GDCprepare done. Saved raw SE to:\n  ", se_raw_path)


# ---------- 4) استخراج counts و نرمال‌سازی ----------

# فرض: assay اصلی counts است
counts_mat <- SummarizedExperiment::assay(se_dlbc)

message("Counts matrix: ",
        nrow(counts_mat), " genes x ",
        ncol(counts_mat), " samples.")

# حذف ژن‌های با شمارش خیلی پایین (noise)
keep_genes <- rowSums(counts_mat) >= 10
counts_filt <- counts_mat[keep_genes, , drop = FALSE]

message("After filtering: ",
        nrow(counts_filt), " genes with sum(counts) >= 10.")

# ساخت DESeqDataSet برای VST
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = SummarizedExperiment::colData(se_dlbc),
  design    = ~ 1
)

dds <- DESeq2::estimateSizeFactors(dds)
vst_obj <- DESeq2::vst(dds, blind = TRUE)
vst_mat <- SummarizedExperiment::assay(vst_obj)

message("VST matrix: ",
        nrow(vst_mat), " genes x ",
        ncol(vst_mat), " samples.")

# ---------- 5) Feature annotation (gene-level) ----------

rd <- SummarizedExperiment::rowData(se_dlbc)
gene_id <- rownames(counts_filt)

# تلاش برای پیدا کردن سمبل ژن
sym_col <- guess_gene_symbol_column_tcga(se_dlbc)
gene_symbol_full <- as.character(rd[[sym_col]])
gene_symbol <- gene_symbol_full[keep_genes]

feature_df <- tibble::tibble(
  gene_id    = gene_id,
  gene_symbol = gene_symbol
)

feature_path <- file.path(metadata_dir, "TCGA_DLBC_feature_annotation.csv")
readr::write_csv(feature_df, feature_path)

# ---------- 6) Sample metadata + QC ----------

# colData را به data.frame تبدیل می‌کنیم
meta_df <- SummarizedExperiment::colData(se_dlbc) %>%
  as.data.frame()

# اگر از قبل ستونی به اسم sample_id هست، اول اسمش را عوض کن که دعوا نشود
if ("sample_id" %in% colnames(meta_df)) {
  meta_df <- meta_df %>%
    dplyr::rename(sample_id_tcga = sample_id)
}

# حالا rownames (که در واقع همون colnames expr/vst هستن) رو به ستون sample_id می‌ریزیم
meta_df <- meta_df %>%
  tibble::rownames_to_column("sample_id")

# QC روی ماتریس VST
qc <- compute_qc_metrics_rnaseq(vst_mat)

meta_df <- meta_df %>%
  dplyr::left_join(qc, by = c("sample_id" = "sample_id"))

meta_path <- file.path(metadata_dir, "TCGA_DLBC_sample_metadata.csv")
readr::write_csv(meta_df, meta_path)


# ---------- 7) ذخیره ماتریس‌ها (counts + VST) ----------

# counts (فقط ژن‌های فیلترشده)
counts_tbl <- counts_filt %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::left_join(feature_df, by = "gene_id") %>%
  dplyr::relocate(gene_id, gene_symbol)

counts_path <- file.path(processed_dir, "TCGA_DLBC_counts.tsv")
readr::write_tsv(counts_tbl, counts_path)

# VST
vst_tbl <- vst_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::left_join(feature_df, by = "gene_id") %>%
  dplyr::relocate(gene_id, gene_symbol)

vst_path <- file.path(processed_dir, "TCGA_DLBC_vst.tsv")
readr::write_tsv(vst_tbl, vst_path)

# ---------- 8) ساخت SummarizedExperiment نهایی برای downstream ----------

se_final <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    counts = counts_filt,
    vst    = vst_mat
  ),
  rowData = S4Vectors::DataFrame(feature_df),
  colData = S4Vectors::DataFrame(meta_df)
)

se_final_path <- file.path(processed_dir, "TCGA_DLBC_SummarizedExperiment_processed.rds")
saveRDS(se_final, se_final_path)

message("Saved processed SE to:\n  ", se_final_path)

# ---------- 9) آپدیت اینونتوری ----------

if ("dataset_id" %in% colnames(dataset_inventory)) {
  idx <- which(dataset_inventory$dataset_id == "TCGA_DLBC")
  if (length(idx) == 1) {
    dataset_inventory$platform[idx]  <- "RNA-seq (TCGA HTSeq-Counts)"
    dataset_inventory$n_samples[idx] <- ncol(counts_mat)
    readr::write_csv(dataset_inventory, inventory_path)
    message("Updated inventory for TCGA_DLBC (platform, n_samples).")
  } else {
    warning("Could not uniquely match TCGA_DLBC in inventory.")
  }
} else {
  warning("Inventory has no 'dataset_id' column. Skipping inventory update.")
}

cat("\n[Phase 2B] TCGA-DLBC RNA-seq processing completed.\n")
cat("Counts: ", counts_path, "\n", sep = "")
cat("VST: ", vst_path, "\n", sep = "")
cat("Feature annotation: ", feature_path, "\n", sep = "")
cat("Sample metadata: ", meta_path, "\n", sep = "")
cat("Processed SE: ", se_final_path, "\n\n", sep = "")
