## ==========================================
## 04_phase2C_rnaseq_canine_GSE130874.R
## Phase 2C - Download & process canine RNA-seq
## Dataset: GSE130874 (Mee 2022, 25 canine lymphoma pre-CHOP, PFS)
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

# ---------- 1) پکیج‌ها ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkg_bioc <- c("GEOquery", "Biobase", "SummarizedExperiment", "S4Vectors", "DESeq2")
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
  library(DESeq2)
})


# اینونتوری فعلی
dataset_inventory <- readr::read_csv(inventory_path, show_col_types = FALSE)

# ---------- 2) تابع QC برای RNA-seq (مثل TCGA) ----------

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

# ---------- 3) دانلود فایل‌های GSE130874 (supplementary) ----------

message("\n==============================")
message("Phase 2C: Canine RNA-seq - GSE130874")
message("==============================")

gse_id <- "GSE130874"
gse_dir <- file.path(raw_geo_dir, gse_id)
dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)

message("Downloading supplementary files for ", gse_id, " ...")

# آرگومان اول positionally = شناسه GEO؛ لازم نیست اسم‌گذاری کنی
GEOquery::getGEOSuppFiles(
  gse_id,
  baseDir = raw_geo_dir
)

# مسیر فایل‌های count و FPKM (طبق GEO)
count_file_gz <- file.path(gse_dir, "GSE130874_read_count_matrix.tsv.gz")
fpkm_file_gz  <- file.path(gse_dir, "GSE130874_FPKM_matrix.tsv.gz")

if (!file.exists(count_file_gz)) {
  stop("Count matrix file not found: ", count_file_gz,
       "\nCheck downloaded files under: ", gse_dir)
}

# ---------- 4) خواندن ماتریس counts (به صورت general) ----------

message("Reading count matrix from:\n  ", count_file_gz)

counts_tbl <- readr::read_tsv(
  file = count_file_gz,
  progress = TRUE,
  show_col_types = FALSE
)

if (ncol(counts_tbl) < 2) {
  stop("Count matrix has <2 columns; check file structure.")
}

# تشخیص خودکار ستون‌های عددی (نمونه‌ها) و غیرعدد (IDها/annotation)
is_num <- sapply(counts_tbl, is.numeric)

if (!any(is_num)) {
  stop("No numeric columns found in count matrix. Check file content.")
}

id_cols <- which(!is_num)

if (length(id_cols) == 0) {
  stop("No non-numeric column to use as gene_id; inspect the file structure.")
}

# اولین ستون غیرعدد را به عنوان gene_id در نظر می‌گیریم
gene_id_col_index <- id_cols[1]
gene_id_col_name  <- colnames(counts_tbl)[gene_id_col_index]

message("Using column as gene_id: ", gene_id_col_name)

gene_ids <- counts_tbl[[gene_id_col_index]]

# فقط ستون‌های عددی را به عنوان counts نگه می‌داریم
counts_numeric_df <- counts_tbl[, is_num, drop = FALSE]

# تبدیل مطمئن به ماتریس عددی
counts_mat <- counts_numeric_df %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(counts_mat) <- gene_ids

message("Counts matrix: ",
        nrow(counts_mat), " genes x ",
        ncol(counts_mat), " samples.")

# ---------- 5) فیلتر ژن‌های خیلی low-count و VST ----------

keep_genes <- rowSums(counts_mat) >= 10
counts_filt <- counts_mat[keep_genes, , drop = FALSE]

message("After filtering: ",
        nrow(counts_filt), " genes with sum(counts) >= 10.")

# ---------- 6) متادیتای نمونه‌ها از خود GEO + مپ هوشمند ----------

message("Fetching sample metadata via getGEO(", gse_id, ") ...")

gse_list <- GEOquery::getGEO(
  gse_id,
  GSEMatrix = TRUE
)

if (length(gse_list) > 1) {
  ns <- sapply(gse_list, function(es) ncol(Biobase::exprs(es)))
  idx <- which.max(ns)
  gset <- gse_list[[idx]]
} else {
  gset <- gse_list[[1]]
}

meta_raw <- Biobase::pData(gset) %>%
  as.data.frame()

# GSM IDs در rownames
meta_raw$gsm_id <- rownames(meta_raw)

# نمونه‌هایی که در ماتریس counts داریم
sample_ids <- colnames(counts_filt)

# سعی می‌کنیم بهترین ستون متادیتا را برای مَچ پیدا کنیم
all_cols <- unique(c("gsm_id", colnames(meta_raw)))

best_col <- NULL
best_match <- 0

for (col in all_cols) {
  vals <- as.character(meta_raw[[col]])
  m <- sum(vals %in% sample_ids)
  if (!is.na(m) && m > best_match) {
    best_match <- m
    best_col <- col
  }
}

if (!is.null(best_col) && best_match > 0) {
  message("Best metadata mapping column: ", best_col,
          " (", best_match, " / ", length(sample_ids), " matches).")
  
  vals <- as.character(meta_raw[[best_col]])
  
  # فقط ردیف‌هایی که واقعاً مچ می‌شوند را نگه می‌داریم
  keep <- vals %in% sample_ids
  meta_df <- meta_raw[keep, , drop = FALSE]
  
  # ترتیب را براساس sample_ids تنظیم می‌کنیم
  meta_df <- meta_df[match(sample_ids, vals[keep]), , drop = FALSE]
  
  meta_df$sample_id <- sample_ids
} else {
  warning("Could not match sample IDs to any metadata column. Using minimal colData with only sample_id.")
  meta_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
}

# چک نهایی: باید تعداد ردیف‌های متادیتا = تعداد نمونه‌ها باشد
if (nrow(meta_df) != length(sample_ids)) {
  stop("After metadata alignment, nrow(meta_df) = ", nrow(meta_df),
       " but number of samples = ", length(sample_ids),
       ". Fix mapping logic for GSE130874.")
}

# ---------- 7) DESeq2 VST ----------

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = as.data.frame(meta_df),
  design    = ~ 1
)

dds <- DESeq2::estimateSizeFactors(dds)
vst_obj <- DESeq2::vst(dds, blind = TRUE)
vst_mat <- SummarizedExperiment::assay(vst_obj)

message("VST matrix: ",
        nrow(vst_mat), " genes x ",
        ncol(vst_mat), " samples.")

# ---------- 8) QC metrics ----------

qc <- compute_qc_metrics_rnaseq(vst_mat)

meta_df <- meta_df %>%
  dplyr::left_join(qc, by = c("sample_id" = "sample_id"))

meta_path <- file.path(metadata_dir, "GSE130874_sample_metadata.csv")
readr::write_csv(meta_df, meta_path)

# ---------- 9) Feature annotation ----------

# در این مرحله فقط gene_id داریم؛ gene_symbol را بعداً با biomaRt (canine Ensembl) پر می‌کنیم
feature_df <- tibble::tibble(
  gene_id     = rownames(counts_filt),
  gene_symbol = NA_character_
)

feature_path <- file.path(metadata_dir, "GSE130874_feature_annotation.csv")
readr::write_csv(feature_df, feature_path)

# ---------- 10) ذخیره counts + VST ----------

counts_tbl_out <- counts_filt %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

counts_path <- file.path(processed_dir, "GSE130874_counts.tsv")
readr::write_tsv(counts_tbl_out, counts_path)

vst_tbl_out <- vst_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

vst_path <- file.path(processed_dir, "GSE130874_vst.tsv")
readr::write_tsv(vst_tbl_out, vst_path)

# ---------- 11) SummarizedExperiment نهایی ----------

se_final <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    counts = counts_filt,
    vst    = vst_mat
  ),
  rowData = S4Vectors::DataFrame(feature_df),
  colData = S4Vectors::DataFrame(meta_df)
)

se_final_path <- file.path(processed_dir, "GSE130874_SummarizedExperiment_processed.rds")
saveRDS(se_final, se_final_path)

message("Saved processed SE to:\n  ", se_final_path)

# ---------- 12) آپدیت اینونتوری ----------

if ("dataset_id" %in% colnames(dataset_inventory)) {
  # اگر قبلاً ردیفی با Mee_2022_RNAseq_25cases یا GSE130874 داشتیم، همان را به‌روز می‌کنیم
  idx <- which(dataset_inventory$dataset_id %in% c("Mee_2022_RNAseq_25cases", "GSE130874"))
  
  if (length(idx) == 0) {
    # اضافه کردن ردیف جدید
    new_row <- tibble::tibble(
      dataset_id         = "GSE130874",
      species            = "dog",
      role               = "discovery_outcome",
      source_type        = "GEO",
      description        = "Canine lymphoma RNA-seq (25 cases pre-CHOP) with PFS (Mee 2022, BMC Res Notes)",
      has_tumor_normal   = FALSE,
      has_outcome        = TRUE,
      outcome_type       = "PFS",
      is_cell_line       = FALSE,
      platform           = "RNA-seq (GSE130874 read counts; GPL22370)",
      n_samples          = ncol(counts_filt),
      include_discovery  = TRUE,
      include_validation = TRUE,
      notes              = "Primary canine outcome cohort (PFS) for module-outcome analysis."
    )
    dataset_inventory <- dplyr::bind_rows(dataset_inventory, new_row)
    message("Added new inventory row for GSE130874.")
  } else {
    # به‌روزرسانی ردیف موجود (اولی را هدف قرار می‌دهیم)
    i <- idx[1]
    dataset_inventory$dataset_id[i]   <- "GSE130874"
    dataset_inventory$species[i]      <- "dog"
    dataset_inventory$role[i]         <- "discovery_outcome"
    dataset_inventory$source_type[i]  <- "GEO"
    dataset_inventory$description[i]  <- "Canine lymphoma RNA-seq (25 cases pre-CHOP) with PFS (Mee 2022, BMC Res Notes)"
    dataset_inventory$has_tumor_normal[i] <- FALSE
    dataset_inventory$has_outcome[i]      <- TRUE
    dataset_inventory$outcome_type[i]     <- "PFS"
    dataset_inventory$is_cell_line[i]     <- FALSE
    dataset_inventory$platform[i]         <- "RNA-seq (GSE130874 read counts; GPL22370)"
    dataset_inventory$n_samples[i]        <- ncol(counts_filt)
    dataset_inventory$include_discovery[i]  <- TRUE
    dataset_inventory$include_validation[i] <- TRUE
    dataset_inventory$notes[i]           <- "Primary canine outcome cohort (PFS) for module-outcome analysis."
    message("Updated existing inventory row for GSE130874.")
  }
  
  readr::write_csv(dataset_inventory, inventory_path)
} else {
  warning("Inventory has no 'dataset_id' column. Skipping inventory update.")
}

cat("\n[Phase 2C] GSE130874 canine RNA-seq processing completed.\n")
cat("Counts: ", counts_path, "\n", sep = "")
cat("VST: ", vst_path, "\n", sep = "")
cat("Feature annotation: ", feature_path, "\n", sep = "")
cat("Sample metadata: ", meta_path, "\n", sep = "")
cat("Processed SE: ", se_final_path, "\n\n", sep = "")
