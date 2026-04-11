## ==========================================
## 01_phase1_dataset_inventory.R
## Phase 1 - Build dataset inventory for DLBCL project
## ==========================================

project_root <- "D:/Research/My Articles/DLBCL drug"
metadata_dir <- file.path(project_root, "metadata")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

inventory_path <- file.path(metadata_dir, "dataset_inventory_phase1.csv")

# پکیج‌ها
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
suppressPackageStartupMessages({
  library(tidyverse)
})

# جدول دیتاست‌ها
datasets_phase1 <- tibble::tibble(
  dataset_id    = c(
    "TCGA_DLBC",                  # human
    "GSE10846",
    "GSE56315",
    "GSE32018",
    "GSE30881",                   # dog
    "Mee_2022_RNAseq_25cases",
    "GSE43664",
    "GSE39365",
    "GSE224867",                  # cell lines / models
    "GSE149926"
  ),
  species      = c(
    "human",
    "human",
    "human",
    "human",
    "dog",
    "dog",
    "dog",
    "dog",
    "human",
    "dog"
  ),
  role         = c(
    "discovery_main",   # TCGA
    "discovery_main",   # GSE10846
    "supporting_tn",    # tumor vs normal
    "supporting_tn",
    "discovery_main",   # canine tumor vs normal
    "discovery_outcome",
    "supporting_other",
    "supporting_other",
    "cell_line_model",
    "cell_line_model"
  ),
  source_type  = c(
    "TCGA",
    "GEO",
    "GEO",
    "GEO",
    "GEO",
    "paper_or_GEO_TBD",
    "GEO",
    "GEO",
    "GEO",
    "GEO"
  ),
  description  = c(
    "TCGA DLBC (human DLBCL cohort with survival/outcome data)",
    "Human DLBCL (R-CHOP) with clinical outcome (e.g., OS/PFS/response)",
    "Human DLBCL vs normal tonsil/lymph node (tumor-normal contrast)",
    "Human DLBCL vs normal / B-cell controls (tumor-normal contrast)",
    "Canine DLBCL vs healthy lymph node; cross-species co-expression study",
    "RNA-seq, ~25 canine lymphoma cases pre-CHOP with PFS/response (Mee 2022)",
    "Canine B-cell lymphoma microarray cohort",
    "Canine B-cell lymphoma microarray cohort",
    "Human DLBCL cell lines ± chemotherapy (response/resistance model)",
    "CLBL-1 canine lymphoma, doxorubicin-resistant vs sensitive (SuperSeries)"
  ),
  has_tumor_normal = c(
    FALSE,  # TCGA
    FALSE,  # GSE10846
    TRUE,   # GSE56315
    TRUE,   # GSE32018
    TRUE,   # GSE30881
    FALSE,  # Mee 2022
    NA,
    NA,
    FALSE,  # cell lines
    FALSE
  ),
  has_outcome = c(
    TRUE,   # TCGA DLBC
    TRUE,   # GSE10846
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    NA,
    NA,
    FALSE,
    FALSE
  ),
  outcome_type = c(
    "OS/PFS",
    "OS/PFS/response",
    NA,
    NA,
    NA,
    "PFS/response",
    NA,
    NA,
    NA,
    NA
  ),
  is_cell_line = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE
  ),
  platform     = NA_character_,
  n_samples    = NA_real_,
  include_discovery  = c(
    TRUE,   # TCGA
    TRUE,   # GSE10846
    TRUE,   # GSE56315
    TRUE,   # GSE32018
    TRUE,   # GSE30881
    TRUE,   # Mee 2022
    FALSE,  # supporting
    FALSE,
    FALSE,  # cell line models
    FALSE
  ),
  include_validation = c(
    TRUE,   # TCGA
    TRUE,   # GSE10846
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  notes        = c(
    "Use TCGAbiolinks to pull expression + clinical data. Check sample types.",
    "Classic R-CHOP DLBCL cohort; key for survival/outcome link.",
    "Tumor-normal contrast; good for defining disease signature.",
    "Additional tumor-normal cohort; check platform and sample size.",
    "Key canine tumor-normal cohort; important for cross-species mapping.",
    "Need exact accession / GEO ID from paper; but conceptually key PFS cohort.",
    "Check original paper for exact design (tumor vs normal? outcome?).",
    "Check original paper; likely tumor-only B-cell lymphoma.",
    "DLBCL cell line response/resistance dataset; useful for checking module.",
    "SuperSeries with CLBL-1 resistant vs sensitive; good for cross-check."
  )
)

readr::write_csv(datasets_phase1, inventory_path)

cat("\n[Phase 1] Dataset inventory written to:\n")
cat(inventory_path, "\n\n")


# تابع اضافه کردن دیتاست جدید در آینده
add_dataset_to_inventory <- function(dataset_id,
                                     species,
                                     role = "supporting_other",
                                     source_type = "GEO",
                                     description = "",
                                     has_tumor_normal = NA,
                                     has_outcome = NA,
                                     outcome_type = NA,
                                     is_cell_line = FALSE,
                                     platform = NA,
                                     n_samples = NA,
                                     include_discovery = FALSE,
                                     include_validation = TRUE,
                                     notes = "") {
  
  if (!file.exists(inventory_path)) {
    stop("Inventory file does not exist yet: ", inventory_path)
  }
  
  inv <- readr::read_csv(inventory_path,
                         show_col_types = FALSE,
                         progress = FALSE)
  
  if (dataset_id %in% inv$dataset_id) {
    warning("Dataset ID already in inventory: ", dataset_id)
  }
  
  new_row <- tibble::tibble(
    dataset_id         = dataset_id,
    species            = species,
    role               = role,
    source_type        = source_type,
    description        = description,
    has_tumor_normal   = has_tumor_normal,
    has_outcome        = has_outcome,
    outcome_type       = outcome_type,
    is_cell_line       = is_cell_line,
    platform           = platform,
    n_samples          = n_samples,
    include_discovery  = include_discovery,
    include_validation = include_validation,
    notes              = notes
  )
  
  inv2 <- dplyr::bind_rows(inv, new_row)
  
  readr::write_csv(inv2, inventory_path)
  cat("Added dataset:", dataset_id, "→", inventory_path, "\n")
  invisible(inv2)
}

cat("Helper function `add_dataset_to_inventory()` is defined in this script.\n\n")
