## 45_phase9_Table1_cohort_overview.R
## Build final Table 1: overview of human/dog cohorts used in the study

message("=== Phase 9: Build Table 1 – Cohort overview ===")

required_pkgs <- c("tibble", "dplyr", "readr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('"%s"', required_pkgs), collapse = ", "),
      ")) و دوباره اسکریپت را اجرا کن."
    )
  }
}

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(readr)
})

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)

message("Project root: ", project_root)

out_dir <- file.path(project_root, "results", "tables", "manuscript")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created manuscript tables directory:\n  ", out_dir)
}

out_tsv <- file.path(out_dir, "Table1_cohorts_overview.tsv")

## --------------------------------------------------------------------
## 2. Hard-coded cohort metadata (based on project logs)
##    - اینجا از عددهایی استفاده می‌کنیم که در آنالیز واقعاً به‌دست آوردیم
## --------------------------------------------------------------------

# Discovery human DLBCL (tumour vs normal)
row_GSE56315 <- tibble(
  Cohort      = "GSE56315",
  Species     = "Human",
  Data_type   = "Bulk tumour biopsy gene expression (microarray)",
  N_samples   = 88L,
  Design      = "55 DLBCL vs 33 non-malignant lymphoid controls",
  Outcome     = "No survival; discovery cohort for disease module",
  Primary_role = "Define human DLBCL differential expression and seed the cross-species module"
)

# Canine tumour vs normal lymph node
row_GSE30881 <- tibble(
  Cohort      = "GSE30881",
  Species     = "Dog",
  Data_type   = "Lymph node gene expression (microarray)",
  N_samples   = 33L,
  Design      = "23 canine DLBCL vs 10 histologically normal lymph nodes",
  Outcome     = "No survival; case–control design",
  Primary_role = "Define canine DLBCL differential expression and build human–dog cross-species module"
)

# Canine CHOP-treated RNA-seq cohort (PFS + hubs)
row_GSE130874 <- tibble(
  Cohort      = "GSE130874",
  Species     = "Dog",
  Data_type   = "Bulk RNA-seq (pretreatment DLBCL under CHOP)",
  N_samples   = 23L,
  Design      = "23 dogs with newly diagnosed DLBCL, sampled before CHOP-based chemotherapy",
  Outcome     = "Progression-free survival under CHOP (21 events; see text)",
  Primary_role = "Validate cross-species module activity, link module scores to PARP1/TOP2A expression and PFS"
)

# Additional canine B-cell lymphoma expression cohort
row_GSE39365 <- tibble(
  Cohort      = "GSE39365",
  Species     = "Dog",
  Data_type   = "B-cell lymphoma gene expression (microarray)",
  N_samples   = 36L,
  Design      = "36 canine B-cell lymphoma biopsies (no matched normal lymph nodes)",
  Outcome     = "No survival; tumour-only",
  Primary_role = "Support canine-side module projection and PPI context for cross-species network"
)

# Human EV-miRNA cohort
row_GSE171272 <- tibble(
  Cohort      = "GSE171272",
  Species     = "Human",
  Data_type   = "Plasma extracellular-vesicle microRNA profiling",
  N_samples   = 15L,
  Design      = "DLBCL patients vs non-malignant controls (EV-miRNA)",
  Outcome     = "No time-to-event; case–control EV-miRNA design",
  Primary_role = "Identify differentially expressed EV-miRNAs targeting the cross-species module and PPI hubs"
)

# Human TCGA-DLBC RNA-seq cohort
row_TCGA_DLBC <- tibble(
  Cohort      = "TCGA-DLBC",
  Species     = "Human",
  Data_type   = "Bulk RNA-seq (TCGA diffuse large B-cell lymphoma project)",
  N_samples   = 47L,  # تعداد کیس‌های با OS قابل‌استفاده در آنالیز ما
  Design      = "De novo DLBCL tumour biopsies; no matched normal tissue",
  Outcome     = "Overall survival (limited events; underpowered)",
  Primary_role = "Evaluate cross-species module scores as a prognostic signal in an independent RNA-seq cohort"
)

# Large human OS validation cohort
row_GSE31312 <- tibble(
  Cohort      = "GSE31312",
  Species     = "Human",
  Data_type   = "Bulk tumour biopsy gene expression (microarray)",
  N_samples   = 498L,
  Design      = "498 de novo DLBCL patients treated with immunochemotherapy",
  Outcome     = "Overall survival available for 470 patients (170 deaths) from companion clinical dataset",
  Primary_role = "Large-scale validation of cross-species module scores as a prognostic marker"
)

## اگر نخواستی GSE39365 را در متن اصلی بیاوری، می‌تونی بعداً این ردیف را حذف کنی؛
## اما اینجا می‌آوریم تا Table 1 واقعاً تصویری کامل از کوهورت‌ها بدهد.

table1 <- bind_rows(
  row_GSE56315,
  row_GSE30881,
  row_GSE130874,
  row_GSE39365,
  row_TCGA_DLBC,
  row_GSE31312,
  row_GSE171272
) %>%
  # ترتیب منطقی: اول کشف انسان/سگ، بعد کوهورت‌های survival، بعد EV-miRNA
  dplyr::mutate(
    Cohort = factor(
      Cohort,
      levels = c("GSE56315", "GSE30881", "GSE130874",
                 "GSE39365", "TCGA-DLBC", "GSE31312", "GSE171272")
    )
  ) %>%
  dplyr::arrange(Cohort) %>%
  dplyr::mutate(
    Cohort = as.character(Cohort)
  )

readr::write_tsv(table1, out_tsv)

message("Saved Table 1 to:")
message("  ", out_tsv)

message("=== Table 1 build completed successfully. ===")
