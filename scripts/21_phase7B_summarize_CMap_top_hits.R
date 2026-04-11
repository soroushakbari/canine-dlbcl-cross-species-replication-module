## =====================================================================
## Phase 7B: Summarize and categorize top CMap hits (cross-species BROAD)
##   - Input:
##       results/tables/Drug/CMap_queryl1k_cross_species_BROAD_drug_top50_negative.tsv
##   - Output:
##       results/tables/Drug/CMap_queryl1k_cross_species_BROAD_drug_top50_negative_annotated.tsv
##       results/tables/Drug/CMap_queryl1k_cross_species_BROAD_drug_top50_negative_for_manuscript.tsv
## =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) file.path(project_root, ...)

drug_dir <- path_proj("results", "tables", "Drug")
dir.create(drug_dir, recursive = TRUE, showWarnings = FALSE)

input_top50 <- file.path(
  drug_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_top50_negative.tsv"
)

if (!file.exists(input_top50)) {
  stop("Top50 file not found: ", input_top50,
       "\nRun 20_phase7_CLUE_query_result_to_drug_rank.R first.")
}

top50 <- readr::read_tsv(input_top50, show_col_types = FALSE)

message("[CMap-summary] Loaded top50: ", nrow(top50), " rows.")

## ---------------------------------------------------------------------
## 1. یک ستون ساده‌شده برای MOA/targets بسازیم
## ---------------------------------------------------------------------

## helper: رشته‌ها رو کمی تمیز کنیم
clean_text <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

top50 <- top50 %>%
  mutate(
    moa      = if (!"moa"     %in% colnames(.)) NA_character_ else clean_text(moa),
    targets  = if (!"targets" %in% colnames(.)) NA_character_ else clean_text(targets)
  )

## ---------------------------------------------------------------------
## 2. دسته‌بندی مکانیکی (MOA_category) با pattern matching ساده
## ---------------------------------------------------------------------

## همه حروف رو lowercase می‌کنیم برای جستجوی راحت‌تر
top50 <- top50 %>%
  mutate(
    moa_lower     = tolower(ifelse(is.na(moa), "", moa)),
    targets_lower = tolower(ifelse(is.na(targets), "", targets))
  )

## تابع برای تشخیص دسته‌ی MOA
assign_moa_category <- function(moa_low, target_low) {
  txt <- paste(moa_low, target_low, sep = " ")
  
  ## BCR / BTK / SYK / PI3K و غیره
  if (str_detect(txt, "btk") | str_detect(txt, "bruton") |
      str_detect(txt, "b-cell receptor") | str_detect(txt, "bcr")) {
    return("BCR/BTK-pathway inhibitor")
  }
  
  if (str_detect(txt, "pi3k") | str_detect(txt, "akt") |
      str_detect(txt, "mtor")) {
    return("PI3K/AKT/mTOR pathway inhibitor")
  }
  
  if (str_detect(txt, "bcl2") | str_detect(txt, "bcl-2") |
      str_detect(txt, "apoptosis regulator")) {
    return("BCL2/apoptosis modulator")
  }
  
  ## epigenetic
  if (str_detect(txt, "hdac") | str_detect(txt, "histone deacetylase")) {
    return("HDAC/epigenetic modulator")
  }
  
  if (str_detect(txt, "dnmt") | str_detect(txt, "dna methyltransferase")) {
    return("DNMT/epigenetic modulator")
  }
  
  ## proteasome
  if (str_detect(txt, "proteasome")) {
    return("Proteasome inhibitor")
  }
  
  ## kinase broad
  if (str_detect(txt, "kinase") | str_detect(txt, "tyrosine kinase") |
      str_detect(txt, "tk inhibitor")) {
    return("Multi-kinase inhibitor")
  }
  
  ## cell cycle / CDK
  if (str_detect(txt, "cdk") | str_detect(txt, "cyclin-dependent kinase")) {
    return("CDK/cell-cycle inhibitor")
  }
  
  ## topoisomerase / classic chemo
  if (str_detect(txt, "topoisomerase") | str_detect(txt, "topo ii") |
      str_detect(txt, "doxorubicin") | str_detect(txt, "etoposide")) {
    return("Topoisomerase/chemotherapy-like agent")
  }
  
  ## immunomodulatory / steroid-like
  if (str_detect(txt, "glucocorticoid") | str_detect(txt, "steroid") |
      str_detect(txt, "immunomodulatory")) {
    return("Immunomodulatory/steroid-related")
  }
  
  ## اگر هیچ کدام تشخیص داده نشد:
  if (nchar(txt) == 0) {
    return("Unknown/NA")
  } else {
    return("Other/unspecified")
  }
}

top50 <- top50 %>%
  rowwise() %>%
  mutate(
    MOA_category = assign_moa_category(moa_lower, targets_lower)
  ) %>%
  ungroup()

## ---------------------------------------------------------------------
## 3. مرتب‌سازی نهایی برای مقاله
## ---------------------------------------------------------------------

## مرتب‌سازی: اول بر اساس min_score (بیشترین منفی)، بعد بر اساس MOA_category
top50_for_manuscript <- top50 %>%
  arrange(min_score) %>%
  mutate(
    rank = row_number()
  ) %>%
  select(
    rank,
    pert_name,
    MOA_category,
    moa,
    targets,
    n_signatures,
    min_score,
    median_score,
    mean_score,
    best_fdr_q_nlog10,
    best_cell
  )

## ---------------------------------------------------------------------
## 4. خروجی‌ها
## ---------------------------------------------------------------------

annotated_path <- file.path(
  drug_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_top50_negative_annotated.tsv"
)

manuscript_path <- file.path(
  drug_dir,
  "CMap_queryl1k_cross_species_BROAD_drug_top50_negative_for_manuscript.tsv"
)

readr::write_tsv(top50, annotated_path)
readr::write_tsv(top50_for_manuscript, manuscript_path)

message("[CMap-summary] Annotated top50 written to: ", annotated_path)
message("[CMap-summary] Manuscript-style top50 table written to: ", manuscript_path)
message("=== [Phase 7B - Summarized & categorized top CMap hits] DONE ===")
