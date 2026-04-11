## 21_phase7C_backfill_CMap_targets_from_CLUE_metadata.R
## هدف:
##  - گشتن در CLUE_raw برای پیدا کردن فایل متادیتایی که هم pert_iname دارد هم target
##  - چسباندن این targets به جدول annotated CMap top50
##  - آپدیت کردن ستون 'targets' در:
##      results/tables/Drug/CMap_queryl1k_cross_species_BROAD_drug_top50_negative_annotated.tsv

message("=== Phase 7C: Backfill CMap targets from CLUE metadata ===")

required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "purrr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(\"", pkg, "\") and then re-run this script."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(purrr)
})

## -------------------------------------------------------------------------
## 1. مسیرها
## -------------------------------------------------------------------------

project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)

message("Project root: ", project_root)

## این همون جاییه که زیپ CLUE رو قبلاً unzip کردی
clue_root <- file.path(
  project_root, "results", "tables", "Drug", "CLUE_raw"
)

if (!dir.exists(clue_root)) {
  stop(
    "CLUE_raw directory not found at:\n", clue_root, "\n\n",
    "Make sure you've unzipped the CLUE download under results/tables/Drug/CLUE_raw/"
  )
}

path_drugs_annot <- file.path(
  project_root, "results", "tables", "Drug",
  "CMap_queryl1k_cross_species_BROAD_drug_top50_negative_annotated.tsv"
)

if (!file.exists(path_drugs_annot)) {
  stop(
    "Annotated CMap top50 file not found at:\n", path_drugs_annot, "\n\n",
    "Expected: CMap_queryl1k_cross_species_BROAD_drug_top50_negative_annotated.tsv"
  )
}

## -------------------------------------------------------------------------
## 2. تلاش اول: فایل‌هایی که اسمشان siginfo دارد
## -------------------------------------------------------------------------

message("Searching CLUE_raw for 'siginfo' metadata files ...")

siginfo_files <- list.files(
  clue_root,
  pattern = "siginfo.*\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

candidate_meta <- NULL
meta_source <- NULL

if (length(siginfo_files) > 0) {
  message("  - Found ", length(siginfo_files), " siginfo file(s). Trying the first one ...")
  sigfile <- siginfo_files[1]
  message("    Using: ", sigfile)
  meta_try <- try(
    readr::read_tsv(sigfile, n_max = 0, show_col_types = FALSE),
    silent = TRUE
  )
  if (!inherits(meta_try, "try-error")) {
    if ("pert_iname" %in% names(meta_try) &&
        any(c("target", "targets", "target_gene", "target_genes", "primary_target") %in% names(meta_try))) {
      candidate_meta <- sigfile
      meta_source <- "siginfo"
    } else {
      message("    First siginfo file does not contain both 'pert_iname' and any target-like column.")
    }
  } else {
    message("    Could not read header of siginfo file: ", sigfile)
  }
}

## -------------------------------------------------------------------------
## 3. اگر siginfo مناسب پیدا نشد، کل txt/tsv ها را برای target/pert_iname می‌گردیم
## -------------------------------------------------------------------------

if (is.null(candidate_meta)) {
  message("Falling back to scanning all *.txt / *.tsv for target + pert_iname ...")
  
  all_text_files <- list.files(
    clue_root,
    pattern = "\\.(txt|tsv)$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(all_text_files) == 0) {
    stop(
      "No .txt or .tsv files found under CLUE_raw. Cannot search for metadata with targets."
    )
  }
  
  find_meta_file <- function(f) {
    hdr <- try(
      readr::read_tsv(f, n_max = 0, show_col_types = FALSE),
      silent = TRUE
    )
    if (inherits(hdr, "try-error")) {
      return(NULL)
    }
    has_pert <- "pert_iname" %in% names(hdr)
    has_target <- any(c("target", "targets", "target_gene", "target_genes", "primary_target") %in% names(hdr))
    if (has_pert && has_target) {
      tibble(file = f)
    } else {
      NULL
    }
  }
  
  meta_hits <- purrr::map_dfr(all_text_files, find_meta_file)
  
  if (nrow(meta_hits) == 0) {
    stop(
      "No metadata file with both 'pert_iname' and any target-like column found under CLUE_raw.\n",
      "Either CLUE download does not include gene-level targets, or the format is very different."
    )
  }
  
  ## اگر چندتا بود، اولویت بده به فایل‌هایی که siginfo یا pert_info در نامشان هست
  meta_hits <- meta_hits %>%
    dplyr::mutate(
      weight = dplyr::case_when(
        stringr::str_detect(basename(file), "siginfo") ~ 3L,
        stringr::str_detect(basename(file), "pert") ~ 2L,
        TRUE ~ 1L
      )
    ) %>%
    dplyr::arrange(dplyr::desc(weight))
  
  candidate_meta <- meta_hits$file[1]
  meta_source <- "scan_all"
  message("  - Using metadata file: ", candidate_meta)
}

if (is.null(candidate_meta)) {
  stop("Failed to identify any CLUE metadata file with pert_iname + target columns.")
}

## -------------------------------------------------------------------------
## 4. خواندن متادیتا و استخراج ستون‌های pert_iname + target
## -------------------------------------------------------------------------

message("Reading candidate metadata file ...")

meta_hdr <- readr::read_tsv(candidate_meta, n_max = 0, show_col_types = FALSE)
meta <- readr::read_tsv(candidate_meta, show_col_types = FALSE)

target_col_candidates <- c("targets", "target", "target_gene", "target_genes", "primary_target")
target_col <- target_col_candidates[target_col_candidates %in% names(meta_hdr)][1]

if (is.na(target_col)) {
  stop(
    "Selected metadata file does not actually contain a usable target column.\n",
    "Checked for: ", paste(target_col_candidates, collapse = ", ")
  )
}

if (!("pert_iname" %in% names(meta_hdr))) {
  stop("Selected metadata file does not contain 'pert_iname' column.")
}

message("  - Using column 'pert_iname' and target column '", target_col, "' from metadata.")
message("  - Metadata source: ", meta_source)

meta_small <- meta %>%
  dplyr::select(pert_iname, !!target_col) %>%
  dplyr::rename(meta_targets_raw = !!target_col) %>%
  dplyr::mutate(
    meta_targets_raw = dplyr::if_else(
      is.na(meta_targets_raw),
      "",
      as.character(meta_targets_raw)
    ),
    meta_targets_raw = stringr::str_trim(meta_targets_raw)
  )

## -------------------------------------------------------------------------
## 5. خواندن annotated CMap top50 و merge تارگت‌ها
## -------------------------------------------------------------------------

message("Reading annotated CMap top50 ...")

drugs_ann <- readr::read_tsv(path_drugs_annot, show_col_types = FALSE)

## تشخیص ستون نام دارو در annotated
drug_name_candidates <- c("drug_name", "pert_name", "pert_iname", "compound", "name")
drug_name_col <- drug_name_candidates[drug_name_candidates %in% names(drugs_ann)][1]

if (is.na(drug_name_col)) {
  stop(
    "Could not detect a drug-name column in the annotated CMap table.\n",
    "Checked for: ", paste(drug_name_candidates, collapse = ", ")
  )
}

message("  - Using '", drug_name_col, "' as drug name column in annotated table.")

## اگر ستون targets وجود ندارد، بسازیمش
if (!("targets" %in% names(drugs_ann))) {
  drugs_ann$targets <- NA_character_
}

## join روی pert_iname vs drug_name_col
drugs_joined <- drugs_ann %>%
  dplyr::mutate(drug_name_tmp = .data[[drug_name_col]]) %>%
  dplyr::left_join(
    meta_small,
    by = c("drug_name_tmp" = "pert_iname")
  )

## اگر بعضی داروها match نشدند، مشکلی نیست؛ فقط بعضی target ندارند
matched_n <- sum(!is.na(drugs_joined$meta_targets_raw) & drugs_joined$meta_targets_raw != "")
message("  - Drugs with non-empty targets from metadata: ", matched_n, " / ", nrow(drugs_joined))

if (matched_n == 0) {
  stop(
    "Metadata file was read, but no non-empty targets matched your annotated CMap top50 by pert_iname.\n",
    "This might mean CLUE metadata doesn't actually have gene-level targets for these perturbagens."
  )
}

## حالا ستون targets را با meta_targets_raw پر می‌کنیم (فقط جایی که قبلاً NA یا خالی بوده)
drugs_updated <- drugs_joined %>%
  dplyr::mutate(
    targets = dplyr::case_when(
      !is.na(meta_targets_raw) & meta_targets_raw != "" ~ meta_targets_raw,
      TRUE ~ targets
    )
  ) %>%
  dplyr::select(-drug_name_tmp)

## کمی گزارش
filled_n <- sum(!is.na(drugs_updated$targets) & drugs_updated$targets != "")
message("  - Drugs with non-empty 'targets' after backfill: ", filled_n, " / ", nrow(drugs_updated))

if (filled_n == 0) {
  stop(
    "After attempting to backfill, 'targets' is still empty for all drugs.\n",
    "Either CLUE truly provided no targets, or the metadata format is very unusual."
  )
}

## -------------------------------------------------------------------------
## 6. نوشتن فایل annotated آپدیت‌شده
## -------------------------------------------------------------------------

readr::write_tsv(drugs_updated, path_drugs_annot)
message("Updated annotated CMap file (with 'targets') written to:\n  ", path_drugs_annot)

message("You can now re-run 23_phase6_map_drugs_to_PPI.R to map drugs onto the PPI using these targets.")
message("=== Phase 7C completed successfully. ===")
