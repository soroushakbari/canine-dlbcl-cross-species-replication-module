## 27_phase7_figures_pathway_axes_summary.R
## هدف:
##  - خلاصه‌کردن Reactome enrichment در چند محور بیولوژیک (Axis)
##  - مقایسه cross-species BROAD vs STRICT
##  - خروجی: یک شکل 3B + یک TSV summary برای محورهای بیولوژیک

message("=== Phase 7 (Figures): Reactome axis summary for BROAD & STRICT (Fig. 3B) ===")

required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "ggplot2", "forcats")

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
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)

message("Project root: ", project_root)

net_dir <- file.path(project_root, "results", "tables", "network")
fig_dir <- file.path(project_root, "results", "figures")

if (!dir.exists(net_dir)) {
  stop("Network / enrichment tables directory not found at:\n  ", net_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

path_broad_react <- file.path(
  net_dir,
  "pathway_enrichment_cross_species_BROAD_all_Reactome.tsv"
)

path_strict_react <- file.path(
  net_dir,
  "pathway_enrichment_cross_species_STRICT_all_Reactome.tsv"
)

if (!file.exists(path_broad_react)) {
  stop("Reactome enrichment file for BROAD not found at:\n  ", path_broad_react)
}
if (!file.exists(path_strict_react)) {
  stop("Reactome enrichment file for STRICT not found at:\n  ", path_strict_react)
}

## --------------------------------------------------------------------
## 2. helper: read + clean + axis mapping
## --------------------------------------------------------------------

# تابع دسته‌بندی pathway به محورهای بیولوژیک
categorize_axis <- function(desc) {
  d <- tolower(trimws(desc %||% ""))
  
  if (d == "") {
    return("Other")
  }
  
  # Glutathione / redox
  if (grepl("glutathione", d, fixed = TRUE) ||
      grepl("redox", d)) {
    return("Glutathione / redox")
  }
  
  # DNA damage / repair / PCNA / PARP / NER / BER
  if (grepl("dna damage", d) ||
      grepl("ap site", d) ||
      grepl("excision repair", d) ||
      grepl("ner ", d) ||
      grepl("base excision", d) ||
      grepl("pcna", d) ||
      grepl("parp", d)) {
    return("DNA damage / repair")
  }
  
  # Replication / cell cycle / checkpoint / mitosis / spindle / chromatid / cohesion / kinetochore
  if (grepl("dna replication", d) ||
      grepl("strand elongation", d) ||
      grepl("pre-replicative", d) ||
      grepl("synthesis of dna", d) ||
      grepl("cell cycle", d) ||
      grepl("checkpoint", d) ||
      grepl("g1/s", d) ||
      grepl("g0 and early g1", d) ||
      grepl("mitotic", d) ||
      grepl("spindle", d) ||
      grepl("chromatid", d) ||
      grepl("cohesion", d) ||
      grepl("kinetochore", d) ||
      grepl("anaphase", d) ||
      grepl("metaphase", d)) {
    return("Cell cycle / replication")
  }
  
  # RHO / cytoskeleton / motor proteins / kinesin / dynein / formin
  if (grepl("rho gtpase", d) ||
      grepl("formin", d) ||
      grepl("motor protein", d) ||
      grepl("kinesin", d) ||
      grepl("dynein", d) ||
      grepl("cytoskeleton", d)) {
    return("Cytoskeleton / RHO / motor")
  }
  
  # اگر هیچ‌کدوم نخورد
  return("Other")
}

# برای پشتیبانی از %||% (مثل rlang) بدون وابستگی
`%||%` <- function(a, b) if (!is.null(a)) a else b

prepare_reactome_with_axis <- function(path, module_name = "BROAD",
                                       padj_cutoff = 0.05) {
  message("Reading Reactome enrichment for ", module_name, " ...")
  df <- readr::read_tsv(path, show_col_types = FALSE)
  
  needed_cols <- c("ID", "Description", "pvalue", "p.adjust",
                   "qvalue", "Count")
  missing <- setdiff(needed_cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Reactome table for ", module_name, " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  df_clean <- df %>%
    dplyr::mutate(
      Module       = module_name,
      Description  = stringr::str_squish(Description),
      negLog10Padj = -log10(p.adjust)
    )
  
  # محدود به pathwayهای معنی‌دار (اگر چنین مسیرهایی وجود داشته باشند)
  df_sig <- df_clean %>%
    dplyr::filter(!is.na(p.adjust))
  
  n_sig <- df_sig %>% dplyr::filter(p.adjust <= padj_cutoff) %>% nrow()
  
  if (n_sig == 0) {
    message("  [", module_name, "] No pathways with p.adjust <= ", padj_cutoff,
            "; using all Reactome terms for axis summary.")
    df_use <- df_sig
  } else {
    df_use <- df_sig %>% dplyr::filter(p.adjust <= padj_cutoff)
  }
  
  if (nrow(df_use) == 0) {
    stop("[", module_name, "] Reactome table is empty after filtering.")
  }
  
  df_use <- df_use %>%
    dplyr::mutate(
      Axis = vapply(Description, categorize_axis, FUN.VALUE = character(1))
    )
  
  # مرتب‌سازی بر اساس p.adjust
  df_use <- df_use %>% dplyr::arrange(p.adjust)
  
  df_use
}

## --------------------------------------------------------------------
## 3. data for BROAD & STRICT
## --------------------------------------------------------------------
broad_df  <- prepare_reactome_with_axis(path_broad_react,  module_name = "BROAD",  padj_cutoff = 0.05)
strict_df <- prepare_reactome_with_axis(path_strict_react, module_name = "STRICT", padj_cutoff = 0.05)

combined <- dplyr::bind_rows(broad_df, strict_df)

## --------------------------------------------------------------------
## 4. axis-level summary
## --------------------------------------------------------------------
axis_summary <- combined %>%
  dplyr::group_by(Module, Axis) %>%
  dplyr::summarise(
    n_pathways      = dplyr::n(),
    max_negLog10Padj = max(negLog10Padj, na.rm = TRUE),
    sum_Count       = sum(Count, na.rm = TRUE),
    .groups         = "drop"
  )

# برای اطمینان، همه‌ی محورهای ممکن را در هر ماژول داشته باشیم (حتی اگر 0)
all_axes <- c(
  "Cell cycle / replication",
  "DNA damage / repair",
  "Glutathione / redox",
  "Cytoskeleton / RHO / motor",
  "Other"
)

axis_summary_full <- tidyr::crossing(
  Module = factor(c("BROAD", "STRICT"), levels = c("BROAD", "STRICT")),
  Axis   = factor(all_axes, levels = all_axes)
) %>%
  dplyr::left_join(axis_summary, by = c("Module", "Axis")) %>%
  dplyr::mutate(
    n_pathways       = tidyr::replace_na(n_pathways, 0L),
    max_negLog10Padj = tidyr::replace_na(max_negLog10Padj, 0),
    sum_Count        = tidyr::replace_na(sum_Count, 0)
  )

# ذخیره‌ی خلاصه‌ی عددی برای استفاده در متن/جدول
out_tsv <- file.path(
  net_dir,
  "Reactome_axis_summary_cross_species_BROAD_STRICT.tsv"
)
readr::write_tsv(axis_summary_full, out_tsv)
message("Saved Reactome axis summary table to:")
message("  ", out_tsv)

## --------------------------------------------------------------------
## 5. Fig. 3B – tile/bubble axis plot
## --------------------------------------------------------------------
message("Building Fig. 3B axis summary plot ...")

axis_summary_for_plot <- axis_summary_full %>%
  dplyr::mutate(
    Module = factor(Module, levels = c("BROAD", "STRICT")),
    Axis   = factor(Axis,   levels = all_axes)
  )

p_axis <- ggplot(axis_summary_for_plot,
                 aes(x = Axis, y = Module)) +
  geom_tile(aes(fill = max_negLog10Padj), color = "grey85") +
  geom_point(aes(size = sum_Count), shape = 21, fill = "white", color = "black") +
  geom_text(
    aes(label = ifelse(n_pathways > 0, n_pathways, "")),
    size = 3
  ) +
  scale_fill_gradient(
    name = expression(max~"-log"[10]~"adj. p"),
    low  = "#f7fbff",
    high = "#08306b",
    limits = c(0, max(axis_summary_for_plot$max_negLog10Padj, na.rm = TRUE) * 1.05)
  ) +
  scale_size_continuous(
    name  = "Total gene hits",
    range = c(2.5, 8)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Biological axes of the cross-species DLBCL modules (Reactome summary)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 25, hjust = 1, vjust = 1, size = 9),
    axis.text.y      = element_text(size = 10, face = "bold"),
    panel.grid       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(hjust = 0, face = "bold", size = 13),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )

out_png <- file.path(fig_dir, "Fig3B_Reactome_axis_summary_BROAD_STRICT.png")
out_pdf <- file.path(fig_dir, "Fig3B_Reactome_axis_summary_BROAD_STRICT.pdf")

ggsave(out_png, p_axis, width = 8.5, height = 4.5, dpi = 400)
ggsave(out_pdf, p_axis, width = 8.5, height = 4.5)

message("Saved Fig. 3B axis summary plot to:")
message("  ", out_png)
message("  ", out_pdf)

message("=== Phase 7 (Figures) – Fig. 3B completed successfully. ===")
