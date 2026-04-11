## 26_phase7_figures_pathway_enrichment.R
## Reactome bubble plot (BROAD vs STRICT) – upgraded

message("=== Phase 7 (Figures): Reactome pathway bubble plots for BROAD & STRICT ===")

required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "forcats", "ggplot2")

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
  library(forcats)
  library(ggplot2)
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
## 2. helper: prepare Reactome df
## --------------------------------------------------------------------
prepare_reactome_df <- function(path, module_name = "BROAD", top_n = 15) {
  message("Reading Reactome enrichment for ", module_name, " ...")
  df <- readr::read_tsv(path, show_col_types = FALSE)
  
  needed_cols <- c("ID", "Description", "pvalue", "p.adjust",
                   "qvalue", "Count", "RichFactor")
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
    ) %>%
    dplyr::arrange(p.adjust)
  
  n_rows <- nrow(df_clean)
  top_n_final <- min(top_n, n_rows)
  
  df_top <- df_clean %>%
    dplyr::slice_head(n = top_n_final) %>%
    dplyr::mutate(
      Description_short = Description,
      Description_short = dplyr::if_else(
        nchar(Description_short) > 60,
        paste0(substr(Description_short, 1, 57), "..."),
        Description_short
      )
    )
  
  df_top
}

## --------------------------------------------------------------------
## 3. prepare BROAD & STRICT
## --------------------------------------------------------------------
broad_df  <- prepare_reactome_df(path_broad_react,  module_name = "BROAD",  top_n = 15)
strict_df <- prepare_reactome_df(path_strict_react, module_name = "STRICT", top_n = 10)

combined <- dplyr::bind_rows(broad_df, strict_df) %>%
  dplyr::mutate(
    Module = factor(Module, levels = c("BROAD", "STRICT"))
  )

## order y per module by RichFactor (or negLog10Padj – هر کدوم را می‌خواهی)
combined <- combined %>%
  dplyr::group_by(Module) %>%
  dplyr::mutate(
    Description_f = forcats::fct_reorder(Description_short, RichFactor, .desc = FALSE)
  ) %>%
  dplyr::ungroup()

## x-range و breaks تمیز
range_x  <- range(combined$RichFactor, na.rm = TRUE)
pad      <- diff(range_x) * 0.08
limits_x <- c(range_x[1] - pad, range_x[2] + pad)
breaks_x <- pretty(range_x, n = 4)

## --------------------------------------------------------------------
## 4. plot
## --------------------------------------------------------------------
message("Building bubble plot figure ...")

p <- ggplot(combined, aes(x = RichFactor, y = Description_f)) +
  geom_point(
    aes(size = Count, fill = negLog10Padj),
    shape = 21,
    color = "grey20",
    alpha = 0.95
  ) +
  scale_size_continuous(
    name  = "Gene count",
    range = c(3, 9)
  ) +
  scale_fill_gradient(
    name = expression(-log[10]~"adj. p"),
    low  = "#dce9f2",
    high = "#084594"
  ) +
  scale_x_continuous(
    name   = "Rich factor",
    limits = limits_x,
    breaks = breaks_x
  ) +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  labs(
    y = NULL,
    title = "Reactome enrichment for cross-species modules"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),
    axis.title.x       = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    plot.title         = element_text(
      hjust = 0,
      face  = "bold",
      size  = 13,
      margin = margin(b = 6)
    ),
    plot.margin        = margin(t = 8, r = 8, b = 8, l = 8)
  )

out_png <- file.path(fig_dir, "Fig3A_Reactome_enrichment_BROAD_STRICT_bubble.png")
out_pdf <- file.path(fig_dir, "Fig3A_Reactome_enrichment_BROAD_STRICT_bubble.pdf")

ggsave(out_png, p, width = 10.5, height = 5.8, dpi = 400)
ggsave(out_pdf, p, width = 10.5, height = 5.8)

message("Saved Reactome bubble plot to:")
message("  ", out_png)
message("  ", out_pdf)

message("=== Phase 7 (Figures) completed successfully. ===")
