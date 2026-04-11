## 43_phase10D_EVmiRNA_summary_dotplot.R
## Fig 6B – EV-miRNA summary vs PPI core (GSE171272)

message("=== Phase 10D: EV-miRNA summary dotplot for PPI core (Fig 6B) ===")

required_pkgs <- c("readr", "dplyr", "tibble", "forcats", "ggplot2", "stringr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Run:\n  install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      "))"
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(forcats)
  library(ggplot2)
  library(stringr)
})

## -------------------------------------------------------------------
## 1. paths
## -------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

tab_dir <- file.path(project_root, "results", "tables", "miRNA")
fig_dir <- file.path(project_root, "results", "figures")

if (!dir.exists(tab_dir)) {
  stop("miRNA tables directory not found at:\n  ", tab_dir)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

summary_path <- file.path(
  tab_dir,
  "GSE171272_EVmiRNA_targets_module_summary_by_miRNA.tsv"
)

if (!file.exists(summary_path)) {
  stop("EV-miRNA summary table not found at:\n  ", summary_path)
}

## -------------------------------------------------------------------
## 2. read + prepare
## -------------------------------------------------------------------
message("--- Step 1: Load EV-miRNA summary table ---")
mir_sum_raw <- readr::read_tsv(summary_path, show_col_types = FALSE)

needed_cols <- c(
  "mature_mirna_id",
  "n_targets_total",
  "n_targets_BROAD",
  "n_targets_PPI_core",
  "n_hits_TOP2A",
  "n_hits_PARP1",
  "hits_any_core_hub",
  "hits_TOP2A_or_PARP1",
  "logFC",
  "adj.P.Val"
)

missing <- setdiff(needed_cols, names(mir_sum_raw))
if (length(missing) > 0) {
  stop("Summary table is missing required columns: ",
       paste(missing, collapse = ", "))
}

mir_sum <- mir_sum_raw %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      logFC >  0 ~ "Up in DLBCL EV",
      logFC <  0 ~ "Down in DLBCL EV",
      TRUE       ~ "No change"
    ),
    negLog10AdjP = -log10(adj.P.Val + 1e-300),
    any_core    = hits_any_core_hub %in% TRUE,
    hits_TOP2A_or_PARP1 = hits_TOP2A_or_PARP1 %in% TRUE
  )

message("--- Step 2: Select miRNAs for plotting ---")

## اگر می‌خواهی خودت لیست را دستی انتخاب کنی، این وکتور را پر کن:
manual_miRNAs <- c(
  # مثال:
  # "hsa-let-7b-5p", "hsa-miR-205-5p", "hsa-miR-21-5p",
  # "hsa-miR-124-3p", "hsa-miR-29a-3p", "hsa-miR-29b-3p",
  # "hsa-miR-181a-5p", "hsa-miR-147a", "hsa-miR-10a-5p"
)

top_n <- 12L

if (length(manual_miRNAs) > 0) {
  mir_plot <- mir_sum %>%
    dplyr::filter(mature_mirna_id %in% manual_miRNAs)
  message("  Using MANUAL miRNA list (n = ", nrow(mir_plot), ").")
} else {
  mir_plot <- mir_sum %>%
    dplyr::filter(any_core) %>%
    dplyr::arrange(dplyr::desc(n_targets_PPI_core)) %>%
    dplyr::slice_head(n = top_n)
  message("  Using top ", top_n,
          " miRNAs by n_targets_PPI_core (any_core == TRUE).")
}

if (nrow(mir_plot) == 0) {
  stop("No miRNAs selected for plotting – check filters/manual list.")
}

## مرتب‌سازی با توجه به تعداد تارگت‌ها در core
mir_plot <- mir_plot %>%
  dplyr::mutate(
    miRNA_label = mature_mirna_id,
    miRNA_label = stringr::str_replace(miRNA_label, "^hsa-", ""),
    miRNA_label = forcats::fct_reorder(miRNA_label,
                                       n_targets_PPI_core,
                                       .desc = TRUE)
  )

message("--- Step 3: Build dotplot (lollipop) ---")

p <- ggplot(mir_plot,
            aes(x = n_targets_PPI_core,
                y = miRNA_label)) +
  ## segment از صفر تا مقدار core targets
  geom_segment(
    aes(x = 0,
        xend = n_targets_PPI_core,
        y = miRNA_label,
        yend = miRNA_label),
    linewidth = 0.6,
    colour = "grey80"
  ) +
  ## نقطه در انتهای segment
  geom_point(
    aes(
      colour   = direction,
      size     = negLog10AdjP,
      shape    = hits_TOP2A_or_PARP1
    ),
    stroke = 0.9
  ) +
  scale_x_continuous(
    name   = "Number of targets in PPI core",
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_discrete(name = NULL) +
  scale_size_continuous(
    name  = expression(-log[10]~"adj. p (EV miRNA DE)"),
    range = c(2.5, 6)
  ) +
  scale_colour_manual(
    name   = "Direction in EV (DLBCL vs control)",
    values = c(
      "Up in DLBCL EV"   = "#b2182b",
      "Down in DLBCL EV" = "#2166ac",
      "No change"        = "grey40"
    )
  ) +
  scale_shape_manual(
    name = "Hits TOP2A or PARP1",
    values = c(`TRUE` = 21, `FALSE` = 16),
    labels = c(`TRUE` = "Yes", `FALSE` = "No")
  ) +
  guides(
    size   = guide_legend(order = 1),
    colour = guide_legend(order = 2),
    shape  = guide_legend(order = 3)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    plot.title         = element_text(
      hjust = 0,
      face  = "bold",
      size  = 13,
      margin = margin(b = 6)
    ),
    plot.margin        = margin(t = 8, r = 10, b = 8, l = 8)
  ) +
  ggtitle("EV-miRNAs converging on the cross-species PPI core")

out_png <- file.path(fig_dir,
                     "Fig6B_GSE171272_EVmiRNA_PPIcore_summary_dotplot.png")
out_pdf <- file.path(fig_dir,
                     "Fig6B_GSE171272_EVmiRNA_PPIcore_summary_dotplot.pdf")

ggsave(out_png, p, width = 7.5, height = 4.5, dpi = 400)
ggsave(out_pdf, p, width = 7.5, height = 4.5)

message("Saved Fig 6B EV-miRNA summary to:")
message("  ", out_png)
message("  ", out_pdf)
message("=== Phase 10D completed successfully. ===")
