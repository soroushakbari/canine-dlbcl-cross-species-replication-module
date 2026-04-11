## 35_phase9D_crossCohort_survival_forest.R
## Fig 6E – Cross-cohort survival association of the cross-species DLBCL module

message("=== Phase 9D: Cross-cohort survival forest plot (Fig 6E) ===")

required_pkgs <- c("tibble", "dplyr", "ggplot2")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) و دوباره اسکریپت را اجرا کن."
    )
  }
}

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(ggplot2)
})

## ------------------------------------------------------------------
## 1. paths
## ------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

fig_dir <- file.path(project_root, "results", "figures")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

## ------------------------------------------------------------------
## 2. Cross-cohort summary (numbers از اسکریپت‌های بقا)
##    TCGA:   12_phase4E_TCGA_survival_analysis.R
##    GSE130874: 34_phase8A_GSE130874_PFS_analysis.R
##    GSE31312:  39_phase9C_GSE31312_OS_analysis_alt_from_pdf.R
## ------------------------------------------------------------------

forest_df <- tibble::tribble(
  ~cohort_label,                                             ~HR,   ~lower95, ~upper95, ~n,   ~events, ~species, ~endpoint,
  "GSE31312 (human, OS) (n = 470, events = 170)",            0.54,  0.32,     0.90,     470L, 170L,    "human",  "OS",
  "TCGA-DLBC (human, OS) (n = 47, events = 9)",              0.584, 0.26,     1.312,    47L,  9L,      "human",  "OS",
  "GSE130874 (dog, PFS) (n = 23, events = 21)",              0.49,  0.26,     0.94,     23L,  21L,     "dog",    "PFS"
)

## ترتیب: اول human، بعد dog
forest_df <- forest_df %>%
  dplyr::mutate(
    cohort_label = factor(
      cohort_label,
      levels = c(
        "GSE31312 (human, OS) (n = 470, events = 170)",
        "TCGA-DLBC (human, OS) (n = 47, events = 9)",
        "GSE130874 (dog, PFS) (n = 23, events = 21)"
      )
    )
  )

## ------------------------------------------------------------------
## 3. Plot
## ------------------------------------------------------------------
message("Building forest plot ...")

p <- ggplot(forest_df,
            aes(x = HR, y = cohort_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95),
                 height = 0.15, linewidth = 0.6) +
  geom_point(size = 3) +
  scale_x_log10(
    name   = "Hazard ratio per 1 SD increase in module score (log scale)",
    limits = c(0.25, 2),
    breaks = c(0.25, 0.5, 1, 2),
    minor_breaks = NULL
  ) +
  labs(
    title = "Cross-cohort survival association of the cross-species DLBCL module",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
  )

out_png <- file.path(fig_dir, "Fig6E_crossCohort_survival_forest.png")
out_pdf <- file.path(fig_dir, "Fig6E_crossCohort_survival_forest.pdf")

ggsave(out_png, p, width = 7.2, height = 3.4, dpi = 400)
ggsave(out_pdf, p, width = 7.2, height = 3.4)

message("Saved forest plot to:")
message("  ", out_png)
message("  ", out_pdf)
message("=== Phase 9D completed successfully. ===")

