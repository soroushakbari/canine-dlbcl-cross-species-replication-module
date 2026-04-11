## 28_phase8_figures_PPI_centrality_and_drugs.R
## هدف:
##  - Fig. 4A: PPI centrality landscape for cross-species BROAD module
##  - Fig. 4B: CMap drug–target landscape with PPI mapping

message("=== Phase 8 (Figures): PPI centrality + CMap drug landscape (Fig. 4A, 4B) ===")

required_pkgs <- c("readr", "dplyr", "tibble", "stringr",
                   "ggplot2", "forcats", "ggrepel")

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
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(ggrepel)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)

message("Project root: ", project_root)

net_dir  <- file.path(project_root, "results", "tables", "network")
drug_dir <- file.path(project_root, "results", "tables", "Drug")
fig_dir  <- file.path(project_root, "results", "figures")

if (!dir.exists(net_dir))  stop("Network tables dir not found at:\n  ", net_dir)
if (!dir.exists(drug_dir)) stop("Drug tables dir not found at:\n  ", drug_dir)
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

nodes_path      <- file.path(net_dir,  "PPI_cross_species_BROAD_nodes.tsv")
hubs_drug_path  <- file.path(net_dir,  "PPI_cross_species_BROAD_hubs_with_drug_hits.tsv")
drug_map_path   <- file.path(drug_dir, "CMap_queryl1k_cross_species_BROAD_drug_top50_with_PPI_mapping.tsv")

if (!file.exists(nodes_path)) {
  stop("PPI nodes table not found at:\n  ", nodes_path)
}
if (!file.exists(hubs_drug_path)) {
  stop("Hub+drug table not found at:\n  ", hubs_drug_path)
}
if (!file.exists(drug_map_path)) {
  stop("Drug–PPI mapping table not found at:\n  ", drug_map_path)
}

## --------------------------------------------------------------------
## 2. Fig. 4A – PPI centrality landscape
## --------------------------------------------------------------------
message("Reading PPI node and hub+drug tables ...")

nodes <- readr::read_tsv(nodes_path, show_col_types = FALSE)
hubs_drug <- readr::read_tsv(hubs_drug_path, show_col_types = FALSE)

needed_nodes <- c("gene", "degree", "betweenness")
missing_nodes <- setdiff(needed_nodes, names(nodes))
if (length(missing_nodes) > 0) {
  stop("PPI nodes table is missing required columns: ",
       paste(missing_nodes, collapse = ", "))
}

# اگر ستون‌های اضافی نباشند، مقداردهی ایمن
if (!"hub_type" %in% names(nodes)) {
  nodes$hub_type <- "non_hub"
}
if (!"logFC_human" %in% names(nodes)) {
  nodes$logFC_human <- NA_real_
}

# hub_type را factor استاندارد می‌کنیم
nodes <- nodes %>%
  mutate(
    hub_type = case_when(
      hub_type %in% c("hub_bottleneck", "hub-bottleneck") ~ "hub_bottleneck",
      hub_type %in% c("bottleneck")                      ~ "bottleneck",
      hub_type %in% c("hub")                             ~ "hub",
      TRUE                                               ~ "non_hub"
    ),
    hub_type = factor(
      hub_type,
      levels = c("non_hub", "hub", "bottleneck", "hub_bottleneck")
    )
  )

# join با جدول hub+drug برای آوردن n_drugs و لیست داروها
nodes2 <- nodes %>%
  left_join(
    hubs_drug %>%
      select(gene, n_drugs, drugs),
    by = "gene"
  ) %>%
  mutate(
    n_drugs = replace_na(n_drugs, 0L),
    has_drug = if_else(n_drugs > 0, "Drug-targeted", "No drug hit"),
    has_drug = factor(has_drug, levels = c("No drug hit", "Drug-targeted"))
  )

# برای readability، فقط نودهای با degree بالا یا hub/bottleneck را نگه می‌داریم
deg_q90 <- quantile(nodes2$degree, 0.90, na.rm = TRUE)

nodes_plot <- nodes2 %>%
  filter(degree >= deg_q90 | hub_type != "non_hub")

message("  - Nodes used for Fig. 4A: ", nrow(nodes_plot),
        " (degree >= 90th percentile or hub/bottleneck)")

# نودهایی که دارو می‌زنند (TOP2A, PARP1) را جدا می‌کنیم
nodes_drug <- nodes_plot %>%
  filter(has_drug == "Drug-targeted")

# رنگ‌ها برای hub_type
hub_cols <- c(
  "non_hub"       = "grey80",
  "hub"           = "#3182bd",
  "bottleneck"    = "#e6550d",
  "hub_bottleneck"= "#9e0168"
)

p4A <- ggplot(nodes_plot, aes(x = degree, y = betweenness)) +
  # پس‌زمینه: همه‌ی این ساب‌ست
  geom_point(aes(color = hub_type),
             size = 2.2,
             alpha = 0.9) +
  scale_color_manual(
    name   = "Node class",
    values = hub_cols,
    breaks = c("hub_bottleneck", "hub", "bottleneck", "non_hub"),
    labels = c("Hub–bottleneck", "Hub", "Bottleneck", "Other")
  ) +
  # حلقه‌ی ضخیم دور نودهایی که دارو می‌زنند
  geom_point(
    data = nodes_drug,
    aes(x = degree, y = betweenness),
    shape = 21,
    size  = 4.5,
    stroke = 1.2,
    colour = "black",
    fill   = NA
  ) +
  # برچسب برای TOP2A و PARP1 (اگر حضور دارند)
  ggrepel::geom_text_repel(
    data = nodes_plot %>%
      filter(gene %in% c("TOP2A", "PARP1")),
    aes(label = gene),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.25,
    point.padding = 0.25,
    max.overlaps = 20,
    min.segment.length = 0
  ) +
  labs(
    x = "Degree centrality",
    y = "Betweenness centrality",
    title = "PPI centrality landscape of the cross-species BROAD module"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(hjust = 0, face = "bold", size = 13)
  )

out4A_png <- file.path(fig_dir, "Fig4A_PPI_centrality_BROAD.png")
out4A_pdf <- file.path(fig_dir, "Fig4A_PPI_centrality_BROAD.pdf")

ggsave(out4A_png, p4A, width = 7.5, height = 5.5, dpi = 400)
ggsave(out4A_pdf, p4A, width = 7.5, height = 5.5)

message("Saved Fig. 4A centrality plot to:")
message("  ", out4A_png)
message("  ", out4A_pdf)

## --------------------------------------------------------------------
## 3. Fig. 4B – CMap drug landscape with PPI mapping
## --------------------------------------------------------------------
message("Reading CMap drug–PPI mapping table ...")

drug_map <- readr::read_tsv(drug_map_path, show_col_types = FALSE)

needed_drug <- c("pert_name", "drug_name", "min_score", "MOA_category",
                 "n_targets_in_module", "n_targets_in_hubs",
                 "n_targets_in_hub_bottleneck",
                 "genes_in_module", "genes_in_hubs", "genes_in_hub_bottleneck")
missing_drug <- setdiff(needed_drug, names(drug_map))
if (length(missing_drug) > 0) {
  stop("Drug mapping table is missing required columns: ",
       paste(missing_drug, collapse = ", "))
}

drug_df <- drug_map %>%
  filter(!is.na(n_targets_in_module),
         n_targets_in_module > 0) %>%
  mutate(
    neg_conn = -min_score,
    # core hits: TOP2A و/یا PARP1
    genes_in_hubs          = genes_in_hubs %||% "",
    genes_in_hub_bottleneck = genes_in_hub_bottleneck %||% "",
    hit_TOP2A = str_detect(genes_in_hubs, "TOP2A") |
      str_detect(genes_in_hub_bottleneck, "TOP2A"),
    hit_PARP1 = str_detect(genes_in_hubs, "PARP1") |
      str_detect(genes_in_hub_bottleneck, "PARP1"),
    core_hit = case_when(
      hit_TOP2A & hit_PARP1 ~ "TOP2A + PARP1",
      hit_TOP2A            ~ "TOP2A",
      hit_PARP1            ~ "PARP1",
      TRUE                 ~ "None"
    ),
    core_hit = factor(core_hit,
                      levels = c("None", "TOP2A", "PARP1", "TOP2A + PARP1")),
    drug_name_plot = if_else(
      is.na(drug_name) | drug_name == "",
      pert_name,
      drug_name
    )
  ) %>%
  arrange(neg_conn) %>%
  mutate(
    drug_name_plot = fct_reorder(drug_name_plot, neg_conn)
  )

message("  - Drugs with at least one target in module: ", nrow(drug_df))

# شکل لالی‌پاپ افقی
p4B <- ggplot(drug_df,
              aes(x = neg_conn, y = drug_name_plot)) +
  geom_segment(
    aes(x = 0, xend = neg_conn, y = drug_name_plot, yend = drug_name_plot),
    colour = "grey85",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      fill  = MOA_category,
      size  = n_targets_in_module,
      shape = core_hit
    ),
    colour = "black",
    alpha = 0.95
  ) +
  scale_size_continuous(
    name  = "Targets in module",
    range = c(2.5, 7)
  ) +
  scale_shape_manual(
    name  = "Core hub overlap",
    values = c(
      "None"          = 21,
      "TOP2A"         = 23,
      "PARP1"         = 24,
      "TOP2A + PARP1" = 22
    )
  ) +
  guides(
    fill = guide_legend(title = "MOA category", override.aes = list(shape = 21, size = 4))
  ) +
  labs(
    x = expression("- connectivity strength (min score)"),
    y = NULL,
    title = "CMap drug landscape with cross-species module and PPI overlap"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 8),
    legend.position  = "right",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(hjust = 0, face = "bold", size = 13),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )

out4B_png <- file.path(fig_dir, "Fig4B_CMap_drug_landscape_PPI_overlap.png")
out4B_pdf <- file.path(fig_dir, "Fig4B_CMap_drug_landscape_PPI_overlap.pdf")

ggsave(out4B_png, p4B, width = 7.8, height = 6.2, dpi = 400)
ggsave(out4B_pdf, p4B, width = 7.8, height = 6.2)

message("Saved Fig. 4B drug landscape plot to:")
message("  ", out4B_png)
message("  ", out4B_pdf)

message("=== Phase 8 (Figures) completed successfully. ===")
