## 28_phase8_figures_PPI_centrality_and_drugs.R
## هدف:
##  - Fig. 3A: PPI centrality landscape for cross-species BROAD module
##  - Fig. 3C: CMap drug–target landscape with PPI mapping
## نکات ریوایز:
##  - all hub–bottlenecks are labelled in Fig. 3A
##  - PARP1 is additionally labelled as a drug-anchored bottleneck
##  - full list of hub–bottleneck and bottleneck genes is exported for Supplementary Table S1

message("=== Phase 8 (Figures): PPI centrality + CMap drug landscape (Fig. 3A, 3C) ===")

required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "ggplot2", "forcats", "ggrepel", "tidyr"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) and rerun the script."
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
  library(tidyr)
})

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath(
  "D:/Research/My Articles/DLBCL drug",
  winslash = "/",
  mustWork = TRUE
)

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

for (p in c(nodes_path, hubs_drug_path, drug_map_path)) {
  if (!file.exists(p)) {
    stop("Required file not found at:\n  ", p)
  }
}

## --------------------------------------------------------------------
## 2. helpers
## --------------------------------------------------------------------
normalize_hub_type <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("hub_bottleneck", "hub-bottleneck", "hub bottleneck") ~ "hub_bottleneck",
    x %in% c("bottleneck") ~ "bottleneck",
    x %in% c("hub") ~ "hub",
    TRUE ~ "non_hub"
  )
}

## --------------------------------------------------------------------
## 3. Fig. 3A – PPI centrality landscape
## --------------------------------------------------------------------
message("Reading PPI node and hub+drug tables ...")

nodes <- readr::read_tsv(nodes_path, show_col_types = FALSE)
hubs_drug <- readr::read_tsv(hubs_drug_path, show_col_types = FALSE)

needed_nodes <- c("gene", "degree", "betweenness")
missing_nodes <- setdiff(needed_nodes, names(nodes))
if (length(missing_nodes) > 0) {
  stop(
    "PPI nodes table is missing required columns: ",
    paste(missing_nodes, collapse = ", ")
  )
}

if (!"hub_type" %in% names(nodes)) {
  nodes$hub_type <- "non_hub"
}
if (!"logFC_human" %in% names(nodes)) {
  nodes$logFC_human <- NA_real_
}

nodes <- nodes %>%
  mutate(
    gene = as.character(gene),
    degree = as.numeric(degree),
    betweenness = as.numeric(betweenness),
    hub_type = normalize_hub_type(hub_type),
    hub_type = factor(
      hub_type,
      levels = c("non_hub", "hub", "bottleneck", "hub_bottleneck")
    )
  )

hubs_drug2 <- hubs_drug %>%
  mutate(gene = as.character(gene)) %>%
  select(any_of(c("gene", "n_drugs", "drugs")))

nodes2 <- nodes %>%
  left_join(hubs_drug2, by = "gene") %>%
  mutate(
    n_drugs = tidyr::replace_na(as.integer(n_drugs), 0L),
    drugs   = tidyr::replace_na(drugs, ""),
    has_drug = n_drugs > 0
  )
## export full central-node list for Supplementary Table S1
central_nodes_tbl <- nodes2 %>%
  filter(hub_type %in% c("hub_bottleneck", "bottleneck")) %>%
  mutate(
    node_class = dplyr::case_when(
      hub_type == "hub_bottleneck" ~ "Hub–bottleneck",
      hub_type == "bottleneck" ~ "Bottleneck",
      TRUE ~ "Other"
    ),
    labelled_in_Fig3A = gene %in% c(
      nodes2 %>% filter(hub_type == "hub_bottleneck") %>% pull(gene),
      "PARP1"
    )
  ) %>%
  arrange(
    factor(node_class, levels = c("Hub–bottleneck", "Bottleneck")),
    desc(betweenness),
    desc(degree),
    gene
  ) %>%
  select(gene, degree, betweenness, node_class, n_drugs, drugs, labelled_in_Fig3A)

central_nodes_out <- file.path(net_dir, "PPI_central_nodes_for_STable_S1.tsv")
readr::write_tsv(central_nodes_tbl, central_nodes_out)

message("Saved central-node table for Supplementary Table S1 to:")
message("  ", central_nodes_out)

## for readability: retain degree >= 90th percentile OR any central node
deg_q90 <- stats::quantile(nodes2$degree, 0.90, na.rm = TRUE)

nodes_plot <- nodes2 %>%
  filter(degree >= deg_q90 | hub_type != "non_hub")

message(
  "  - Nodes used for Fig. 3A: ", nrow(nodes_plot),
  " (degree >= 90th percentile or hub/bottleneck)"
)

## nodes with direct drug hits: black ring
nodes_drug <- nodes_plot %>%
  filter(has_drug)

## labels:
##  - all hub–bottlenecks
##  - PARP1 additionally, even though it is bottleneck not hub–bottleneck
label_df <- nodes_plot %>%
  filter(hub_type == "hub_bottleneck" | gene == "PARP1") %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(betweenness), desc(degree), gene)

message("Genes labelled in Fig. 3A:")
cat("  ", paste(label_df$gene, collapse = ", "), "\n")

hub_cols <- c(
  "non_hub"        = "grey80",
  "hub"            = "#3182bd",
  "bottleneck"     = "#e6550d",
  "hub_bottleneck" = "#9e0168"
)

p3A <- ggplot(nodes_plot, aes(x = degree, y = betweenness)) +
  geom_point(
    aes(color = hub_type),
    size = 2.2,
    alpha = 0.9
  ) +
  scale_color_manual(
    name   = "Node class",
    values = hub_cols,
    breaks = c("hub_bottleneck", "hub", "bottleneck", "non_hub"),
    labels = c("Hub–bottleneck", "Hub", "Bottleneck", "Other")
  ) +
  ## black ring around drug-targeted nodes
  geom_point(
    data = nodes_drug,
    aes(x = degree, y = betweenness),
    inherit.aes = FALSE,
    shape  = 21,
    size   = 4.5,
    stroke = 1.2,
    colour = "black",
    fill   = NA
  ) +
  ## labels for all hub–bottlenecks + PARP1
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.3,
    fontface = "bold",
    box.padding = 0.28,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 12345,
    segment.color = "grey40"
  ) +
  labs(
    x = "Degree centrality",
    y = "Betweenness centrality"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9)
  )

out3A_png <- file.path(fig_dir, "Fig3A_PPI_centrality_BROAD.png")
out3A_pdf <- file.path(fig_dir, "Fig3A_PPI_centrality_BROAD.pdf")

ggsave(out3A_png, p3A, width = 7.5, height = 5.5, dpi = 400)
ggsave(out3A_pdf, p3A, width = 7.5, height = 5.5)

message("Saved Fig. 3A centrality plot to:")
message("  ", out3A_png)
message("  ", out3A_pdf)

## --------------------------------------------------------------------
## 4. Fig. 3C – CMap drug landscape with PPI mapping
## --------------------------------------------------------------------
message("Reading CMap drug–PPI mapping table ...")

drug_map <- readr::read_tsv(drug_map_path, show_col_types = FALSE)

needed_drug <- c(
  "pert_name", "drug_name", "min_score", "MOA_category",
  "n_targets_in_module", "n_targets_in_hubs",
  "n_targets_in_hub_bottleneck",
  "genes_in_module", "genes_in_hubs", "genes_in_hub_bottleneck"
)
missing_drug <- setdiff(needed_drug, names(drug_map))
if (length(missing_drug) > 0) {
  stop(
    "Drug mapping table is missing required columns: ",
    paste(missing_drug, collapse = ", ")
  )
}

drug_df <- drug_map %>%
  mutate(
    genes_in_hubs = tidyr::replace_na(genes_in_hubs, ""),
    genes_in_hub_bottleneck = tidyr::replace_na(genes_in_hub_bottleneck, ""),
    drug_name_plot = dplyr::if_else(
      is.na(drug_name) | drug_name == "",
      pert_name,
      drug_name
    )
  ) %>%
  filter(!is.na(n_targets_in_module), n_targets_in_module > 0) %>%
  mutate(
    neg_conn = -min_score,
    hit_TOP2A = str_detect(genes_in_hubs, "TOP2A") |
      str_detect(genes_in_hub_bottleneck, "TOP2A"),
    hit_PARP1 = str_detect(genes_in_hubs, "PARP1") |
      str_detect(genes_in_hub_bottleneck, "PARP1"),
    core_hit = case_when(
      hit_TOP2A & hit_PARP1 ~ "TOP2A + PARP1",
      hit_TOP2A             ~ "TOP2A",
      hit_PARP1             ~ "PARP1",
      TRUE                  ~ "None"
    ),
    core_hit = factor(
      core_hit,
      levels = c("None", "TOP2A", "PARP1", "TOP2A + PARP1")
    )
  ) %>%
  arrange(neg_conn) %>%
  mutate(
    drug_name_plot = forcats::fct_reorder(drug_name_plot, neg_conn)
  )

message("  - Drugs with at least one target in module: ", nrow(drug_df))

p3C <- ggplot(drug_df, aes(x = neg_conn, y = drug_name_plot)) +
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
    name   = "Core hub overlap",
    values = c(
      "None"          = 21,
      "TOP2A"         = 23,
      "PARP1"         = 24,
      "TOP2A + PARP1" = 22
    )
  ) +
  guides(
    fill = guide_legend(
      title = "MOA category",
      override.aes = list(shape = 21, size = 4)
    )
  ) +
  labs(
    x = expression("- connectivity strength (min score)"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 8),
    axis.text.x      = element_text(size = 8),
    legend.position  = "right",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )

out3C_png <- file.path(fig_dir, "Fig3C_CMap_drug_landscape_PPI_overlap.png")
out3C_pdf <- file.path(fig_dir, "Fig3C_CMap_drug_landscape_PPI_overlap.pdf")

ggsave(out3C_png, p3C, width = 7.8, height = 6.2, dpi = 400)
ggsave(out3C_pdf, p3C, width = 7.8, height = 6.2)

message("Saved Fig. 3C drug landscape plot to:")
message("  ", out3C_png)
message("  ", out3C_pdf)

message("=== Phase 8 (Figures) completed successfully. ===")
