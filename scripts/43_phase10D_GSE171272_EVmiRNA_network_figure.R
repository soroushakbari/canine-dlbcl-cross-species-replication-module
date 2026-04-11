## 43_phase10D_GSE171272_EVmiRNA_network_figure.R
## EV-miRNA -> PPI core network figure (GSE171272, TOP2A/PARP1)

message("=== Phase 10D: EV-miRNA -> PPI core network figure (GSE171272) ===")

## ----------------------------------------------------------------------
## 0. packages
## ----------------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tibble", "stringr", "ggplot2", "forcats")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      ")) and re-run this script."
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

## ----------------------------------------------------------------------
## 1. paths
## ----------------------------------------------------------------------
project_root <- normalizePath(file.path(".."),
                              winslash = "/",
                              mustWork = TRUE)
message("Project root: ", project_root)

mi_dir  <- file.path(project_root, "results", "tables", "miRNA")
net_dir <- file.path(project_root, "results", "tables", "network")
fig_dir <- file.path(project_root, "results", "figures")

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

path_targets   <- file.path(mi_dir,
                            "GSE171272_EVmiRNA_targets_in_BROAD_module_PPI.tsv"
)
path_summary   <- file.path(mi_dir,
                            "GSE171272_EVmiRNA_targets_module_summary_by_miRNA.tsv"
)
path_ppi_nodes <- file.path(net_dir, "PPI_cross_species_BROAD_nodes.tsv")

for (p in c(path_targets, path_summary, path_ppi_nodes)) {
  if (!file.exists(p)) {
    stop("Required file not found:\n  ", p)
  }
}

## ----------------------------------------------------------------------
## 2. read data
## ----------------------------------------------------------------------
targets    <- readr::read_tsv(path_targets, show_col_types = FALSE)
summary_mi <- readr::read_tsv(path_summary, show_col_types = FALSE)
ppi_nodes  <- readr::read_tsv(path_ppi_nodes, show_col_types = FALSE)

needed_tgt_cols <- c("mature_mirna_id", "target_symbol", "in_PPI_core")
missing_tgt     <- setdiff(needed_tgt_cols, names(targets))
if (length(missing_tgt) > 0) {
  stop("Targets table is missing columns: ",
       paste(missing_tgt, collapse = ", "))
}

needed_sum_cols <- c("mature_mirna_id", "n_targets_PPI_core",
                     "n_targets_BROAD", "logFC", "adj.P.Val")
missing_sum     <- setdiff(needed_sum_cols, names(summary_mi))
if (length(missing_sum) > 0) {
  stop("Summary table is missing columns: ",
       paste(missing_sum, collapse = ", "))
}

if (!"gene" %in% names(ppi_nodes)) {
  stop("PPI node table must contain a column named 'gene'.")
}
if (!"hub_type" %in% names(ppi_nodes)) {
  # fallback: اگر قبلاً hub_type ساخته نشده
  ppi_nodes <- ppi_nodes %>%
    dplyr::mutate(
      hub_type = dplyr::if_else(gene %in% c("TOP2A", "PARP1"),
                                "hub_core", "core")
    )
}

## ----------------------------------------------------------------------
## 3. choose EV-miRNAs to plot
## ----------------------------------------------------------------------
mi_sel <- summary_mi %>%
  dplyr::arrange(dplyr::desc(n_targets_PPI_core), adj.P.Val) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::mutate(
    mirna_dir = dplyr::case_when(
      logFC > 0 ~ "miRNA up (DLBCL EVs)",
      logFC < 0 ~ "miRNA down (DLBCL EVs)",
      TRUE      ~ "miRNA (neutral)"
    )
  )

message("--- Selected EV-miRNAs for plotting ---")
print(
  mi_sel %>%
    dplyr::select(
      mature_mirna_id, n_targets_PPI_core,
      n_targets_BROAD, logFC, adj.P.Val
    )
)

## ----------------------------------------------------------------------
## 4. edge list restricted to PPI core genes
## ----------------------------------------------------------------------
edges_core <- targets %>%
  dplyr::filter(
    mature_mirna_id %in% mi_sel$mature_mirna_id,
    in_PPI_core == TRUE
  ) %>%
  dplyr::select(mature_mirna_id, target_symbol) %>%
  dplyr::left_join(
    mi_sel %>%
      dplyr::select(mature_mirna_id, mirna_dir),
    by = "mature_mirna_id"
  )

## ----------------------------------------------------------------------
## 5. node tables
## ----------------------------------------------------------------------
mi_order <- mi_sel %>%
  dplyr::arrange(dplyr::desc(mirna_dir), mature_mirna_id) %>%
  dplyr::pull(mature_mirna_id)

n_mir <- length(mi_order)

nodes_mir <- tibble(
  node_id = mi_order,
  x       = 0,
  y       = seq(from = n_mir, to = 1)
) %>%
  dplyr::left_join(
    mi_sel %>%
      dplyr::select(mature_mirna_id, mirna_dir),
    by = c("node_id" = "mature_mirna_id")
  ) %>%
  dplyr::mutate(
    node_group = mirna_dir
  )

genes_core <- edges_core %>%
  dplyr::pull(target_symbol) %>%
  unique()

nodes_gene <- tibble(node_id = genes_core) %>%
  dplyr::left_join(
    ppi_nodes %>% dplyr::select(gene, hub_type),
    by = c("node_id" = "gene")
  ) %>%
  dplyr::mutate(
    hub_type = dplyr::if_else(is.na(hub_type), "core", hub_type),
    node_group = dplyr::if_else(
      node_id %in% c("TOP2A", "PARP1"),
      "PPI hub (TOP2A/PARP1)",
      "PPI core gene"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(node_group == "PPI hub (TOP2A/PARP1)"),
    node_id
  ) %>%
  dplyr::mutate(
    x = 1,
    y = seq(from = dplyr::n(), to = 1)
  )

## ----------------------------------------------------------------------
## 6. attach coordinates to edges
## ----------------------------------------------------------------------
edges_plot <- edges_core %>%
  dplyr::inner_join(
    nodes_mir %>%
      dplyr::select(mature_mirna_id = node_id, x_mir = x, y_mir = y),
    by = "mature_mirna_id"
  ) %>%
  dplyr::inner_join(
    nodes_gene %>%
      dplyr::select(target_symbol = node_id, x_gene = x, y_gene = y),
    by = "target_symbol"
  )

## ----------------------------------------------------------------------
## 7. combined nodes + factors
## ----------------------------------------------------------------------
nodes_all <- dplyr::bind_rows(
  nodes_mir %>%
    dplyr::transmute(node_id, x, y, node_group),
  nodes_gene %>%
    dplyr::transmute(node_id, x, y, node_group)
) %>%
  dplyr::mutate(
    node_group = factor(
      node_group,
      levels = c(
        "miRNA down (DLBCL EVs)",
        "miRNA up (DLBCL EVs)",
        "PPI core gene",
        "PPI hub (TOP2A/PARP1)"
      )
    )
  )

## ----------------------------------------------------------------------
## 8. plot
## ----------------------------------------------------------------------
message("Building EV-miRNA network plot ...")

edge_col_scale <- scale_colour_manual(
  values = c(
    "miRNA down (DLBCL EVs)" = "#6baed6",
    "miRNA up (DLBCL EVs)"   = "#fb6a4a",
    "miRNA (neutral)"        = "grey70"
  ),
  guide = "none"
)

node_shape_scale <- scale_shape_manual(
  name = NULL,
  values = c(
    "miRNA down (DLBCL EVs)" = 21,
    "miRNA up (DLBCL EVs)"   = 21,
    "PPI core gene"          = 22,
    "PPI hub (TOP2A/PARP1)"  = 22
  )
)

node_fill_scale <- scale_fill_manual(
  name = NULL,
  values = c(
    "miRNA down (DLBCL EVs)" = "#2b8cbe",
    "miRNA up (DLBCL EVs)"   = "#de2d26",
    "PPI core gene"          = "#31a354",
    "PPI hub (TOP2A/PARP1)"  = "#ff7f00"
  )
)

p <- ggplot() +
  geom_curve(
    data = edges_plot,
    aes(
      x = x_mir, y = y_mir,
      xend = x_gene, yend = y_gene,
      colour = mirna_dir
    ),
    curvature = 0.25,
    size      = 0.3,
    alpha     = 0.6
  ) +
  edge_col_scale +
  geom_point(
    data = nodes_all,
    aes(x = x, y = y, fill = node_group, shape = node_group),
    size   = 4,
    colour = "black",
    stroke = 0.4
  ) +
  node_shape_scale +
  node_fill_scale +
  geom_text(
    data = nodes_mir,
    aes(x = x - 0.03, y = y, label = node_id),
    hjust = 1,
    size  = 3.1
  ) +
  geom_text(
    data = nodes_gene,
    aes(x = x + 0.03, y = y, label = node_id),
    hjust = 0,
    size  = 3.1
  ) +
  scale_x_continuous(
    limits = c(-0.4, 1.2),
    breaks = c(0, 1),
    labels = c("EV-miRNAs (GSE171272)", "PPI module genes")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Network of EV-miRNAs targeting cross-species PPI core (TOP2A/PARP1)",
    x     = NULL,
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.x      = element_text(size = 10),
    plot.title       = element_text(
      hjust = 0,
      face  = "bold",
      size  = 13,
      margin = margin(b = 8)
    ),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.text      = element_text(size = 9),
    plot.margin      = margin(t = 10, r = 12, b = 10, l = 210),
    panel.border     = element_rect(colour = "black", fill = NA)
  ) +
  coord_cartesian(clip = "off")

## ----------------------------------------------------------------------
## 9. save
## ----------------------------------------------------------------------
out_png <- file.path(fig_dir, "Fig10D_EVmiRNA_PPIcore_network_GSE171272.png")
out_pdf <- file.path(fig_dir, "Fig10D_EVmiRNA_PPIcore_network_GSE171272.pdf")

ggsave(out_png, p, width = 9.5, height = 5.8, dpi = 400)
ggsave(out_pdf, p, width = 9.5, height = 5.8)

message("Saved EV-miRNA network figure to:")
message("  ", out_png)
message("  ", out_pdf)
message("=== Phase 10D completed successfully. ===")
