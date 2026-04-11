## 43_phase7_fig_PPI_network_with_CMap_hits.R
## Network-anchored CMap hits within the cross-species BROAD PPI module

message("=== Phase 7 (Figures): Network-anchored CMap hits in cross-species BROAD PPI ===")

required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "igraph", "ggraph", "ggplot2"
)

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
  library(igraph)
  library(ggraph)
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

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

path_edges      <- file.path(net_dir, "PPI_cross_species_BROAD_edges_filtered.tsv")
path_nodes      <- file.path(net_dir, "PPI_cross_species_BROAD_nodes.tsv")
path_hubs_drugs <- file.path(net_dir, "PPI_cross_species_BROAD_hubs_with_drug_hits.tsv")

for (p in c(path_edges, path_nodes, path_hubs_drugs)) {
  if (!file.exists(p)) {
    stop("Required network file not found at:\n  ", p)
  }
}

## --------------------------------------------------------------------
## 2. read edges / nodes / hubs+drugs
## --------------------------------------------------------------------
message("--- Step 1: Read PPI edges & node attributes ---")

edges <- readr::read_tsv(path_edges, show_col_types = FALSE)
nodes <- readr::read_tsv(path_nodes, show_col_types = FALSE)
hubs_hits <- readr::read_tsv(path_hubs_drugs, show_col_types = FALSE)

message("  Raw edges rows:  ", nrow(edges))
message("  Raw nodes rows:  ", nrow(nodes))
message("  Hubs-with-drugs: ", length(unique(hubs_hits$gene)))
message("  Core drugs:      ",
        paste(sort(unique(unlist(strsplit(hubs_hits$drugs, ";")))),
              collapse = ", "))

## detect edge gene columns
from_candidates <- intersect(
  c("geneA", "preferredName_A", "node1", "protein1", "from", "source"),
  names(edges)
)
to_candidates <- intersect(
  c("geneB", "preferredName_B", "node2", "protein2", "to", "target"),
  names(edges)
)

if (length(from_candidates) == 0 || length(to_candidates) == 0) {
  stop(
    "Could not detect edge gene columns in PPI edges file.\n",
    "Columns موجود:\n  ", paste(names(edges), collapse = ", ")
  )
}

col_from <- from_candidates[1]
col_to   <- to_candidates[1]

edges2 <- edges %>%
  dplyr::transmute(
    geneA = .data[[col_from]],
    geneB = .data[[col_to]]
  ) %>%
  dplyr::filter(!is.na(geneA), !is.na(geneB))

## ensure node table has 'gene'
if (!"gene" %in% names(nodes)) {
  if ("name" %in% names(nodes)) {
    nodes <- nodes %>% dplyr::rename(gene = name)
  } else {
    stop(
      "Node table must contain a 'gene' (or 'name') column.\n",
      "Columns موجود در nodes:\n  ", paste(names(nodes), collapse = ", ")
    )
  }
}

## if hub_type missing, create it to avoid coalesce error
if (!"hub_type" %in% names(nodes)) {
  nodes <- nodes %>% dplyr::mutate(hub_type = NA_character_)
}

hubs_hits2 <- hubs_hits %>%
  dplyr::select(gene, hub_type, n_drugs, drugs) %>%
  dplyr::mutate(
    hub_type = stringr::str_replace_all(hub_type, "-", "_")
  )

## --------------------------------------------------------------------
## 3. annotate nodes
## --------------------------------------------------------------------
message("--- Step 2: Annotate nodes with hub/bottleneck and CMap hits ---")

nodes_annot <- nodes %>%
  dplyr::left_join(
    hubs_hits2,
    by = "gene",
    suffix = c("", "_from_hubs")
  ) %>%
  dplyr::mutate(
    hub_type = dplyr::coalesce(.data$hub_type, .data$hub_type_from_hubs, "other"),
    hub_class = dplyr::case_when(
      hub_type %in% c("hub_bottleneck", "hub_bottleneck_core", "hub_bottleneck_module") ~ "Hub–bottleneck",
      hub_type == "hub" ~ "Hub",
      hub_type == "bottleneck" ~ "Bottleneck",
      TRUE ~ "Other"
    ),
    hub_class = factor(
      hub_class,
      levels = c("Hub–bottleneck", "Hub", "Bottleneck", "Other")
    ),
    has_drug = !is.na(n_drugs) & n_drugs > 0
  ) %>%
  dplyr::select(-hub_type_from_hubs)

## make sure all genes in edges are present in vertices
genes_in_edges <- sort(unique(c(edges2$geneA, edges2$geneB)))
missing_genes <- setdiff(genes_in_edges, nodes_annot$gene)

if (length(missing_genes) > 0) {
  message("  [Warning] ", length(missing_genes),
          " genes from edges not found in node table; اضافه می‌کنم به صورت 'Other'.")
  nodes_annot <- dplyr::bind_rows(
    nodes_annot,
    tibble(
      gene       = missing_genes,
      hub_type   = "other",
      hub_class  = factor("Other",
                          levels = c("Hub–bottleneck", "Hub", "Bottleneck", "Other")),
      has_drug   = FALSE,
      n_drugs    = NA_integer_,
      drugs      = NA_character_
    )
  )
}

## --------------------------------------------------------------------
## 4. build igraph
## --------------------------------------------------------------------
message("--- Step 3: Build igraph object ---")

g <- igraph::graph_from_data_frame(
  d = edges2 %>% dplyr::rename(from = geneA, to = geneB),
  directed = FALSE,
  vertices = nodes_annot
)

message("  Graph: ", igraph::vcount(g), " nodes, ",
        igraph::ecount(g), " edges")

## NOTE:
## igraph/ggraph expose vertex name as column `name` in node data.
## ما gene را به عنوان vertex attribute داریم، ولی برای label و فیلتر
## امن‌تر است از `name` استفاده کنیم که با edges سازگار است.

## --------------------------------------------------------------------
## 5. plot with ggraph
## --------------------------------------------------------------------
message("--- Step 4: Build network figure ---")

set.seed(1234)

p <- ggraph(g, layout = "fr") +
  geom_edge_link(
    edge_colour = "grey85",
    edge_alpha  = 0.25,
    edge_width  = 0.3
  ) +
  # پایه نودها: کلاس شبکه
  geom_node_point(
    aes(fill = hub_class),
    shape  = 21,
    size   = 3,
    colour = "grey20",
    stroke = 0.25,
    alpha  = 0.95
  ) +
  # نودهایی که دارو دارند: حلقه‌ی مشکی
  geom_node_point(
    data = function(df) dplyr::filter(df, has_drug),
    shape  = 21,
    size   = 4,
    fill   = NA,
    colour = "black",
    stroke = 0.9
  ) +
  # highlight ویژه برای TOP2A / PARP1 (بر اساس name)
  geom_node_point(
    data = function(df) dplyr::filter(df, name %in% c("TOP2A", "PARP1")),
    shape  = 21,
    size   = 4.5,
    fill   = "white",
    colour = "black",
    stroke = 1.2
  ) +
  geom_node_text(
    data = function(df) dplyr::filter(df, name %in% c("TOP2A", "PARP1")),
    aes(label = name),
    size     = 3,
    fontface = "bold",
    vjust    = -1,
    family   = ""
  ) +
  scale_fill_manual(
    name   = "Node class",
    values = c(
      "Hub–bottleneck" = "#b2182b",
      "Hub"            = "#ef8a62",
      "Bottleneck"     = "#67a9cf",
      "Other"          = "#d9d9d9"
    )
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    colour = "none"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 9),
    plot.title      = element_text(
      hjust  = 0.5,
      face   = "bold",
      size   = 14,
      margin = margin(b = 8)
    ),
    plot.margin     = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  ggtitle("Network-anchored CMap hits in the cross-species DLBCL module")

out_png <- file.path(fig_dir, "Fig5A_PPI_network_core_CMap_hits.png")
out_pdf <- file.path(fig_dir, "Fig5A_PPI_network_core_CMap_hits.pdf")

ggsave(out_png, p, width = 7.5, height = 6.5, dpi = 450)
ggsave(out_pdf, p, width = 7.5, height = 6.5)

message("Saved network figure to:")
message("  ", out_png)
message("  ", out_pdf)
message("=== Phase 7 (network-anchored CMap hits) completed successfully. ===")
