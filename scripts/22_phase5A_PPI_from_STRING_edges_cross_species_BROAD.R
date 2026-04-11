## =====================================================================
## Phase 5A: Build PPI network & hubs from STRING edges
##   - Uses cross-species BROAD human gene set
##   - Input:
##       results/tables/network/STRING_input_human_cross_species_BROAD_full.tsv
##       results/tables/network/STRING_network_human_cross_species_BROAD_edges.tsv
##   - Output:
##       results/tables/network/PPI_human_cross_species_BROAD_nodes_with_centrality.tsv
##       results/tables/network/PPI_human_cross_species_BROAD_top_hubs.tsv
## =====================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(igraph)
  library(stringr)
})

project_root <- "D:/Research/My Articles/DLBCL drug"

path_proj <- function(...) file.path(project_root, ...)

net_dir <- path_proj("results", "tables", "network")
dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 0. Paths
## ---------------------------------------------------------------------

nodes_path_full <- file.path(
  net_dir,
  "STRING_input_human_cross_species_BROAD_full.tsv"
)

edges_path <- file.path(
  net_dir,
  "STRING_network_human_cross_species_BROAD_edges.tsv"
)

if (!file.exists(nodes_path_full)) {
  stop("[PPI] Node file not found: ", nodes_path_full)
}

if (!file.exists(edges_path)) {
  stop("[PPI] STRING edges file not found: ", edges_path,
       "\nExport it from STRING and save as this path.")
}

## ---------------------------------------------------------------------
## 1. Load node table (genes + DE info)
## ---------------------------------------------------------------------

nodes <- readr::read_tsv(nodes_path_full, show_col_types = FALSE)

if (!"human_symbol" %in% colnames(nodes)) {
  stop("[PPI] 'human_symbol' column not found in node file.")
}

nodes <- nodes %>%
  mutate(
    gene       = human_symbol,
    gene_upper = toupper(human_symbol)
  )

message("[PPI] Loaded node table: ", nrow(nodes), " genes.")

## ---------------------------------------------------------------------
## 2. Load STRING edges
## ---------------------------------------------------------------------

edges_raw <- readr::read_tsv(edges_path, show_col_types = FALSE)

message("[PPI] Loaded edges: ", nrow(edges_raw), " rows.")
message("[PPI] Edge columns (first 15): ",
        paste(head(colnames(edges_raw), 15), collapse = ", "))

cols <- colnames(edges_raw)

## تشخیص ستون‌های نام گره‌ها
nameA_col <- dplyr::case_when(
  "preferredName_A" %in% cols ~ "preferredName_A",
  "protein1"        %in% cols ~ "protein1",
  "#node1"          %in% cols ~ "#node1",
  "node1"           %in% cols ~ "node1",
  TRUE ~ NA_character_
)

nameB_col <- dplyr::case_when(
  "preferredName_B" %in% cols ~ "preferredName_B",
  "protein2"        %in% cols ~ "protein2",
  "node2"           %in% cols ~ "node2",
  TRUE ~ NA_character_
)

if (is.na(nameA_col) || is.na(nameB_col)) {
  stop("[PPI] Could not find any suitable node columns in edges.\n",
       "Available columns: ", paste(cols, collapse = ", "))
}

## ستون score اگر وجود داشته باشد
score_col <- if ("combined_score" %in% cols) {
  "combined_score"
} else if ("experimental" %in% cols) {
  "experimental"
} else if ("experimentally_determined_interaction" %in% cols) {
  "experimentally_determined_interaction"
} else {
  NA_character_
}

edges <- edges_raw %>%
  transmute(
    geneA = .data[[nameA_col]],
    geneB = .data[[nameB_col]],
    score = if (!is.na(score_col)) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_
  ) %>%
  filter(!is.na(geneA), !is.na(geneB)) %>%
  mutate(
    geneA_upper = toupper(geneA),
    geneB_upper = toupper(geneB)
  )

message("[PPI] Clean edges: ", nrow(edges), " rows.")

if (!all(is.na(edges$score))) {
  message("[PPI] Score column used: ", score_col)
  score_range <- range(edges$score, na.rm = TRUE)
  message("[PPI] Score range: [", score_range[1], ", ", score_range[2], "]")
} else {
  message("[PPI] No numeric score column found; proceeding without score threshold.")
}

## ---------------------------------------------------------------------
## 3. بدون threshold روی score — فقط به module genes محدود می‌کنیم
## ---------------------------------------------------------------------

genes_upper <- unique(nodes$gene_upper)

edges_mod <- edges %>%
  filter(
    geneA_upper %in% genes_upper &
      geneB_upper %in% genes_upper
  )

message("[PPI] Edges within module genes (no score threshold): ",
        nrow(edges_mod), " rows.")

if (nrow(edges_mod) == 0L) {
  stop("[PPI] No edges remain after restricting to module genes.")
}

## ---------------------------------------------------------------------
## 4. Build igraph network
## ---------------------------------------------------------------------

g <- igraph::graph_from_data_frame(
  d = edges_mod %>%
    select(geneA_upper, geneB_upper) %>%
    rename(from = geneA_upper, to = geneB_upper),
  directed = FALSE
)

message("[PPI] Graph: ", igraph::gorder(g), " nodes, ",
        igraph::gsize(g), " edges.")

## ---------------------------------------------------------------------
## 5. Centrality measures
## ---------------------------------------------------------------------

deg  <- igraph::degree(g, mode = "all")
bet  <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
clos <- igraph::closeness(g, normalized = TRUE)

centrality_tbl <- tibble::tibble(
  gene_upper  = names(deg),
  degree      = as.numeric(deg),
  betweenness = as.numeric(bet),
  closeness   = as.numeric(clos)
)

## community detection (Louvain)
cl_louvain <- igraph::cluster_louvain(g)
community_tbl <- tibble::tibble(
  gene_upper = names(igraph::membership(cl_louvain)),
  community  = as.integer(igraph::membership(cl_louvain))
)

## merge به node table
nodes_cent <- nodes %>%
  left_join(centrality_tbl, by = "gene_upper") %>%
  left_join(community_tbl,  by = "gene_upper")

## ---------------------------------------------------------------------
## 6. Define hubs (۸۰امین پرسن‌تایل degree)
## ---------------------------------------------------------------------

deg_non_na <- nodes_cent$degree[!is.na(nodes_cent$degree)]
if (length(deg_non_na) == 0L) {
  warning("[PPI] All degree values are NA.")
  deg_cutoff <- NA_real_
} else {
  deg_cutoff <- stats::quantile(deg_non_na, probs = 0.80, na.rm = TRUE)
}

nodes_cent <- nodes_cent %>%
  mutate(
    is_hub = ifelse(!is.na(degree) & degree >= deg_cutoff, TRUE, FALSE)
  )

message("[PPI] Degree 80th percentile (hub cutoff): ", deg_cutoff)

## ---------------------------------------------------------------------
## 7. Write outputs
## ---------------------------------------------------------------------

out_nodes_path <- file.path(
  net_dir,
  "PPI_human_cross_species_BROAD_nodes_with_centrality.tsv"
)

readr::write_tsv(nodes_cent, out_nodes_path)
message("[PPI] Node+centrality table written to: ", out_nodes_path)

top_n_hubs <- 50L
hubs_tbl <- nodes_cent %>%
  filter(!is.na(degree)) %>%
  arrange(desc(degree), desc(betweenness)) %>%
  slice(1:top_n_hubs)

out_hubs_path <- file.path(
  net_dir,
  "PPI_human_cross_species_BROAD_top_hubs.tsv"
)

readr::write_tsv(hubs_tbl, out_hubs_path)
message("[PPI] Top ", top_n_hubs, " hubs written to: ", out_hubs_path)

message("=== [Phase 5A - PPI network & hubs for cross-species BROAD] DONE ===")
