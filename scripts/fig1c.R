
## fig1c.R
## Cross-species log2FC scatter for the BROAD module
## Final revised version for manuscript / reviewer response

message("=== Phase 8 (Figures): Cross-species logFC scatter ===")

required_pkgs <- c(
  "readr", "dplyr", "tibble", "stringr",
  "ggplot2", "ggrepel", "tidyr"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Please run: install.packages(c(",
      paste(sprintf('\"%s\"', required_pkgs), collapse = ", "),
      "))"
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
project_root <- normalizePath(
  "D:/Research/My Articles/DLBCL drug",
  winslash = "/",
  mustWork = TRUE
)
message("Project root: ", project_root)

sig_dir <- file.path(project_root, "results", "tables", "signatures")
net_dir <- file.path(project_root, "results", "tables", "network")
fig_dir <- file.path(project_root, "results", "figures")

for (d in c(sig_dir, net_dir)) {
  if (!dir.exists(d)) {
    stop("Expected directory not found:\n  ", d)
  }
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
  message("Created figures directory:\n  ", fig_dir)
}

path_broad  <- file.path(sig_dir, "cross_species_module_BROAD_table.tsv")
path_strict <- file.path(sig_dir, "cross_species_module_STRICT_core_table.tsv")
path_nodes  <- file.path(net_dir,  "PPI_cross_species_BROAD_nodes.tsv")

for (p in c(path_broad, path_strict, path_nodes)) {
  if (!file.exists(p)) {
    stop("Required file not found:\n  ", p)
  }
}

## --------------------------------------------------------------------
## 2. helpers
## --------------------------------------------------------------------
get_col <- function(df, candidates, label, required = TRUE) {
  present <- candidates[candidates %in% names(df)]
  if (length(present) == 0) {
    if (required) {
      stop(
        "Could not find any of ",
        paste(candidates, collapse = ", "),
        " in ", label, "."
      )
    } else {
      return(NA_character_)
    }
  }
  present[1]
}

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
## 3. load BROAD & STRICT tables
## --------------------------------------------------------------------
message("Reading cross-species BROAD & STRICT tables ...")

broad  <- readr::read_tsv(path_broad,  show_col_types = FALSE)
strict <- readr::read_tsv(path_strict, show_col_types = FALSE)

human_col <- get_col(
  broad,
  c("human_symbol", "human_gene", "Human_symbol", "human"),
  "BROAD cross-species table"
)
dog_col <- get_col(
  broad,
  c("dog_symbol", "dog_gene", "Dog_symbol", "dog"),
  "BROAD cross-species table"
)
hlogFC_col <- get_col(
  broad,
  c("human_logFC", "logFC_human", "human_logFC_tumor_vs_normal"),
  "BROAD cross-species table"
)
dlogFC_col <- get_col(
  broad,
  c("dog_logFC", "logFC_dog", "dog_logFC_tumor_vs_normal"),
  "BROAD cross-species table"
)

hFDR_col <- get_col(
  broad,
  c("human_FDR", "human_adjP", "human_FDR_tumor_vs_normal", "human_adj.P.Val"),
  "BROAD cross-species table",
  required = FALSE
)
dFDR_col <- get_col(
  broad,
  c("dog_FDR", "dog_adjP", "dog_FDR_tumor_vs_normal", "dog_adj.P.Val"),
  "BROAD cross-species table",
  required = FALSE
)

broad2 <- broad %>%
  mutate(
    human_symbol = as.character(.data[[human_col]]),
    dog_symbol   = as.character(.data[[dog_col]]),
    human_logFC  = as.numeric(.data[[hlogFC_col]]),
    dog_logFC    = as.numeric(.data[[dlogFC_col]])
  ) %>%
  select(human_symbol, dog_symbol, human_logFC, dog_logFC, everything())

strict2 <- strict %>%
  mutate(
    human_symbol = as.character(.data[[get_col(
      strict,
      c("human_symbol", "human_gene", "Human_symbol", "human"),
      "STRICT cross-species table"
    )]]),
    dog_symbol = as.character(.data[[get_col(
      strict,
      c("dog_symbol", "dog_gene", "Dog_symbol", "dog"),
      "STRICT cross-species table"
    )]])
  ) %>%
  select(human_symbol, dog_symbol, everything())

strict_pairs <- strict2 %>%
  mutate(pair_id = paste0(human_symbol, "|", dog_symbol)) %>%
  pull(pair_id) %>%
  unique()

plot_df <- broad2 %>%
  mutate(
    pair_id = paste0(human_symbol, "|", dog_symbol),
    membership = if_else(pair_id %in% strict_pairs, "STRICT core", "BROAD-only")
  )

## --------------------------------------------------------------------
## 4. join PPI nodes
## --------------------------------------------------------------------
message("Reading PPI node table ...")
nodes <- readr::read_tsv(path_nodes, show_col_types = FALSE)

gene_col_nodes <- get_col(
  nodes,
  c("gene", "human_symbol", "symbol"),
  "PPI node table"
)

hub_col <- get_col(
  nodes,
  c("hub_type", "node_class", "class"),
  "PPI node table",
  required = FALSE
)

deg_col <- get_col(
  nodes,
  c("degree", "degree_centrality", "deg"),
  "PPI node table",
  required = FALSE
)

btw_col <- get_col(
  nodes,
  c("betweenness", "betweenness_centrality", "btw"),
  "PPI node table",
  required = FALSE
)

nodes2 <- nodes %>%
  mutate(
    gene = as.character(.data[[gene_col_nodes]]),
    hub_type = if (!is.na(hub_col)) normalize_hub_type(.data[[hub_col]]) else "non_hub",
    degree_use = if (!is.na(deg_col)) as.numeric(.data[[deg_col]]) else NA_real_,
    betweenness_use = if (!is.na(btw_col)) as.numeric(.data[[btw_col]]) else NA_real_
  ) %>%
  select(gene, hub_type, degree_use, betweenness_use)

plot_df <- plot_df %>%
  left_join(nodes2, by = c("human_symbol" = "gene")) %>%
  mutate(
    hub_type = tidyr::replace_na(hub_type, "non_hub"),
    hub_type = factor(
      hub_type,
      levels = c("non_hub", "hub", "bottleneck", "hub_bottleneck")
    ),
    membership = factor(
      membership,
      levels = c("STRICT core", "BROAD-only")
    )
  )

if (!is.na(hFDR_col) && !is.na(dFDR_col)) {
  plot_df <- plot_df %>%
    mutate(
      human_FDR = as.numeric(.data[[hFDR_col]]),
      dog_FDR   = as.numeric(.data[[dFDR_col]]),
      mean_negLog10FDR = -log10(pmax(human_FDR, dog_FDR, na.rm = TRUE))
    )
} else {
  plot_df$mean_negLog10FDR <- NA_real_
}

plot_df <- plot_df %>%
  filter(!is.na(human_logFC), !is.na(dog_logFC)) %>%
  mutate(
    max_abs_logFC = pmax(abs(human_logFC), abs(dog_logFC), na.rm = TRUE),
    hub_rank = case_when(
      hub_type == "hub_bottleneck" ~ 1L,
      hub_type == "hub"            ~ 2L,
      hub_type == "bottleneck"     ~ 3L,
      TRUE                         ~ 4L
    )
  )

message("  - N ortholog pairs in BROAD module (with logFC): ", nrow(plot_df))
message("  - STRICT core pairs: ", sum(plot_df$membership == "STRICT core", na.rm = TRUE))

## --------------------------------------------------------------------
## 5. choose labelled genes
## --------------------------------------------------------------------
## Fixed core story genes:
priority_genes <- c("TOP2A", "PARP1")
preferred_extra_genes <- c("AURKA", "BUB1B", "BRCA1")

plot_df <- plot_df %>%
  mutate(
    degree_rank = if (all(is.na(degree_use))) {
      NA_integer_
    } else {
      dplyr::dense_rank(dplyr::desc(tidyr::replace_na(degree_use, -Inf)))
    },
    betweenness_rank = if (all(is.na(betweenness_use))) {
      NA_integer_
    } else {
      dplyr::dense_rank(dplyr::desc(tidyr::replace_na(betweenness_use, -Inf)))
    }
  ) %>%
  mutate(
    centrality_rank = case_when(
      !is.na(degree_rank) & !is.na(betweenness_rank) ~ degree_rank + betweenness_rank,
      !is.na(degree_rank) &  is.na(betweenness_rank) ~ degree_rank,
      is.na(degree_rank)  & !is.na(betweenness_rank) ~ betweenness_rank,
      TRUE ~ 9999L
    )
  )

priority_df <- plot_df %>%
  filter(human_symbol %in% priority_genes) %>%
  distinct(human_symbol, .keep_all = TRUE)

preferred_extra_df <- plot_df %>%
  filter(human_symbol %in% preferred_extra_genes) %>%
  distinct(human_symbol, .keep_all = TRUE)

n_missing_extra <- 3 - nrow(preferred_extra_df)

fallback_extra_df <- plot_df %>%
  filter(hub_type %in% c("hub", "bottleneck", "hub_bottleneck")) %>%
  filter(!human_symbol %in% c(priority_genes, preferred_extra_genes)) %>%
  arrange(hub_rank, centrality_rank, desc(max_abs_logFC), human_symbol) %>%
  distinct(human_symbol, .keep_all = TRUE) %>%
  slice_head(n = max(0, n_missing_extra))

label_df <- bind_rows(
  priority_df,
  preferred_extra_df,
  fallback_extra_df
) %>%
  distinct(human_symbol, .keep_all = TRUE)

## desired order in output
label_df <- label_df %>%
  mutate(
    label_order = match(human_symbol, c(priority_genes, preferred_extra_genes))
  ) %>%
  arrange(label_order, centrality_rank, desc(max_abs_logFC), human_symbol) %>%
  select(-label_order)

message("Genes labelled in Fig. 1C:")
cat("  ", paste(label_df$human_symbol, collapse = ", "), "\n")

label_out <- file.path(fig_dir, "Fig1C_labelled_genes.tsv")
label_df %>%
  select(
    human_symbol, dog_symbol, human_logFC, dog_logFC, membership,
    hub_type, degree_use, betweenness_use, centrality_rank, max_abs_logFC
  ) %>%
  readr::write_tsv(label_out)
message("Saved labelled-gene table to:")
message("  ", label_out)

## --------------------------------------------------------------------
## 6. aesthetics
## --------------------------------------------------------------------
range_all <- range(c(plot_df$human_logFC, plot_df$dog_logFC), na.rm = TRUE)
pad <- diff(range_all) * 0.08
limits_xy <- c(range_all[1] - pad, range_all[2] + pad)

col_membership <- c(
  "STRICT core" = "#08519c",
  "BROAD-only"  = "#9ecae1"
)

shape_hub <- c(
  "non_hub"        = 21,  # circle
  "hub"            = 22,  # square
  "bottleneck"     = 24,  # triangle
  "hub_bottleneck" = 23   # diamond
)

## --------------------------------------------------------------------
## 7. build scatter plot
## --------------------------------------------------------------------
p_scatter <- ggplot(plot_df, aes(x = human_logFC, y = dog_logFC)) +
  geom_point(
    aes(shape = hub_type, fill = membership),
    colour = "grey30",
    size   = 2.15,
    alpha  = 0.72,
    stroke = 0.25
  ) +
  geom_point(
    data = label_df,
    aes(
      x = human_logFC,
      y = dog_logFC,
      shape = hub_type,
      fill = membership
    ),
    colour = "black",
    size   = 3.2,
    alpha  = 1,
    stroke = 0.55,
    inherit.aes = FALSE
  ) +
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "dashed",
    colour    = "grey70",
    linewidth = 0.4
  ) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(
      x = human_logFC,
      y = dog_logFC,
      label = human_symbol
    ),
    inherit.aes        = FALSE,
    size               = 3,
    fontface           = "bold",
    label.size         = 0.2,
    label.padding      = grid::unit(0.15, "lines"),
    min.segment.length = 0,
    segment.colour     = "grey40",
    segment.size       = 0.3,
    max.overlaps       = Inf,
    box.padding        = 0.25,
    point.padding      = 0.2,
    seed               = 12345
  ) +
  scale_shape_manual(
    values = shape_hub,
    name   = "Network class",
    labels = c(
      non_hub = "Other",
      hub = "Hub",
      bottleneck = "Bottleneck",
      hub_bottleneck = "Hub–bottleneck"
    )
  ) +
  scale_fill_manual(
    values = col_membership,
    name   = "Module membership"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape  = 21,
        colour = "grey30",
        size   = 4,
        alpha  = 1,
        stroke = 0.35
      ),
      order = 1
    ),
    shape = guide_legend(
      override.aes = list(
        fill   = "white",
        colour = "grey30",
        size   = 4,
        alpha  = 1,
        stroke = 0.6
      ),
      order = 2
    )
  ) +
  coord_equal(xlim = limits_xy, ylim = limits_xy, expand = FALSE) +
  labs(
    x = "Human log2 fold-change (DLBCL vs normal)",
    y = "Canine log2 fold-change (DLBCL vs normal)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    legend.key.size  = grid::unit(0.5, "lines"),
    legend.box       = "vertical"
  )

## --------------------------------------------------------------------
## 8. save
## --------------------------------------------------------------------
fig_file_png <- file.path(fig_dir, "Fig1C_crossSpecies_logFC_scatter.png")
fig_file_pdf <- file.path(fig_dir, "Fig1C_crossSpecies_logFC_scatter.pdf")

ggsave(fig_file_png, p_scatter, width = 6.5, height = 5.5, units = "in", dpi = 300)
ggsave(fig_file_pdf, p_scatter, width = 6.5, height = 5.5, units = "in")

message("Saved cross-species logFC scatter to:")
message("  - ", fig_file_png)
message("  - ", fig_file_pdf)
message("=== Phase 8 (Figures): Cross-species logFC scatter completed successfully. ===")
