## 31_phase8_figures_crossSpecies_logFC_scatter.R
## هدف:
##  - Fig. 3C: cross-species log2FC scatter (human vs canine) برای ماژول BROAD
##  - نمایش:
##      * STRICT core vs BROAD-only (Module membership)
##      * Hub / Bottleneck / Hub–bottleneck / Other (Network class)
##      * برچسب تمیز TOP2A و PARP1 با ggrepel

message("=== Phase 8 (Figures): Cross-species logFC scatter (Fig. 3C) ===")

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
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------------------------------------------------------------------
## 1. paths
## --------------------------------------------------------------------
## به‌صورت صریح روت پروژه DLBCL رو ست می‌کنیم
project_root <- normalizePath("D:/Research/My Articles/DLBCL drug",
                              winslash = "/",
                              mustWork = TRUE)
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

if (!file.exists(path_broad)) {
  stop("BROAD cross-species table not found at:\n  ", path_broad)
}
if (!file.exists(path_strict)) {
  stop("STRICT core cross-species table not found at:\n  ", path_strict)
}
if (!file.exists(path_nodes)) {
  stop("PPI nodes table not found at:\n  ", path_nodes)
}

## --------------------------------------------------------------------
## 2. helpers
## --------------------------------------------------------------------
get_col <- function(df, candidates, label) {
  present <- candidates[candidates %in% names(df)]
  if (length(present) == 0) {
    stop(
      "Could not find any of ",
      paste(candidates, collapse = ", "),
      " in ", label, "."
    )
  }
  present[1]
}

## --------------------------------------------------------------------
## 3. load BROAD & STRICT tables
## --------------------------------------------------------------------
message("Reading cross-species BROAD & STRICT tables ...")

broad  <- readr::read_tsv(path_broad,  show_col_types = FALSE)
strict <- readr::read_tsv(path_strict, show_col_types = FALSE)

human_col  <- get_col(
  broad,
  c("human_symbol", "human_gene", "Human_symbol", "human"),
  "BROAD cross-species table"
)
dog_col    <- get_col(
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

hFDR_col <- intersect(
  c("human_FDR", "human_adjP", "human_FDR_tumor_vs_normal"),
  names(broad)
)
dFDR_col <- intersect(
  c("dog_FDR", "dog_adjP", "dog_FDR_tumor_vs_normal"),
  names(broad)
)

broad2 <- broad %>%
  mutate(
    human_symbol = .data[[human_col]],
    dog_symbol   = .data[[dog_col]],
    human_logFC  = as.numeric(.data[[hlogFC_col]]),
    dog_logFC    = as.numeric(.data[[dlogFC_col]])
  ) %>%
  select(
    human_symbol, dog_symbol,
    human_logFC, dog_logFC,
    everything()
  )

strict2 <- strict %>%
  mutate(
    human_symbol = .data[[get_col(
      strict,
      c("human_symbol", "human_gene", "Human_symbol", "human"),
      "STRICT cross-species table"
    )]],
    dog_symbol   = .data[[get_col(
      strict,
      c("dog_symbol", "dog_gene", "Dog_symbol", "dog"),
      "STRICT cross-species table"
    )]]
  ) %>%
  select(human_symbol, dog_symbol, everything())

strict_pairs <- strict2 %>%
  mutate(pair_id = paste0(human_symbol, "|", dog_symbol)) %>%
  pull(pair_id) %>%
  unique()

plot_df <- broad2 %>%
  mutate(
    pair_id    = paste0(human_symbol, "|", dog_symbol),
    membership = if_else(pair_id %in% strict_pairs,
                         "STRICT core",
                         "BROAD-only")
  )

## --------------------------------------------------------------------
## 4. join PPI nodes (hub_type)
## --------------------------------------------------------------------
nodes <- readr::read_tsv(path_nodes, show_col_types = FALSE)

if (!"gene" %in% names(nodes)) {
  stop("PPI nodes table lacks 'gene' column.")
}
if (!"hub_type" %in% names(nodes)) {
  nodes$hub_type <- "non_hub"
}

nodes2 <- nodes %>%
  mutate(
    hub_type = case_when(
      hub_type %in% c("hub_bottleneck", "hub-bottleneck") ~ "hub_bottleneck",
      hub_type %in% c("bottleneck")                       ~ "bottleneck",
      hub_type %in% c("hub")                              ~ "hub",
      TRUE                                                ~ "non_hub"
    ),
    hub_type = factor(
      hub_type,
      levels = c("non_hub", "hub", "bottleneck", "hub_bottleneck")
    )
  ) %>%
  select(gene, hub_type)

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

# اگر FDRها باشند (اختیاری، فعلاً در plot استفاده نمی‌کنیم)
if (length(hFDR_col) == 1 && length(dFDR_col) == 1) {
  plot_df <- plot_df %>%
    mutate(
      human_FDR = as.numeric(.data[[hFDR_col]]),
      dog_FDR   = as.numeric(.data[[dFDR_col]]),
      mean_negLog10FDR = -log10(pmax(human_FDR, dog_FDR, na.rm = TRUE))
    )
} else {
  plot_df$mean_negLog10FDR <- NA_real_
}

# فیلتر جفت‌هایی که logFC ندارند
plot_df <- plot_df %>%
  filter(!is.na(human_logFC), !is.na(dog_logFC))

message("  - N ortholog pairs in BROAD module (with logFC): ", nrow(plot_df))
message("  - STRICT core pairs: ",
        sum(plot_df$membership == "STRICT core", na.rm = TRUE))

## --------------------------------------------------------------------
## 5. aesthetics
## --------------------------------------------------------------------
range_all <- range(c(plot_df$human_logFC, plot_df$dog_logFC), na.rm = TRUE)
pad       <- diff(range_all) * 0.08
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
## 6. labels برای TOP2A و PARP1
## --------------------------------------------------------------------
label_genes <- c("TOP2A", "PARP1")
df_labels   <- plot_df %>%
  filter(human_symbol %in% label_genes)

if (nrow(df_labels) == 0) {
  warning("No TOP2A / PARP1 rows found in plot_df – labels will be skipped.")
}

## --------------------------------------------------------------------
## 7. build scatter plot (Fig. 3C)
## --------------------------------------------------------------------
p_scatter <- ggplot(
  plot_df,
  aes(x = human_logFC, y = dog_logFC)
) +
  # همه ژن‌ها
  geom_point(
    aes(
      shape = hub_type,
      fill  = membership
    ),
    colour = "grey30",
    size   = 2.3,
    alpha  = 0.8,
    stroke = 0.25
  ) +
  # خط y = x
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "dashed",
    colour    = "grey70",
    size      = 0.4
  ) +
  # لیبل تمیز برای TOP2A و PARP1
  {
    if (nrow(df_labels) > 0) {
      geom_label_repel(
        data  = df_labels,
        aes(label = human_symbol),
        size          = 3,
        fontface      = "bold",
        label.size    = 0.2,
        label.padding = unit(0.15, "lines"),
        nudge_x       = c(0.4, -0.4),
        nudge_y       = c(-0.3, 0.3),
        min.segment.length = 0,
        segment.colour     = "grey40",
        segment.size       = 0.3,
        max.overlaps       = Inf
      )
    }
  } +
  scale_shape_manual(
    values = shape_hub,
    name   = "Network class"
  ) +
  scale_fill_manual(
    values = col_membership,
    name   = "Module membership"
  ) +
  coord_equal(xlim = limits_xy, ylim = limits_xy, expand = FALSE) +
  labs(
    x = "Human log2 fold-change (DLBCL vs normal)",
    y = "Canine log2 fold-change (DLBCL vs normal)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", size = 0.3),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.5, "lines"),
    legend.box       = "vertical",
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5)
  )

## --------------------------------------------------------------------
## 8. save
## --------------------------------------------------------------------
fig_file_png <- file.path(fig_dir, "Fig3C_crossSpecies_logFC_scatter.png")
fig_file_pdf <- file.path(fig_dir, "Fig3C_crossSpecies_logFC_scatter.pdf")

ggsave(fig_file_png, p_scatter, width = 6.5, height = 5.5, units = "in", dpi = 300)
ggsave(fig_file_pdf, p_scatter, width = 6.5, height = 5.5, units = "in")

message("Saved Fig. 3C scatter to:")
message("  - ", fig_file_png)
message("  - ", fig_file_pdf)
message("=== Phase 8 (Figures): Cross-species logFC scatter completed successfully. ===")