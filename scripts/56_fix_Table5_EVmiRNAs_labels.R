## 56_fix_Table5_EVmiRNAs_labels.R
project_root <- "D:/Research/My Articles/DLBCL drug"

tab_path <- file.path(
  project_root,
  "results", "tables",
  "Table5_EVmiRNAs_targeting_PPI_core_selected.tsv"
)

library(readr)
library(dplyr)

tab <- read_tsv(tab_path, show_col_types = FALSE)

tab_fixed <- tab %>%
  mutate(
    regulation_EV = if_else(logFC > 0, "Up in DLBCL EVs", "Down in DLBCL EVs")
  )

out_path <- file.path(
  project_root,
  "results", "tables",
  "Table5_EVmiRNAs_targeting_PPI_core_selected_fixed.tsv"
)

write_tsv(tab_fixed, out_path)

message("Saved fixed Table 5 to:\n  ", out_path)
