# canine-dlbcl-cross-species-replication-module

Curated analysis scripts and supporting metadata for the final manuscript:

**Cross-species transcriptomic and network analysis of human and canine diffuse large B-cell lymphoma identifies a conserved replication–DNA damage module for drug repurposing**

## Overview

This repository contains the curated R scripts and supporting metadata used for the final analyses reported in the manuscript. The study integrates public human and canine DLBCL transcriptomic datasets, pathway enrichment, STRING protein–protein interaction (PPI) networks, Connectivity Map/LINCS drug prioritisation, survival modelling, and extracellular-vesicle microRNA (EV-miRNA) analyses.

The repository is intentionally curated rather than exhaustive: only scripts and small derived files needed to understand and reproduce the final manuscript workflow are included.

## Repository structure

- `scripts/`  
  Curated R scripts used for preprocessing, module derivation, pathway analysis, PPI/CMap analyses, survival modelling, EV-miRNA analysis, and final figure/table generation.

- `metadata/`  
  Dataset metadata and key intermediate tables required for cross-species module definition and survival analyses.

- `derived_tables/`  
  Small derived result tables supporting PPI, drug-prioritisation, and EV-miRNA analyses.

- `environment/`  
  Software environment information, including `sessionInfo.txt` and package versions.

## Public datasets and resources used

- GEO: `GSE56315`, `GSE30881`, `GSE130874`, `GSE31312`, `GSE171272`
- TCGA / Genomic Data Commons: `TCGA-DLBC`
- STRING
- CMap/LINCS (CLUE)
- Reactome
- KEGG

## High-level workflow

1. Preprocess human and canine expression datasets
2. Perform tumour-versus-normal differential expression
3. Map one-to-one human–dog orthologs
4. Define cross-species BROAD / STRICT / Tier1 modules
5. Perform Reactome/KEGG enrichment
6. Build STRING PPI network and identify hubs / bottlenecks
7. Integrate CMap/LINCS drug connectivity and target overlap
8. Compute module scores in validation cohorts
9. Perform canine and human survival analyses
10. Analyse EV-miRNA differential expression and validated target convergence on the PPI core
11. Generate final figures and manuscript tables

## Scripts used for key figures and tables

### Main figures
- `44_phase3E_crossSpecies_UMAP_BROAD.R` → Figure 1A
- `33_phase4_fig2A_GSE56315_module_scores_tumor_vs_normal.R` → Figure 1B
- `fig1c.R` → Figure 1C
- `26_phase7_figures_pathway_enrichment.R` → Figure 2A
- `27_phase7_figures_pathway_axes_summary.R` → Figure 2B
- `28_phase8_figures_PPI_centrality_and_drugs.R` → Figure 3A–3C
- `30_phase9_figures_module_scores_canine_cross_species_BROAD.R` → Figure 4A
- `32_phase7_fig6B_GSE130874_module_vs_hubs.R` → Figure 4B
- `34_phase8A_GSE130874_PFS_analysis.R` → Figure 5A
- `39_phase9C_GSE31312_OS_analysis.R` or `40_phase9C_GSE31312_OS_from_pdf.R` → Figure 5B
- `40_phase9D_crossCohort_survival_forest.R` → Figure 5C
- `43_phase10D_GSE171272_EVmiRNA_network_figure.R` → Figure 6A
- `43_phase10D_EVmiRNA_summary_dotplot.R` → Figure 6B

### Tables
- `45_phase9_Table1_cohort_overview.R` → Table 1
- `56_build_Table3_core_CMap_drugs_forManuscript.R` → Table 2
- `55_build_Table2_crossSpecies_modules_network_summary_from_tables.R` → Supplementary Table S1
- `56_fix_Table4_survival_associations_crossSpecies_module.R` → Supplementary Table S2
- `56_fix_Table5_EVmiRNAs_labels.R` → Supplementary Table S3

## Software environment

Analyses were run in **R 4.5.2** on Windows. See:
- `environment/sessionInfo.txt`
- `environment/package_versions.tsv`

## Reproducibility note

Raw public data files are not redistributed in this repository. All primary datasets were obtained from public repositories listed above. This repository instead provides the curated scripts and compact derived metadata required to document and reproduce the final analytical workflow reported in the manuscript.

## Contact

For questions regarding the repository or manuscript workflow, please contact the corresponding author.
