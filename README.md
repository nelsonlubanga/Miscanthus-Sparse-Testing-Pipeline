# Miscanthus Sparse-Testing Pipeline

Reproducible analysis pipeline for "Sparse testing and multi-environment
genomic prediction improve selection efficiency in Miscanthus breeding."
Covers phenotype preparation, BLUE/heritability estimation, genomic
relationship matrix construction, sparse-testing designs (D1-D6), CV1/CV2
cross-validation, variance-component decomposition, and the line-count-matched
supplementary analysis (Reviewer 2, Major Comment 3).

## Requirements

R packages: `readxl`, `dplyr`, `tidyr`, `stringr`, `rrBLUP`, `adegenet`,
`BGLR`, `parallel`, `ggplot2`, `asreml`.

`asreml` (ASReml-R 4.2, used in scripts 01-02) requires a commercial license
from VSNi and is not on CRAN. All other packages are open source.

Scripts use `mclapply` for parallelism (fork-based; Linux/macOS only) and set
`N_CORES <- 4` by default -- adjust to the available machine.

Each script begins with `setwd()` to an absolute path; edit this to match
where the pipeline is located before running.

## Required raw inputs (not included in this folder)

Script `00_prepare_combined_phenotypes.R` expects two raw files one directory
above this pipeline folder:

- `Original phenotypes.xlsx` (sheet `ABR33Data_ForNelson`) -- Aberystwyth
- `JKI_Ger.csv` -- Braunschweig

Place both alongside (one level above) this folder, or edit `RAW_UK_XLSX`
and `RAW_JKI_CSV` at the top of script 00 to point to their actual location.

## Directory structure

- `data/` -- inputs and intermediate data (marker matrix, genomic relationship
  matrix, population-group assignments, sparse-testing design definitions)
- `scripts/` -- numbered pipeline scripts, run in order (see below)
- `results/` -- all outputs (BLUEs, heritability, predictions, figures,
  variance components)

## Run order

| Script | Produces | Manuscript reference |
|---|---|---|
| `00_prepare_combined_phenotypes.R` | `data/combined_phenotypes_long.csv` | -- |
| `00a_derive_population_groups.R` | `data/popgroups_125.csv` | population structure (DAPC) |
| `01_estimate_BLUEs.R` | `results/BLUEs/BLUEs_all_environments.csv` | Equation 1 |
| `02_estimate_heritability.R` | `results/heritability/heritability_by_environment.csv` | Equation 2, Table with H² |
| `03_generate_sparse_testing_designs.R` | `data/sparse_testing_designs.rds` | Table 1 (D1-D6), CV1/CV2 folds |
| `04_build_Gmatrix.R` | `data/markers/Gmatrix_125.rda` | Equation 6 (G matrix) |
| `05_fit_sparse_testing_models.R` | `results/predictions/predictive_ability.csv` | Figure 4 |
| `06_plot_predictive_ability.R` | `results/predictions/predictive_ability_boxplot.png` | Figure 4 |
| `07_fit_cv1_cv2_models.R` | `results/predictions/predictive_ability_cv1_cv2.csv` | Figure 3a/3b |
| `08_plot_cv1_cv2.R` | `..._CV1_boxplot.png`, `..._CV2_boxplot.png` | Figure 3a/3b |
| `09_variance_components.R` | `results/variance_components/variance_components.csv` | Figure 2 |
| `10_plot_variance_components.R` | `results/variance_components/variance_components_stacked.png` | Figure 2 |
| `11_fit_linecount_matched_designs.R` | `linecount_matched_predictive_ability_k<K_ENVS>.csv` | Supplementary Figure/Table S1 |
| `12_summarize_linecount_matched.R` | `linecount_matched_supplementary/summary_table.csv` | Supplementary Table S1 |
| `13_prepare_supplementary_materials.R` | `Supplementary_Figure_S1.png`, `Supplementary_Table_S1.csv` | Supplementary Figure/Table S1 |

Script 11 takes `K_ENVS` (2, 3, or 4) as a parameter at the top of the file
and must be run once per value before scripts 12-13.

Scripts 03-13 depend on the outputs of earlier scripts in this table and
must be run in order on a first pass.
