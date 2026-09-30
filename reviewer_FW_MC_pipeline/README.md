# Within-population-group analysis for FW and MC, and robustness checks

Extends `reviewer_DM_pipeline` (objective iii; DM) to fresh weight (FW) and
moisture content (MC), and adds the checks reported in the revised Methods.
Scripts take the trait and/or an output folder as command-line arguments.

| Script | Purpose |
|---|---|
| `00_pooled_BLUEs_FW_MC.R` | One BLUE per genotype per trial for FW and MC (ASReml; Equations 9a/9b) |
| `01_run_trait_reviewer_analysis.R` | CV1/CV2 with population group unadjusted/adjusted phenotypes (`FW`/`MC`) |
| `02_validate_results.R` | Checks masks, training-only adjustment, scores |
| `03`–`05` | Figures: within-group predictive ability (Figures 5–6 in the manuscript) |
| `06_chain_variability.R` | Monte Carlo error: extra chains for M2/M3 (all traits) |
| `07_plot_BLUE_correlations.R` | Site-year BLUE correlation matrix (supplementary figure) |
| `08_representative_chain_check.R` | Rank-normalised split R-hat for FW and MC (8,000 and 80,000 iterations) |
| `09_build_corrected_G.R` | G from correctly coded (-1/0/1) markers: `data/markers/Gmatrix_125_corrected.rda` |
| `10_corrected_G_impact.R`, `11_corrected_G_siteyear_check.R` | Effect of the corrected G on predictive ability |
| `12_variance_components_corrected_G.R`, `13_plot_variance_corrected_G.R` | Variance components and Figure 2 with the corrected G |

MCMC checkpoints and chain-check `.rds` files are not tracked (see `.gitignore`);
rerun the scripts to regenerate them.
