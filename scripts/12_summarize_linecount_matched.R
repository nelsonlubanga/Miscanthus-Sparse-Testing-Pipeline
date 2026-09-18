## =============================================================================
## 12 -- Combine and summarize the line-count-matched analyses (k=2,3,4)
## -----------------------------------------------------------------------------
## Reads the three per-k outputs of 11_fit_linecount_matched_designs.R and
## produces a combined summary table, for direct comparison against the
## original (confounded) D1-D6 sweep (script 05), which correlates
## unique-line-count with predictive ability at r=0.985. The corresponding
## figure is built from this table by script 13.
##
## Output: results/predictions/linecount_matched_supplementary/summary_table.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages(library(dplyr))

out_dir <- "../results/predictions/linecount_matched_supplementary"
dir.create(out_dir, showWarnings = FALSE)

all_k <- bind_rows(lapply(2:4, function(k) {
  f <- sprintf("../results/predictions/linecount_matched_predictive_ability_k%d.csv", k)
  read.csv(f, stringsAsFactors = FALSE)
}))

## Per-design mean PA, per k
by_design <- all_k |>
  group_by(k_envs, design, s, n, unique_lines, total_plots) |>
  summarise(mean_PA = mean(predictive_ability, na.rm = TRUE),
            sd_PA = sd(predictive_ability, na.rm = TRUE), .groups = "drop") |>
  arrange(k_envs, s)

## Correlation between overlap (s) and PA, per k (design-level, n=5 each)
cor_by_k <- by_design |>
  group_by(k_envs) |>
  summarise(r_overlap_vs_PA = cor(s, mean_PA), .groups = "drop")

summary_table <- by_design |> left_join(cor_by_k, by = "k_envs")
write.csv(summary_table, file.path(out_dir, "summary_table.csv"), row.names = FALSE)

cat("=== Summary table ===\n")
print(as.data.frame(summary_table))
cat("\n=== Correlation (overlap vs PA at fixed line-count), by k, vs. original confounded r=0.985 ===\n")
print(as.data.frame(cor_by_k))
cat("\nSaved:", file.path(out_dir, "summary_table.csv"), "\n")
