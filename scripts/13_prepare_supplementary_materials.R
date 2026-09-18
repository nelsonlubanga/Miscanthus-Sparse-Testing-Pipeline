## =============================================================================
## 13 -- Prepare final supplementary table and figure for the line-count-matched
##        analysis (k=2,3,4), addressing Reviewer 2 Major Comment 3.
## -----------------------------------------------------------------------------
## Builds one multi-panel supplementary figure (facet_grid: k rows x trait
## columns) directly from the per-observation data (scripts 11/12), and
## formats summary_table.csv (script 12) into a clean, publication-ready
## supplementary table with %NOL.
##
## Output: results/predictions/linecount_matched_supplementary/
##   Supplementary_Figure_S1.png
##   Supplementary_Table_S1.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(ggplot2)
})

out_dir <- "../results/predictions/linecount_matched_supplementary"

## -----------------------------------------------------------------------------
## Supplementary Figure: one combined boxplot, k (rows) x trait (columns)
## -----------------------------------------------------------------------------
all_k <- bind_rows(lapply(2:4, function(k) {
  f <- sprintf("../results/predictions/linecount_matched_predictive_ability_k%d.csv", k)
  read.csv(f, stringsAsFactors = FALSE)
}))

all_k$pct_NOL <- round(100 * (all_k$unique_lines - all_k$s) / all_k$unique_lines)
all_k$pct_NOL <- factor(all_k$pct_NOL, levels = sort(unique(all_k$pct_NOL), decreasing = TRUE))
all_k$trait <- factor(all_k$trait, levels = c("DM", "FW", "MC"))
all_k$model <- factor(all_k$model, levels = c("M1", "M2", "M3"))
all_k$k_label <- factor(paste0("k = ", all_k$k_envs, " environments"),
                         levels = paste0("k = ", 2:4, " environments"))

p <- ggplot(all_k, aes(x = pct_NOL, y = predictive_ability, fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.65,
               outlier.size = 0.5, linewidth = 0.3) +
  facet_grid(k_label ~ trait) +
  scale_fill_manual(values = c(M1 = "#F08A7E", M2 = "#4FAE62", M3 = "#5B9BD5"),
                     name = "Model") +
  labs(x = "Percentage of non-overlapping lines (unique training-line count held fixed within each k)",
       y = "Predictive ability") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey85"),
        panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave(file.path(out_dir, "Supplementary_Figure_S1.png"), p, width = 9.5, height = 8.5, dpi = 300)
cat("Saved:", file.path(out_dir, "Supplementary_Figure_S1.png"), "\n")

## -----------------------------------------------------------------------------
## Supplementary Table: clean, publication-formatted version of summary_table.csv
## -----------------------------------------------------------------------------
summ <- read.csv(file.path(out_dir, "summary_table.csv"), stringsAsFactors = FALSE)

tab <- summ |>
  mutate(
    `Environments (k)` = k_envs,
    Design = design,
    `Shared lines (s)` = s,
    `Lines unique per environment (n)` = n,
    `Unique training lines` = unique_lines,
    `Total plots` = total_plots,
    `% NOL` = round(100 * (unique_lines - s) / unique_lines),
    `Mean predictive ability` = round(mean_PA, 3),
    `SD` = round(sd_PA, 3),
    `Correlation (overlap vs. PA, fixed line count)` = round(r_overlap_vs_PA, 2)
  ) |>
  select(`Environments (k)`, Design, `Shared lines (s)`, `Lines unique per environment (n)`,
         `Unique training lines`, `Total plots`, `% NOL`, `Mean predictive ability`, `SD`,
         `Correlation (overlap vs. PA, fixed line count)`) |>
  arrange(`Environments (k)`, desc(`% NOL`))

write.csv(tab, file.path(out_dir, "Supplementary_Table_S1.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "Supplementary_Table_S1.csv"), "\n")
print(tab, row.names = FALSE)
