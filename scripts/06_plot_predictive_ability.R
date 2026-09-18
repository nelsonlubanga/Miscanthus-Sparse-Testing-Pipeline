## =============================================================================
## 06 -- Plot predictive ability across sparse-testing designs (D1-D6),
##        matching the manuscript's original figure format but as boxplots
##        (distribution across the 10 reps) instead of barplots.
## -----------------------------------------------------------------------------
## Facets: trait (DM, FW, MC). X-axis: percentage of non-overlapping (NOL)
## genotypes per environment (100, 80, 60, 40, 20, 0 <- D1..D6). Fill: model
## (M1 = E+L, M2 = E+L+G, M3 = E+L+G+GxE; all with PopGroup as a fixed effect).
##
## Output: results/predictions/predictive_ability_boxplot.png
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(ggplot2)
})

pa <- read.csv("../results/predictions/predictive_ability.csv", stringsAsFactors = FALSE)

nol_pct <- c(D1 = 100, D2 = 80, D3 = 60, D4 = 40, D5 = 20, D6 = 0)
pa$NOL_pct <- factor(nol_pct[pa$design], levels = c(100, 80, 60, 40, 20, 0))
pa$trait <- factor(pa$trait, levels = c("DM", "FW", "MC"))
pa$model <- factor(pa$model, levels = c("M1", "M2", "M3"))

p <- ggplot(pa, aes(x = NOL_pct, y = predictive_ability, fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.65,
               outlier.size = 0.8, linewidth = 0.3) +
  facet_wrap(~ trait, nrow = 1) +
  scale_fill_manual(values = c(M1 = "#F08A7E", M2 = "#4FAE62", M3 = "#5B9BD5"),
                     name = "Model") +
  labs(x = "Percentage of non-overlapping genotypes per environment",
       y = "Predictive ability") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey85"),
        panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave("../results/predictions/predictive_ability_boxplot.png", p,
       width = 9, height = 3.2, dpi = 300)
cat("Saved: ../results/predictions/predictive_ability_boxplot.png\n")
