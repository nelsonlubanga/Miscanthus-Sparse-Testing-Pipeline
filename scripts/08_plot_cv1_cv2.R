## =============================================================================
## 08 -- Plot predictive ability under CV1 and CV2 (manuscript Figure 3a/3b),
##        as boxplots (distribution across the 10 reps x 5 folds) instead of
##        barplots. CV1 and CV2 are saved as two SEPARATE figures, matching
##        the manuscript's Figure 3a / Figure 3b split.
## -----------------------------------------------------------------------------
## X-axis: environment. Facets: trait. Fill: model (M1 = E+L, M2 = E+L+G,
## M3 = E+L+G+GxE; all with PopGroup as a fixed effect).
##
## Output: results/predictions/predictive_ability_CV1_boxplot.png
##         results/predictions/predictive_ability_CV2_boxplot.png
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(ggplot2)
})

pa <- read.csv("../results/predictions/predictive_ability_cv1_cv2.csv", stringsAsFactors = FALSE)
pa <- pa |> filter(!is.na(predictive_ability))

env_short <- c("Aber-2015" = "Aber15", "Aber-2016" = "Aber16",
               "Brau-2014" = "Brau14", "Brau-2015" = "Brau15", "Brau-2016" = "Brau16")
pa$Environment <- factor(env_short[pa$Environment], levels = unname(env_short))
pa$trait <- factor(pa$trait, levels = c("DM", "FW", "MC"))
pa$model <- factor(pa$model, levels = c("M1", "M2", "M3"))

plot_one_cv <- function(cv_scheme, out_file) {
  sub <- pa |> filter(CV == cv_scheme)
  p <- ggplot(sub, aes(x = Environment, y = predictive_ability, fill = model)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                 outlier.size = 0.6, linewidth = 0.3) +
    facet_wrap(~ trait, nrow = 1) +
    scale_fill_manual(values = c(M1 = "#F08A7E", M2 = "#4FAE62", M3 = "#5B9BD5"),
                       name = "Model") +
    labs(x = NULL, y = "Predictive ability") +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey85"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
  ggsave(out_file, p, width = 9, height = 3.4, dpi = 300)
  cat("Saved:", out_file, "\n")
}

plot_one_cv("CV1", "../results/predictions/predictive_ability_CV1_boxplot.png")
plot_one_cv("CV2", "../results/predictions/predictive_ability_CV2_boxplot.png")

## --- summary tables ---
agg <- pa |>
  group_by(CV, trait, model) |>
  summarise(mean_PA = round(mean(predictive_ability), 3),
            min_PA = round(min(predictive_ability), 3),
            max_PA = round(max(predictive_ability), 3), .groups = "drop")
cat("\n=== Mean / range predictive ability by CV x trait x model ===\n")
print(as.data.frame(agg))

cv_overall <- pa |> group_by(CV) |>
  summarise(mean_PA = round(mean(predictive_ability), 3), .groups = "drop")
cat("\n=== Overall mean by CV scheme ===\n")
print(as.data.frame(cv_overall))
