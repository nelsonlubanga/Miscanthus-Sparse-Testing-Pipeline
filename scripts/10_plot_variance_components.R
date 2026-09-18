## =============================================================================
## 10 -- 100% stacked barplot of variance components (E, G, GxE, L, R) by
##        model (M1/M2/M3) and trait, matching the manuscript's Figure 2
##        format. Source: results/variance_components/variance_components.csv
##        (revised analysis, PopGroup fitted as a fixed effect throughout).
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

vc <- read.csv("../results/variance_components/variance_components.csv", stringsAsFactors = FALSE)

long <- vc |>
  select(trait, model, E = sigma2_E, G = sigma2_G, GxE = sigma2_GxE,
         L = sigma2_L, R = sigma2_e) |>
  pivot_longer(cols = c(E, G, GxE, L, R), names_to = "component", values_to = "variance") |>
  filter(!is.na(variance)) |>
  group_by(trait, model) |>
  mutate(pct = 100 * variance / sum(variance)) |>
  ungroup()

long$trait <- factor(long$trait, levels = c("DM", "FW", "MC"))
long$model <- factor(long$model, levels = c("M1", "M2", "M3"))
long$component <- factor(long$component, levels = c("R", "L", "G", "GxE", "E"))

long <- long |>
  arrange(trait, model, component) |>
  group_by(trait, model) |>
  mutate(ymax = pmin(cumsum(pct), 100), ymin = pmax(ymax - pct, 0)) |>
  ungroup()

w <- 0.35
long$xnum <- as.numeric(long$model)

p <- ggplot(long) +
  geom_rect(aes(xmin = xnum - w, xmax = xnum + w, ymin = ymin, ymax = ymax, fill = component),
            colour = "white", linewidth = 0.3) +
  scale_x_continuous(breaks = 1:3, labels = c("M1", "M2", "M3")) +
  facet_wrap(~ trait, nrow = 1) +
  scale_fill_manual(
    values = c(E = "#F4867A", G = "#9B9B2E", GxE = "#4CAF6D", L = "#4FA3D9", R = "#E066D6"),
    breaks = c("E", "G", "GxE", "L", "R"),
    name = "Variance component"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage variance") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey85"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("../results/variance_components/variance_components_stacked.png", p,
       width = 8.5, height = 4.5, dpi = 300)
cat("Saved: ../results/variance_components/variance_components_stacked.png\n")
