## Figure 1: (a) PCA of the 125 genotypes (existing panel); (b) genotype BLUEs by population group,
## per trait and environment, with Tukey HSD letters from one-way ANOVA (BLUE ~ population group).
## Run from the pipeline root: Rscript scripts/14_figure1_pca_and_blues.R
suppressMessages({library(dplyr); library(ggplot2); library(patchwork); library(multcompView)})
blues <- read.csv("results/BLUEs/BLUEs_all_environments.csv", stringsAsFactors = FALSE)
pg <- read.csv("data/popgroups_125.csv", stringsAsFactors = FALSE)
d <- blues |> inner_join(pg, by = "Genotype") |> filter(is.finite(BLUE))
d$PopGroup <- factor(d$PopGroup, levels = c("M. sacchariflorus", "M. sinensis", "M. x gig"),
                     labels = c("M. sacchariflorus", "M. sinensis", "Interspecific hybrids"))
d$Environment <- sub("^(Aber|Brau)-20", "\\1", d$Environment)
d$trait <- factor(d$trait, levels = c("DM", "FW", "MC"), labels = c("DM (kg)", "FW (kg)", "MC (%)"))

tests <- d |> group_by(trait, Environment) |> group_modify(function(s, k) {
  fit <- aov(BLUE ~ PopGroup, data = s)
  lt <- multcompLetters4(fit, TukeyHSD(fit))$PopGroup$Letters
  data.frame(PopGroup = names(lt), letter = unname(lt), P = summary(fit)[[1]][["Pr(>F)"]][1],
             y = max(s$BLUE) + 0.08 * diff(range(s$BLUE)))
}) |> ungroup()
tests$PopGroup <- factor(tests$PopGroup, levels = levels(d$PopGroup))
write.csv(tests, "results/figure1/Figure1b_ANOVA_Tukey.csv", row.names = FALSE)

b <- ggplot(d, aes(PopGroup, BLUE, fill = PopGroup)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.45) +
  geom_text(data = tests, aes(PopGroup, y, label = letter), inherit.aes = FALSE, size = 3) +
  facet_grid(trait ~ Environment, scales = "free_y") +
  scale_x_discrete(labels = c(expression(italic("M. sacchariflorus")), expression(italic("M. sinensis")), "Interspecific hybrids")) +
  scale_fill_manual(values = c("grey35", "#2ca02c", "#e41a1c")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Genotype BLUE") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey90"),
        axis.text.x = element_text(angle = 30, hjust = 1))
a <- wrap_elements(full = grid::rasterGrob(png::readPNG("results/figure1/Figure1a_PCA.png")))
p <- (a / b) + plot_layout(heights = c(1, 1.3)) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave("results/figure1/Figure1_PCA_and_BLUEs.png", p, width = 8, height = 11, dpi = 300, bg = "white")
ggsave("results/figure1/Figure1_PCA_and_BLUEs.pdf", p, width = 8, height = 11, bg = "white")
ggsave("results/figure1/Figure1b_BLUEs_by_population_group.png", b, width = 9, height = 6.5, dpi = 300, bg = "white")
print(tests |> distinct(trait, Environment, P) |> mutate(P = signif(P, 2)), n = 20)
