## =============================================================================
## 02 -- Broad-sense heritability, PopGroup added as a FIXED effect
## -----------------------------------------------------------------------------
## Directly addresses Reviewer 2, Major comment 1: "A fixed effect of
## population should also be added in your model (2) for heritability
## computation because it is not so relevant to inflate heritability
## estimates with effects of population structure."
##
## y = mu + PopGroup(fixed) + Genotype(random) + [design terms, random] + e
## H2 = sigma2_G / (sigma2_G + sigma2_e / r_bar), entry-mean basis, per
## environment x trait, exactly as in the original manuscript's Equation 3 --
## the changes are the added PopGroup fixed effect and the design terms:
## Aberystwyth fits Block+Row+Col (Block independently recorded); at
## Braunschweig, Block is reconstructed from Row, so only Block is fitted
## there (Row/Col would duplicate the same spatial information).
##
## Output: results/heritability/heritability_by_environment.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(asreml)
  library(dplyr)
})

d <- read.csv("../data/combined_phenotypes_long.csv", stringsAsFactors = FALSE)
d$Genotype <- factor(d$Genotype)
d$PopGroup <- factor(d$PopGroup)
d$Block <- factor(d$Block)
d$Row <- factor(d$Row)
d$Col <- factor(d$Col)

traits <- unique(d$trait)
envs <- unique(d$Environment)

fit_h2_one <- function(sub, has_row_col) {
  sub <- sub |> filter(!is.na(value))
  sub$Genotype <- droplevels(sub$Genotype)
  sub$PopGroup <- droplevels(sub$PopGroup)
  sub$Block <- droplevels(sub$Block)
  sub$Row <- droplevels(sub$Row)
  sub$Col <- droplevels(sub$Col)
  if (n_distinct(sub$Genotype) < 5 || nlevels(sub$PopGroup) < 2) return(NULL)

  ## Aberystwyth: Block+Row+Col (Block independently recorded, not
  ## confounded with Row/Col). Braunschweig: Block only (its Block is
  ## reconstructed from Row, so adding Row/Col would duplicate the same
  ## spatial information -- see 01_estimate_BLUEs.R and
  ## 00_prepare_combined_phenotypes.R for the reconstruction).
  random_form <- if (has_row_col) ~ Genotype + Block + Row + Col else ~ Genotype + Block
  fit <- tryCatch(
    asreml(fixed = value ~ PopGroup, random = random_form, data = sub,
           trace = FALSE, maxit = 30),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      asreml(fixed = value ~ PopGroup, random = ~ Genotype, data = sub,
             trace = FALSE, maxit = 30),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) return(NULL)
  for (i in 1:8) {
    if (isTRUE(fit$converge)) break
    fit <- tryCatch(update.asreml(fit), error = function(e) fit)
  }

  vc <- summary(fit)$varcomp
  if (!("Genotype" %in% rownames(vc))) return(NULL)
  sigma2_G <- max(vc["Genotype", "component"], 0)
  sigma2_e <- vc["units!R", "component"]
  r_bar <- mean(table(sub$Genotype))
  H2 <- sigma2_G / (sigma2_G + sigma2_e / r_bar)

  data.frame(sigma2_G = sigma2_G, sigma2_e = sigma2_e, r_bar = r_bar,
             H2 = H2, n_genotypes = n_distinct(sub$Genotype),
             converged = isTRUE(fit$converge))
}

results <- list()
for (env in envs) {
  has_row_col <- grepl("Aber-", env, fixed = TRUE)
  for (tr in traits) {
    sub <- d |> filter(Environment == env, trait == tr)
    if (nrow(sub) == 0) next
    out <- fit_h2_one(sub, has_row_col)
    if (is.null(out)) { cat(env, tr, ": FAILED (insufficient replication/groups)\n"); next }
    out$Environment <- env
    out$trait <- tr
    results[[paste(env, tr)]] <- out
    cat(env, tr, ": H2=", round(out$H2, 3), " r_bar=", round(out$r_bar, 2),
        " converged=", out$converged, "\n")
  }
}

h2_tab <- bind_rows(results) |>
  select(Environment, trait, n_genotypes, r_bar, sigma2_G, sigma2_e, H2, converged)
write.csv(h2_tab, "../results/heritability/heritability_by_environment.csv", row.names = FALSE)
cat("\nSaved: ../results/heritability/heritability_by_environment.csv\n")
print(h2_tab, digits = 3)
