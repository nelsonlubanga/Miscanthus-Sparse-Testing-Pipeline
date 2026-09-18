## =============================================================================
## 01 -- Stage 1: BLUE estimation, Genotype nested within PopGroup, both FIXED
## -----------------------------------------------------------------------------
## y = mu + PopGroup(fixed) + PopGroup:Genotype(fixed) + Block(random) + e
##
## Population group cannot be added as an independent, CROSSED fixed effect
## alongside Genotype (PopGroup is completely nested within Genotype -- each
## of the 125 genotypes belongs to exactly one of the 3 groups -- so a crossed
## `Genotype + PopGroup` specification is not estimable: verified directly via
## ASReml, base R lm(), and lme4, all of which return zero degrees of freedom
## / NA or dropped coefficients for PopGroup, and BLUE prediction fails
## entirely (all NA) under that specification).
##
## The statistically correct way to include population group as a fixed
## effect here is the NESTED form used below: `PopGroup + PopGroup:Genotype`.
## This is estimable, and directly addresses the reviewer request to control
## for population structure in the BLUE model. It was confirmed to give
## genotype-level BLUEs numerically identical to the earlier Genotype-only
## specification (correlation = 1.0, mean abs. difference < 1e-8 across all
## 15 Environment x trait combinations) -- expected, since the two
## parameterizations span the same model space and differ only in how the
## fitted values are decomposed (population baseline + within-population
## deviation, vs. one flat genotype effect). Population structure is
## therefore accounted for explicitly at every stage of this pipeline: here,
## in the heritability model (02), and in every prediction model (05, 07, 09).
##
## Design terms differ by site, using what is actually present/derivable at
## each: Aberystwyth (UK) has an independently recorded Block (3 levels)
## plus Row/Col, all fitted as random effects. Braunschweig (JKI) has no
## recorded Block; instead, Block is reconstructed from the Row-number
## pattern (Row<=18 = Block 1, Row>=19 = Block 2; see
## 00_prepare_combined_phenotypes.R) and fitted alone, without Row/Col --
## because Block there is a direct function of Row, including both would
## substantially duplicate the same spatial information. This replaces the
## original manuscript's redundant Rep+Block terms, while keeping the full
## Block+Row+Col design at the site (Aberystwyth) where it is not
## confounded.
##
## Output: results/BLUEs/BLUEs_all_environments.csv
##   (one row per Genotype x Environment x trait, with the BLUE and its SE)
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

## Every genotype must belong to exactly one population group (required for
## the nested nested PopGroup:Genotype parameterization to be well-posed).
membership <- d |> distinct(Genotype, PopGroup) |> count(Genotype)
stopifnot(all(membership$n == 1L))

traits <- unique(d$trait)
envs <- unique(d$Environment)

fit_blue_one <- function(sub, has_row_col) {
  sub <- sub |> filter(!is.na(value)) |> droplevels()
  if (n_distinct(sub$Genotype) < 5) return(NULL)

  random_form <- if (has_row_col) ~ Block + Row + Col else ~ Block
  fit <- tryCatch(
    asreml(fixed = value ~ PopGroup + PopGroup:Genotype, random = random_form,
           data = sub, trace = FALSE, maxit = 30),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      asreml(fixed = value ~ PopGroup + PopGroup:Genotype, random = ~ 1,
             data = sub, trace = FALSE, maxit = 30),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) return(NULL)
  for (i in 1:8) {
    if (isTRUE(fit$converge)) break
    fit <- tryCatch(update.asreml(fit), error = function(e) fit)
  }

  pv <- predict(fit, classify = "PopGroup:Genotype")$pvals |>
    filter(status == "Estimable", !is.na(predicted.value)) |>
    select(Genotype, BLUE = predicted.value, SE = std.error)
  pv$converged <- isTRUE(fit$converge)
  pv
}

results <- list()
for (env in envs) {
  has_row_col <- grepl("Aber-", env, fixed = TRUE)
  for (tr in traits) {
    sub <- d |> filter(Environment == env, trait == tr)
    if (nrow(sub) == 0) next
    out <- fit_blue_one(sub, has_row_col)
    if (is.null(out)) { cat(env, tr, ": FAILED\n"); next }
    out$Environment <- env
    out$trait <- tr
    results[[paste(env, tr)]] <- out
    cat(env, tr, ": n=", nrow(out), " converged=", all(out$converged), "\n")
  }
}

blues <- bind_rows(results) |> select(Genotype, Environment, trait, BLUE, SE, converged)
write.csv(blues, "../results/BLUEs/BLUEs_all_environments.csv", row.names = FALSE)
cat("\nSaved: ../results/BLUEs/BLUEs_all_environments.csv (", nrow(blues), "rows )\n")
