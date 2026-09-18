## =============================================================================
## 09 -- Variance components for M1/M2/M3 (PopGroup fixed), fit once per
##        trait on the full (unmasked) data -- mirrors the original
##        manuscript's Figure 2 variance-decomposition analysis, but with
##        PopGroup now included as a fixed effect throughout.
## -----------------------------------------------------------------------------
## Extracts sigma_E^2 (ENV), sigma_L^2 (LINE), sigma_g^2 (GENO, M2/M3),
## sigma_gE^2 (GXE, M3) and sigma_e^2 (residual) directly from BGLR's
## reported variance components (ETA[[k]]$varB for BRR terms,
## ETA[[k]]$varU for RKHS terms, fit$varE for the residual).
##
## Output: results/variance_components/variance_components.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(BGLR)
  library(parallel)
})
N_CORES <- 4

dir.create("../results/variance_components", showWarnings = FALSE)

TRAITS <- c("DM", "FW", "MC")
ENVS <- c("Aber-2015", "Aber-2016", "Brau-2014", "Brau-2015", "Brau-2016")

blues <- read.csv("../results/BLUEs/BLUEs_all_environments.csv", stringsAsFactors = FALSE)
load("../data/markers/Gmatrix_125.rda")
popgroups <- read.csv("../data/popgroups_125.csv", stringsAsFactors = FALSE)

GENOTYPES <- rownames(G)
pop_lookup <- setNames(popgroups$PopGroup, popgroups$Genotype)

obs <- expand.grid(Genotype = GENOTYPES, Environment = ENVS, stringsAsFactors = FALSE)
obs$PopGroup <- pop_lookup[obs$Genotype]

Zg <- model.matrix(~ factor(obs$Genotype, levels = GENOTYPES) - 1)
Ze <- model.matrix(~ factor(obs$Environment, levels = ENVS) - 1)
Xpop_full <- model.matrix(~ factor(obs$PopGroup))
Xpop <- Xpop_full[, -1, drop = FALSE]

K_obs <- Zg %*% G %*% t(Zg)
K_gxe <- K_obs * (Ze %*% t(Ze))

build_y_full <- function(trait) {
  sub <- blues |> filter(trait == !!trait) |> select(Genotype, Environment, BLUE)
  wide <- obs |> select(Genotype, Environment) |>
    left_join(sub, by = c("Genotype", "Environment"))
  wide$BLUE
}

fit_full_model <- function(y, model_name, n_iter = 8000, burn_in = 2000, thin = 5) {
  ETA <- list(
    POP  = list(X = Xpop, model = "FIXED"),
    ENV  = list(X = Ze,   model = "BRR"),
    LINE = list(X = Zg,   model = "BRR")
  )
  if (model_name %in% c("M2", "M3")) ETA$GENO <- list(K = K_obs, model = "RKHS")
  if (model_name == "M3")            ETA$GXE  <- list(K = K_gxe, model = "RKHS")
  save_dir <- tempfile("bglr_vc_")
  dir.create(save_dir)
  fit <- BGLR(y = y, ETA = ETA, nIter = n_iter, burnIn = burn_in, thin = thin,
              verbose = FALSE, saveAt = file.path(save_dir, "bg_"))
  unlink(save_dir, recursive = TRUE)
  fit
}

extract_vc <- function(fit, model_name) {
  vc <- c(
    sigma2_E   = fit$ETA$ENV$varB,
    sigma2_L   = fit$ETA$LINE$varB,
    sigma2_G   = if (model_name %in% c("M2", "M3")) fit$ETA$GENO$varU else NA_real_,
    sigma2_GxE = if (model_name == "M3") fit$ETA$GXE$varU else NA_real_,
    sigma2_e   = fit$varE
  )
  total <- sum(vc, na.rm = TRUE)
  data.frame(model = model_name, t(vc), total_variance = total,
             pct_E = 100 * vc["sigma2_E"] / total,
             pct_e = 100 * vc["sigma2_e"] / total,
             row.names = NULL)
}

N_REPS <- 5
tasks <- expand.grid(trait = TRAITS, model = c("M1", "M2", "M3"), rep = seq_len(N_REPS),
                      stringsAsFactors = FALSE)
cat("Total tasks:", nrow(tasks), " | parallel workers:", N_CORES, "\n")

y_full_by_trait <- setNames(lapply(TRAITS, build_y_full), TRAITS)

task_results <- mclapply(seq_len(nrow(tasks)), function(i) {
  tr <- tasks$trait[i]; m <- tasks$model[i]; r <- tasks$rep[i]
  fit <- fit_full_model(y_full_by_trait[[tr]], m)
  vc <- extract_vc(fit, m)
  vc$trait <- tr
  vc$rep <- r
  vc
}, mc.cores = N_CORES)

failed <- sapply(task_results, function(x) inherits(x, "try-error"))
cat("Failed tasks:", sum(failed), "of", length(task_results), "\n")

raw <- bind_rows(task_results[!failed]) |>
  select(trait, model, rep, sigma2_E, sigma2_L, sigma2_G, sigma2_GxE, sigma2_e,
         total_variance, pct_E, pct_e)
write.csv(raw, "../results/variance_components/variance_components_raw_reps.csv", row.names = FALSE)

out <- raw |>
  group_by(trait, model) |>
  summarise(across(c(sigma2_E, sigma2_L, sigma2_G, sigma2_GxE, sigma2_e,
                      total_variance, pct_E, pct_e), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")
write.csv(out, "../results/variance_components/variance_components.csv", row.names = FALSE)

cat("\n=== Relative to M1: residual and environmental variance reduction (%), averaged over", N_REPS, "reps ===\n")
for (tr in TRAITS) {
  sub <- out[out$trait == tr, ]
  m1e <- sub$sigma2_e[sub$model == "M1"]; m3e <- sub$sigma2_e[sub$model == "M3"]
  m1E <- sub$sigma2_E[sub$model == "M1"]; m3E <- sub$sigma2_E[sub$model == "M3"]
  cat(tr, ": residual variance change M1->M3 =", round(100 * (m3e - m1e) / m1e, 1),
      "%  |  environmental variance change M1->M3 =", round(100 * (m3E - m1E) / m1E, 1), "%\n")
}

cat("\n=== Per-rep reduction estimates (to gauge stability) ===\n")
wide_e <- raw |> select(trait, model, rep, sigma2_e) |>
  pivot_wider(names_from = model, values_from = sigma2_e)
wide_E <- raw |> select(trait, model, rep, sigma2_E) |>
  pivot_wider(names_from = model, values_from = sigma2_E)
for (tr in TRAITS) {
  re <- wide_e[wide_e$trait == tr, ]; rE <- wide_E[wide_E$trait == tr, ]
  cat(tr, "residual reduction per rep (%):", round(100 * (re$M3 - re$M1) / re$M1, 1), "\n")
  cat(tr, "environmental reduction per rep (%):", round(100 * (rE$M3 - rE$M1) / rE$M1, 1), "\n")
}

cat("\nSaved: ../results/variance_components/variance_components.csv (mean of", N_REPS, "reps)\n")
cat("Saved: ../results/variance_components/variance_components_raw_reps.csv (all reps)\n")
