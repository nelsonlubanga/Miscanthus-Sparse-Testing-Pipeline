## Variance components for M1/M2/M3 with the corrected G and population group fitted as a
## random (BRR) effect, so its share of variance can be shown. Otherwise identical to script 12.
## Output: reviewer_FW_MC_pipeline/results/corrected_G/variance_components_popgroup_random/

rm(list = ls())
setwd(file.path("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus Sparse Testing Pipeline","scripts"))
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(BGLR)
  library(parallel)
})
N_CORES <- 10

OUT <- "../reviewer_FW_MC_pipeline/results/corrected_G/variance_components_popgroup_random"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c("DM", "FW", "MC")
ENVS <- c("Aber-2015", "Aber-2016", "Brau-2014", "Brau-2015", "Brau-2016")

blues <- read.csv("../results/BLUEs/BLUEs_all_environments.csv", stringsAsFactors = FALSE)
load("../data/markers/Gmatrix_125_corrected.rda")
popgroups <- read.csv("../data/popgroups_125.csv", stringsAsFactors = FALSE)

GENOTYPES <- rownames(G)
pop_lookup <- setNames(popgroups$PopGroup, popgroups$Genotype)

obs <- expand.grid(Genotype = GENOTYPES, Environment = ENVS, stringsAsFactors = FALSE)
obs$PopGroup <- pop_lookup[obs$Genotype]

Zg <- model.matrix(~ factor(obs$Genotype, levels = GENOTYPES) - 1)
Ze <- model.matrix(~ factor(obs$Environment, levels = ENVS) - 1)
Xpop <- model.matrix(~ factor(obs$PopGroup) - 1)

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
    POP  = list(X = Xpop, model = "BRR"),
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
    sigma2_PopGroup = fit$ETA$POP$varB,
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
  select(trait, model, rep, sigma2_PopGroup, sigma2_E, sigma2_L, sigma2_G, sigma2_GxE, sigma2_e,
         total_variance, pct_E, pct_e)
write.csv(raw, file.path(OUT, "variance_components_raw_reps.csv"), row.names = FALSE)

out <- raw |>
  group_by(trait, model) |>
  summarise(across(c(sigma2_PopGroup, sigma2_E, sigma2_L, sigma2_G, sigma2_GxE, sigma2_e,
                      total_variance, pct_E, pct_e), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")
write.csv(out, file.path(OUT, "variance_components.csv"), row.names = FALSE)

cat("Saved:", OUT, "\n")
