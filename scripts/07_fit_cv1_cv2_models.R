## =============================================================================
## 07 -- Predictive ability under CV1 and CV2 (manuscript Figure 3a/3b),
##        with PopGroup as a FIXED effect throughout.
## -----------------------------------------------------------------------------
## Distinct from 05_fit_sparse_testing_models.R (the D1-D6 NOL/OL sparse-
## testing design sweep). Here, per the manuscript's "Model assessment using
## cross-validation schemes" section: "A random five-fold partitioning of the
## entire population was used in both scenarios, and predictive ability was
## calculated as the correlation between predicted and observed values within
## the same environment."
##
## CV1 (Table 2a): mask an entire genotype across ALL 5 environments --
##   5-fold partition of the 125 genotypes; predicts newly developed material.
## CV2 (Table 2b): mask individual genotype x environment CELLS -- 5-fold
##   partition of the actual observed genotype x environment records;
##   predicts sparse testing (genotype seen in >=1 other environment).
##
## Repeated 10 times (independent random fold assignments) for both schemes,
## consistent with the repetition used elsewhere in this pipeline, so that
## each Environment x trait x model x CV combination yields a distribution
## (10 reps x 5 folds = 50 values) suitable for a boxplot.
##
## Output: results/predictions/predictive_ability_cv1_cv2.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(BGLR)
  library(parallel)
})

N_REPS <- 10
N_FOLDS <- 5
N_CORES <- 4
TRAITS <- c("DM", "FW", "MC")
ENVS <- c("Aber-2015", "Aber-2016", "Brau-2014", "Brau-2015", "Brau-2016")

## -----------------------------------------------------------------------------
## Shared setup (same construction as 05_fit_sparse_testing_models.R)
## -----------------------------------------------------------------------------
blues <- read.csv("../results/BLUEs/BLUEs_all_environments.csv", stringsAsFactors = FALSE)
load("../data/markers/Gmatrix_125.rda")        # G, M_centered
popgroups <- read.csv("../data/popgroups_125.csv", stringsAsFactors = FALSE)

GENOTYPES <- rownames(G)
stopifnot(length(GENOTYPES) == 125)
pop_lookup <- setNames(popgroups$PopGroup, popgroups$Genotype)
stopifnot(all(GENOTYPES %in% names(pop_lookup)))

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

fit_one_model <- function(y_train, model_name, n_iter = 1500, burn_in = 300, thin = 3) {
  ETA <- list(
    POP  = list(X = Xpop, model = "FIXED"),
    ENV  = list(X = Ze,   model = "BRR"),
    LINE = list(X = Zg,   model = "BRR")
  )
  if (model_name %in% c("M2", "M3")) ETA$GENO <- list(K = K_obs, model = "RKHS")
  if (model_name == "M3")            ETA$GXE  <- list(K = K_gxe, model = "RKHS")
  save_dir <- tempfile("bglr_cv_")
  dir.create(save_dir)
  fit <- tryCatch(
    BGLR(y = y_train, ETA = ETA, nIter = n_iter, burnIn = burn_in, thin = thin,
         verbose = FALSE, saveAt = file.path(save_dir, "bg_")),
    error = function(e) { message("BGLR failed: ", e$message); NULL }
  )
  unlink(save_dir, recursive = TRUE)
  if (is.null(fit)) return(rep(NA_real_, length(y_train)))
  fit$yHat
}

pa_per_env <- function(yHat, y_full, mask_idx) {
  sapply(ENVS, function(e) {
    idx <- mask_idx & (obs$Environment == e)
    yy <- y_full[idx]; pp <- yHat[idx]
    ok <- !is.na(yy) & !is.na(pp)
    if (sum(ok) < 5) return(NA_real_)
    suppressWarnings(cor(pp[ok], yy[ok]))
  })
}

## -----------------------------------------------------------------------------
## Fold assignment (deterministic given rep_id/trait, independent of
## execution order -- safe to recompute inside any parallel worker).
## -----------------------------------------------------------------------------
cv1_mask <- function(rep_id, trait, fold_k) {
  set.seed(1000 * rep_id + nchar(trait))
  g_perm <- sample(GENOTYPES)
  fold_of_g <- setNames(rep(1:N_FOLDS, length.out = length(g_perm)), g_perm)
  obs$Genotype %in% names(fold_of_g)[fold_of_g == fold_k]
}

cv2_mask <- function(rep_id, trait, fold_k, y_full) {
  set.seed(2000 * rep_id + nchar(trait))
  real_idx <- which(!is.na(y_full))
  perm <- sample(real_idx)
  fold_of_row <- rep(1:N_FOLDS, length.out = length(perm))
  fold_assign <- setNames(fold_of_row, perm)
  rows_k <- as.integer(names(fold_assign)[fold_assign == fold_k])
  mask_idx <- rep(FALSE, length(y_full)); mask_idx[rows_k] <- TRUE
  mask_idx
}

## -----------------------------------------------------------------------------
## Main loop -- flattened to one task per (trait, model, CV, rep, fold),
## fit in parallel (fork-based mclapply).
## -----------------------------------------------------------------------------
tasks <- expand.grid(trait = TRAITS, model = c("M1", "M2", "M3"),
                      CV = c("CV1", "CV2"), rep = 1:N_REPS, fold = 1:N_FOLDS,
                      stringsAsFactors = FALSE)
cat("Total tasks:", nrow(tasks), " | parallel workers:", N_CORES, "\n")

y_full_by_trait <- setNames(lapply(TRAITS, build_y_full), TRAITS)

t0 <- Sys.time()
task_results <- mclapply(seq_len(nrow(tasks)), function(i) {
  tr <- tasks$trait[i]; m <- tasks$model[i]; cv <- tasks$CV[i]
  r <- tasks$rep[i]; k <- tasks$fold[i]
  y_full <- y_full_by_trait[[tr]]

  mask_idx <- if (cv == "CV1") cv1_mask(r, tr, k) else cv2_mask(r, tr, k, y_full)
  y_train <- ifelse(mask_idx, NA_real_, y_full)
  yHat <- fit_one_model(y_train, m)
  pa <- pa_per_env(yHat, y_full, mask_idx)

  data.frame(Environment = ENVS, predictive_ability = pa, rep = r, fold = k,
             trait = tr, model = m, CV = cv)
}, mc.cores = N_CORES)

failed <- sapply(task_results, function(x) inherits(x, "try-error"))
cat("Failed tasks:", sum(failed), "of", length(task_results), "\n")
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

out <- bind_rows(task_results[!failed])
write.csv(out, "../results/predictions/predictive_ability_cv1_cv2.csv", row.names = FALSE)
cat("\nSaved: ../results/predictions/predictive_ability_cv1_cv2.csv (", nrow(out), "rows )\n")
