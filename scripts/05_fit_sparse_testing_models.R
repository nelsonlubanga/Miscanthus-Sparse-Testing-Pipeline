## =============================================================================
## 05 -- Fit sparse-testing prediction models M1/M2/M3 (manuscript Equations
##        4, 6, 7), with PopGroup added as a FIXED effect throughout, per
##        Reviewer 2's mandate and the user's explicit instruction ("Fit the
##        sparse testing model with population structure as fixed effect").
## -----------------------------------------------------------------------------
## M1: Environment(random) + Genotype(random)                         (E+L)
## M2: M1 + Genomic random effect g_i (RKHS on VanRaden G)             (E+L+G)
## M3: M2 + Genomic x Environment interaction (Jarquin et al. 2014
##     Hadamard-product kernel: K_gxe = (Zg G Zg') o (Ze Ze'))         (E+L+G+GxE)
## PopGroup is fitted as an additional FIXED effect in every model (not part
## of the original manuscript equations, added per reviewer requirement).
##
## For each (design D1-D6) x (rep 1-10) x (trait DM/FW/MC) x (model M1-M3):
##   - Calibration set = the design's training genotypes per environment
##     (25/environment, per Table 1), using their BLUEs.
##   - Prediction set = the remaining ~100 genotypes/environment (masked to
##     NA for fitting), predicted from the model and compared to their BLUE.
##   - Predictive ability (per environment) = Pearson cor(predicted, BLUE)
##     over the prediction-set genotypes with an actual observed BLUE.
##   - Reported per rep as the across-environment mean predictive ability
##     (matching the manuscript's "repeated 10 times" procedure), so each
##     design x trait x model combination yields 10 values (one per rep) --
##     the natural input to a boxplot.
##
## Output: results/predictions/predictive_ability.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(BGLR)
  library(parallel)
})

TEST_MODE <- FALSE  # set TRUE to run only D1/rep1/DM as a fast sanity check
N_CORES <- 4        # fork-based parallelism (mclapply); moderate to avoid
                     # worsening memory pressure from other running apps

## -----------------------------------------------------------------------------
## Load inputs
## -----------------------------------------------------------------------------
blues <- read.csv("../results/BLUEs/BLUEs_all_environments.csv", stringsAsFactors = FALSE)
load("../data/markers/Gmatrix_125.rda")        # G, M_centered
popgroups <- read.csv("../data/popgroups_125.csv", stringsAsFactors = FALSE)
all_designs <- readRDS("../data/sparse_testing_designs.rds")

GENOTYPES <- rownames(G)                        # canonical 125-genotype order
stopifnot(length(GENOTYPES) == 125)
ENVS <- c("Aber-2015", "Aber-2016", "Brau-2014", "Brau-2015", "Brau-2016")
TRAITS <- c("DM", "FW", "MC")

pop_lookup <- setNames(popgroups$PopGroup, popgroups$Genotype)
stopifnot(all(GENOTYPES %in% names(pop_lookup)))
PopGroup <- factor(pop_lookup[GENOTYPES])

## -----------------------------------------------------------------------------
## Observation-level scaffold: 125 genotypes x 5 environments = 625 rows,
## in a fixed order (genotype-major) reused for every fit.
## -----------------------------------------------------------------------------
obs <- expand.grid(Genotype = GENOTYPES, Environment = ENVS, stringsAsFactors = FALSE)
obs$PopGroup <- pop_lookup[obs$Genotype]

## Design matrices shared across all fits (built once)
Zg <- model.matrix(~ factor(obs$Genotype, levels = GENOTYPES) - 1)   # 625 x 125
Ze <- model.matrix(~ factor(obs$Environment, levels = ENVS) - 1)     # 625 x 5
Xpop_full <- model.matrix(~ factor(obs$PopGroup))                    # intercept + (k-1) dummies
Xpop <- Xpop_full[, -1, drop = FALSE]                                 # drop intercept (BGLR supplies mu)

K_obs <- Zg %*% G %*% t(Zg)                          # 625 x 625, genomic covariance at obs level
K_gxe <- K_obs * (Ze %*% t(Ze))                       # Hadamard product with same-environment indicator

## -----------------------------------------------------------------------------
## Wide BLUE lookup per trait: Genotype x Environment matrix of observed BLUEs
## (NA where no phenotype exists at all, independent of any sparse-testing design)
## -----------------------------------------------------------------------------
build_y_full <- function(trait) {
  sub <- blues |> filter(trait == !!trait) |> select(Genotype, Environment, BLUE)
  wide <- obs |> select(Genotype, Environment) |>
    left_join(sub, by = c("Genotype", "Environment"))
  wide$BLUE
}

## -----------------------------------------------------------------------------
## Fit one model (M1/M2/M3) for one masked y vector, return yHat for all obs.
## -----------------------------------------------------------------------------
fit_one_model <- function(y_train, model_name, n_iter = 1500, burn_in = 300, thin = 3) {
  ETA <- list(
    POP  = list(X = Xpop, model = "FIXED"),
    ENV  = list(X = Ze,   model = "BRR"),
    LINE = list(X = Zg,   model = "BRR")
  )
  if (model_name %in% c("M2", "M3")) {
    ETA$GENO <- list(K = K_obs, model = "RKHS")
  }
  if (model_name == "M3") {
    ETA$GXE <- list(K = K_gxe, model = "RKHS")
  }
  save_dir <- tempfile("bglr_")
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

## -----------------------------------------------------------------------------
## Predictive ability for one (design, rep, trait, model): mask training set,
## fit, evaluate per environment on the prediction set, average across envs.
## -----------------------------------------------------------------------------
run_one <- function(design_rep_key, trait, model_name) {
  dr <- all_designs[[design_rep_key]]
  y_full <- build_y_full(trait)

  is_training <- rep(FALSE, nrow(obs))
  for (e in ENVS) {
    is_training <- is_training | (obs$Environment == e & obs$Genotype %in% dr$training[[e]])
  }
  y_train <- ifelse(is_training, y_full, NA_real_)

  yHat <- fit_one_model(y_train, model_name)

  pa_by_env <- sapply(ENVS, function(e) {
    idx <- obs$Environment == e & !is_training
    yy <- y_full[idx]; pp <- yHat[idx]
    ok <- !is.na(yy) & !is.na(pp)
    if (sum(ok) < 5) return(NA_real_)
    suppressWarnings(cor(pp[ok], yy[ok]))
  })
  mean(pa_by_env, na.rm = TRUE)
}

## -----------------------------------------------------------------------------
## Main loop -- flattened into one task list, fit in parallel (fork-based
## mclapply; shares the read-only G/K matrices via copy-on-write, so memory
## overhead per worker is modest).
## -----------------------------------------------------------------------------
design_rep_keys <- names(all_designs)
if (TEST_MODE) design_rep_keys <- design_rep_keys[grepl("^D1_rep1$", design_rep_keys)]

tasks <- expand.grid(key = design_rep_keys, trait = TRAITS, model = c("M1", "M2", "M3"),
                      stringsAsFactors = FALSE)
cat("Total tasks:", nrow(tasks), " | parallel workers:", N_CORES, "\n")

t0 <- Sys.time()
task_results <- mclapply(seq_len(nrow(tasks)), function(i) {
  key <- tasks$key[i]; tr <- tasks$trait[i]; m <- tasks$model[i]
  dr <- all_designs[[key]]
  pa <- run_one(key, tr, m)
  data.frame(design = dr$design, rep = dr$rep, n_NOL = dr$n_NOL, n_OL = dr$n_OL,
             trait = tr, model = m, predictive_ability = pa)
}, mc.cores = N_CORES)

## Progress/failure summary (mclapply doesn't stream progress like a serial
## loop; report after the fact instead)
failed <- sapply(task_results, function(x) inherits(x, "try-error"))
cat("Failed tasks:", sum(failed), "of", length(task_results), "\n")
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

out <- bind_rows(task_results[!failed])
write.csv(out, "../results/predictions/predictive_ability.csv", row.names = FALSE)
cat("\nSaved: ../results/predictions/predictive_ability.csv (", nrow(out), "rows )\n")
