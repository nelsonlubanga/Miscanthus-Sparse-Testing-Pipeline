## =============================================================================
## 11 -- Line-count-matched designs: isolating NOL/OL from unique-line-count
## -----------------------------------------------------------------------------
## Directly addresses Reviewer 2, Major Comment 3: the D1-D6 sweep (script 05)
## fixes total plots (125) and lets unique training-line count vary (125 in
## D1 down to 25 in D6), confounding NOL/OL proportion with line count. This
## script instead fixes UNIQUE LINE COUNT and lets total plots vary, using
## fewer than all 5 environments at a time (per the reviewer's suggestion),
## averaged over all combinations of K_ENVS environments and multiple reps.
##
## For a given combination of K_ENVS environments: s shared (OL) lines are
## trained in ALL K_ENVS environments; n lines are trained ONLY in each single
## environment (n distinct NOL lines per environment, not shared with any
## other environment or with the shared set). Total unique lines =
## s + K_ENVS*n, held fixed across all 5 designs for a given K_ENVS; total
## plots = K_ENVS*(s+n), which varies as s increases from 0 to TOTAL_UNIQUE.
##
## Design tables (5 evenly-spaced designs per K_ENVS, chosen so unique lines
## is held exactly constant and n is always a whole number):
##   K_ENVS=2: total unique=40, n in {20,15,10,5,0}, s = 40-2n
##   K_ENVS=3: total unique=36, n in {12,9,6,3,0},  s = 36-3n
##   K_ENVS=4: total unique=32, n in {8,6,4,2,0},   s = 32-4n
##
## Prediction set = the remaining (125 - unique lines) genotypes not in any
## of the K_ENVS environments' training sets, evaluated within each
## environment and averaged, exactly as in script 05.
##
## Output: results/predictions/linecount_matched_predictive_ability_k<K_ENVS>.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(BGLR)
  library(parallel)
})

K_ENVS <- 2   # <-- set to 2, 3, or 4 and re-run for each
N_REPS <- 5
N_CORES <- 4
TRAITS <- c("DM", "FW", "MC")
ALL_ENVS <- c("Aber-2015", "Aber-2016", "Brau-2014", "Brau-2015", "Brau-2016")
ENV_COMBOS <- combn(ALL_ENVS, K_ENVS, simplify = FALSE)

design_specs <- list(
  `2` = list(total_unique = 40, n = c(20, 15, 10, 5, 0)),
  `3` = list(total_unique = 36, n = c(12, 9, 6, 3, 0)),
  `4` = list(total_unique = 32, n = c(8, 6, 4, 2, 0))
)
spec <- design_specs[[as.character(K_ENVS)]]
designs <- data.frame(
  design = c("A", "B", "C", "D", "E"),
  n = spec$n,
  s = spec$total_unique - K_ENVS * spec$n
)
stopifnot(all(designs$s + K_ENVS * designs$n == spec$total_unique))
stopifnot(all(designs$s >= 0), all(designs$n >= 0))
cat("K_ENVS =", K_ENVS, " | number of environment combinations:", length(ENV_COMBOS), "\n")
print(designs)

blues <- read.csv("../results/BLUEs/BLUEs_all_environments.csv", stringsAsFactors = FALSE)
load("../data/markers/Gmatrix_125.rda")        # G, M_centered
popgroups <- read.csv("../data/popgroups_125.csv", stringsAsFactors = FALSE)

GENOTYPES <- rownames(G)
stopifnot(length(GENOTYPES) == 125)
pop_lookup <- setNames(popgroups$PopGroup, popgroups$Genotype)

build_y_full_combo <- function(trait, envs) {
  obs <- expand.grid(Genotype = GENOTYPES, Environment = envs, stringsAsFactors = FALSE)
  sub <- blues |> filter(trait == !!trait, Environment %in% envs) |>
    select(Genotype, Environment, BLUE)
  wide <- obs |> left_join(sub, by = c("Genotype", "Environment"))
  list(obs = obs, y_full = wide$BLUE)
}

fit_one_model_combo <- function(y_train, model_name, Xpop, Ze, Zg, K_obs, K_gxe,
                                 n_iter = 1500, burn_in = 300, thin = 3) {
  ETA <- list(
    POP  = list(X = Xpop, model = "FIXED"),
    ENV  = list(X = Ze,   model = "BRR"),
    LINE = list(X = Zg,   model = "BRR")
  )
  if (model_name %in% c("M2", "M3")) ETA$GENO <- list(K = K_obs, model = "RKHS")
  if (model_name == "M3")            ETA$GXE  <- list(K = K_gxe, model = "RKHS")
  save_dir <- tempfile("bglr_lc_")
  dir.create(save_dir)
  fit <- tryCatch(
    BGLR(y = y_train, ETA = ETA, nIter = n_iter, burnIn = burn_in, thin = thin,
         verbose = FALSE, saveAt = file.path(save_dir, "bg_")),
    error = function(e) NULL
  )
  unlink(save_dir, recursive = TRUE)
  if (is.null(fit)) return(rep(NA_real_, length(y_train)))
  fit$yHat
}

## Build one design's training-genotype assignment for a given combination of
## K_ENVS environments and replicate seed: s shared genotypes (trained in
## ALL K_ENVS environments), n genotypes unique to EACH environment (K_ENVS
## distinct sets of n genotypes, none shared with each other or with the
## shared set).
build_training <- function(envs, s, n, seed) {
  set.seed(seed)
  shared <- if (s > 0) sample(GENOTYPES, s) else character(0)
  remaining <- setdiff(GENOTYPES, shared)
  k <- length(envs)
  nol_pool <- if (n > 0) sample(remaining, k * n) else character(0)
  training <- list()
  for (i in seq_along(envs)) {
    nol_i <- if (n > 0) nol_pool[((i - 1) * n + 1):(i * n)] else character(0)
    training[[envs[i]]] <- c(shared, nol_i)
  }
  training
}

tasks <- expand.grid(combo_idx = seq_along(ENV_COMBOS), design = designs$design,
                      rep = seq_len(N_REPS), trait = TRAITS, model = c("M1", "M2", "M3"),
                      stringsAsFactors = FALSE)
cat("Total tasks:", nrow(tasks), " | parallel workers:", N_CORES, "\n")

t0 <- Sys.time()
task_results <- mclapply(seq_len(nrow(tasks)), function(i) {
  envs <- ENV_COMBOS[[tasks$combo_idx[i]]]
  des <- designs[designs$design == tasks$design[i], ]
  tr <- tasks$trait[i]; m <- tasks$model[i]; r <- tasks$rep[i]

  built <- build_y_full_combo(tr, envs)
  obs <- built$obs; y_full <- built$y_full
  obs$PopGroup <- pop_lookup[obs$Genotype]

  seed <- 5000 * tasks$combo_idx[i] + 100 * r + as.integer(factor(tasks$design[i], levels = designs$design))
  trn <- build_training(envs, des$s, des$n, seed)

  is_training <- rep(FALSE, nrow(obs))
  for (e in envs) {
    is_training <- is_training | (obs$Environment == e & obs$Genotype %in% trn[[e]])
  }
  y_train <- ifelse(is_training, y_full, NA_real_)

  Zg <- model.matrix(~ factor(obs$Genotype, levels = GENOTYPES) - 1)
  Ze <- model.matrix(~ factor(obs$Environment, levels = envs) - 1)
  Xpop_full <- model.matrix(~ factor(obs$PopGroup))
  Xpop <- Xpop_full[, -1, drop = FALSE]
  K_obs <- Zg %*% G %*% t(Zg)
  K_gxe <- K_obs * (Ze %*% t(Ze))

  yHat <- fit_one_model_combo(y_train, m, Xpop, Ze, Zg, K_obs, K_gxe)

  pa_by_env <- sapply(envs, function(e) {
    idx <- obs$Environment == e & !is_training
    yy <- y_full[idx]; pp <- yHat[idx]
    ok <- !is.na(yy) & !is.na(pp)
    if (sum(ok) < 5) return(NA_real_)
    suppressWarnings(cor(pp[ok], yy[ok]))
  })
  pa <- mean(pa_by_env, na.rm = TRUE)

  data.frame(combo = paste(envs, collapse = "|"), k_envs = length(envs),
             design = tasks$design[i], s = des$s, n = des$n,
             unique_lines = des$s + length(envs) * des$n,
             total_plots = length(envs) * (des$s + des$n),
             rep = r, trait = tr, model = m, predictive_ability = pa)
}, mc.cores = N_CORES)

failed <- sapply(task_results, function(x) inherits(x, "try-error"))
cat("Failed tasks:", sum(failed), "of", length(task_results), "\n")
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

out <- bind_rows(task_results[!failed])
out_file <- sprintf("../results/predictions/linecount_matched_predictive_ability_k%d.csv", K_ENVS)
write.csv(out, out_file, row.names = FALSE)
cat("\nSaved:", out_file, "(", nrow(out), "rows )\n")

cat("\n=== Mean predictive ability by design (averaged over all env pairs, reps, traits, models) ===\n")
agg <- out |> group_by(design, s, n, unique_lines, total_plots) |>
  summarise(mean_PA = mean(predictive_ability, na.rm = TRUE), .groups = "drop") |>
  arrange(s)
print(as.data.frame(agg))

cat("\n=== By trait ===\n")
agg2 <- out |> group_by(design, trait) |>
  summarise(mean_PA = mean(predictive_ability, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = trait, values_from = mean_PA)
print(as.data.frame(agg2))

cat("\nCorrelation between overlap level (s) and predictive ability, at FIXED unique-line count:\n")
cat(round(cor(agg$s, agg$mean_PA), 4), "\n")
