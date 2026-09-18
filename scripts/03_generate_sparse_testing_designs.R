## =============================================================================
## 03 -- Generate sparse-testing designs D1-D6 (NOL/OL) and CV1/CV2 folds
## -----------------------------------------------------------------------------
## Reproduces the manuscript's Table 1 design exactly:
##   Each of the 5 environments always has 25 "training" genotypes and 100
##   held-out ("prediction set") genotypes. What varies across D1-D6 is how
##   many of the 25 training genotypes per environment are NOL (unique to
##   that one environment) vs OL (the SAME shared set of genotypes, present
##   in the training set of every environment):
##     D1: 25 NOL + 0  OL  (125 unique genotypes overall)
##     D2: 20 NOL + 5  OL  (105 unique genotypes overall)
##     D3: 15 NOL + 10 OL  ( 85 unique genotypes overall)
##     D4: 10 NOL + 15 OL  ( 65 unique genotypes overall)
##     D5:  5 NOL + 20 OL  ( 45 unique genotypes overall)
##     D6:  0 NOL + 25 OL  ( 25 unique genotypes overall)
## Repeated for N_REPS independent random assignments (manuscript uses 10).
##
## CV1: within a design's training genotypes, 5-fold CV masking whole
##      genotypes (all 5 environments) -- predicting unphenotyped new material.
## CV2: within a design's training genotypes, 5-fold CV masking individual
##      genotype x environment cells -- sparse testing proper.
##
## Output: data/sparse_testing_designs.rds -- a list, one element per
##   (design, rep), each containing the training-set genotype assignment
##   per environment and the CV1/CV2 fold assignments.
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages(library(dplyr))

N_REPS <- 10
N_FOLDS <- 5
ENVS <- c("Aber-2015", "Aber-2016", "Brau-2014", "Brau-2015", "Brau-2016")
N_TRAIN_PER_ENV <- 25

designs <- data.frame(
  design = paste0("D", 1:6),
  n_NOL = c(25, 20, 15, 10, 5, 0),
  n_OL  = c(0, 5, 10, 15, 20, 25)
)
stopifnot(all(designs$n_NOL + designs$n_OL == N_TRAIN_PER_ENV))

popgroups <- read.csv("../data/popgroups_125.csv", stringsAsFactors = FALSE)
genotypes_125 <- popgroups$Genotype
stopifnot(length(genotypes_125) == 125)

set.seed(2024)

build_one_design_rep <- function(n_NOL, n_OL, genotypes) {
  ## OL genotypes: one shared set of n_OL genotypes, present in the
  ## training set of every environment.
  ol_genotypes <- if (n_OL > 0) sample(genotypes, n_OL) else character(0)
  remaining <- setdiff(genotypes, ol_genotypes)

  ## NOL genotypes: n_NOL * 5 DISTINCT genotypes needed (unique per
  ## environment); manuscript requires enough genotypes to cover this
  ## without reuse across environments for the NOL portion.
  need <- n_NOL * length(ENVS)
  stopifnot(need <= length(remaining))
  nol_pool <- sample(remaining, need)

  training <- list()
  for (i in seq_along(ENVS)) {
    nol_i <- if (n_NOL > 0) nol_pool[((i - 1) * n_NOL + 1):(i * n_NOL)] else character(0)
    training[[ENVS[i]]] <- c(nol_i, ol_genotypes)
  }
  list(training = training, ol_genotypes = ol_genotypes)
}

assign_cv1_folds <- function(training_genotypes_all) {
  ## training_genotypes_all: the UNION of genotypes appearing in any
  ## environment's training set for this design/rep. CV1 masks a whole
  ## genotype (all environments) per fold.
  g <- sample(training_genotypes_all)
  fold <- rep(1:N_FOLDS, length.out = length(g))
  setNames(fold, g)
}

assign_cv2_folds <- function(training_list) {
  ## CV2 masks individual genotype x environment CELLS. Build the long
  ## list of (genotype, environment) training pairs and fold-assign each
  ## cell independently.
  cells <- bind_rows(lapply(names(training_list), function(e) {
    data.frame(Environment = e, Genotype = training_list[[e]], stringsAsFactors = FALSE)
  }))
  cells <- cells[sample(nrow(cells)), ]
  cells$fold <- rep(1:N_FOLDS, length.out = nrow(cells))
  cells
}

all_designs <- list()
for (d in seq_len(nrow(designs))) {
  for (r in seq_len(N_REPS)) {
    key <- paste0(designs$design[d], "_rep", r)
    built <- build_one_design_rep(designs$n_NOL[d], designs$n_OL[d], genotypes_125)
    union_genotypes <- unique(unlist(built$training))
    cv1 <- assign_cv1_folds(union_genotypes)
    cv2 <- assign_cv2_folds(built$training)
    all_designs[[key]] <- list(
      design = designs$design[d], rep = r,
      n_NOL = designs$n_NOL[d], n_OL = designs$n_OL[d],
      training = built$training, ol_genotypes = built$ol_genotypes,
      n_unique_genotypes = length(union_genotypes),
      cv1_folds = cv1, cv2_cells = cv2
    )
  }
}

saveRDS(all_designs, "../data/sparse_testing_designs.rds")

## --- Validate against manuscript Table 1 exactly ---
cat("=== Validation against manuscript Table 1 ===\n")
val <- bind_rows(lapply(all_designs[grepl("_rep1$", names(all_designs))], function(x) {
  data.frame(design = x$design, n_NOL = x$n_NOL, n_OL = x$n_OL,
             total_plots = sum(sapply(x$training, length)),
             n_unique_genotypes = x$n_unique_genotypes)
}))
print(val)
cat("\nExpected n_unique_genotypes by design: D1=125, D2=105, D3=85, D4=65, D5=45, D6=25\n")
cat("Expected total_plots: 125 for every design\n")
stopifnot(all(val$total_plots == 125))
stopifnot(all(val$n_unique_genotypes == c(125, 105, 85, 65, 45, 25)))
cat("PASSED: matches Table 1 exactly.\n")
