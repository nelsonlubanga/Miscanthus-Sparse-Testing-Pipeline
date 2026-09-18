## =============================================================================
## 04 -- Build the VanRaden genomic relationship matrix (G) for the 125
##        genotypes used in the sparse-testing prediction models.
## -----------------------------------------------------------------------------
## G ~ XX'/p (manuscript Equation 6 context), computed on the centred,
## imputed marker matrix -- the same imputed matrix already validated in
## 00a_derive_population_groups.R (rrBLUP::A.mat), so PopGroup and G are
## built from one consistent marker source.
##
## Output: data/markers/Gmatrix_125.rda  (objects: G, M_centered)
##   G: 125 x 125 genomic relationship matrix, rownames/colnames = Genotype
##   M_centered: 125 x p centred (not standardized) marker matrix, used for
##     the M2 marker-regression (BRR) formulation of g_i, equivalent to G.
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages(library(rrBLUP))

load("../data/markers/SNPs_125genotypes.rda")   # loads `SNPs`
stopifnot(nrow(SNPs) == 125)

snp.imp <- A.mat(SNPs, return.imputed = TRUE)
G <- snp.imp$A
M_centered <- scale(snp.imp$imputed, center = TRUE, scale = FALSE)

rownames(G) <- colnames(G) <- rownames(SNPs)
rownames(M_centered) <- rownames(SNPs)

## Ensure positive-definiteness for RKHS use in BGLR (small ridge if needed)
eig <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
if (min(eig) <= 1e-8) {
  cat("G is not strictly PD (min eigenvalue =", min(eig), "); adding small ridge.\n")
  G <- G + diag(1e-6, nrow(G))
}

save(G, M_centered, file = "../data/markers/Gmatrix_125.rda")
cat("Saved ../data/markers/Gmatrix_125.rda\n")
cat("dim(G):", dim(G), "\n")
cat("range(diag(G)):", range(diag(G)), "\n")
