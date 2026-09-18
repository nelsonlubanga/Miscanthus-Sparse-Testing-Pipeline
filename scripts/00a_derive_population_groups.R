## =============================================================================
## 00a -- Derive population-structure group membership (DAPC clustering)
## -----------------------------------------------------------------------------
## Directly addresses Reviewer 2, Major comment 1: population structure must
## be an explicit, reproducible part of the analysis, not an inherited/
## untraced column. This reproduces the manuscript's own stated method
## ("genetic clusters inferred from the DAPC ... using the adegenet package",
## Jombart, 2008) end to end from the marker data, rather than reusing a
## population-group column of unclear provenance.
##
## Method: find.clusters (adegenet) on the additive-imputed marker matrix
## (rrBLUP::A.mat), forced to 3 clusters (M. sacchariflorus, M. sinensis,
## interspecific hybrids), matching the original clustering exactly
## (same seed, same n.pca).
##
## Validation: an earlier, partial run of this same clustering (91 of the
## 125 genotypes, saved before this pipeline existed) is used as a ground
## truth. This script re-derives the clustering for all 125 genotypes and
## hard-fails (stopifnot) unless it reproduces 100% agreement with those 91
## known labels -- i.e., the pipeline will not silently proceed on a
## different/inconsistent clustering.
##
## Input:  data/markers/SNPs_125genotypes.rda        (125 x 13,687 marker matrix)
##         data/markers/pca_validation_91genotypes.csv (known labels, validation only)
## Output: data/popgroups_125.csv
## =============================================================================

rm(list = ls())
setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")
suppressMessages({
  library(rrBLUP)
  library(adegenet)
})

load("../data/markers/SNPs_125genotypes.rda")   # loads object `SNPs`
stopifnot(nrow(SNPs) == 125)

snp.imp <- A.mat(SNPs, return.imputed = TRUE)
Geno <- snp.imp$imputed

set.seed(165)
grp <- find.clusters(Geno, max.n.clust = 10, n.pca = 50, n.clust = 3, choose.n.clust = FALSE)
d <- data.frame(sample.id = names(grp$grp), grp_new = as.character(grp$grp),
                 stringsAsFactors = FALSE)

## find.clusters assigns arbitrary numeric cluster IDs that are not
## guaranteed to match between runs; this permutation was determined by
## validating against the known 91-genotype labels below and must be
## re-derived if the marker input changes.
map <- c("1" = "3", "2" = "1", "3" = "2")
d$group3 <- map[d$grp_new]
d$PopGroup <- c("1" = "M. sacchariflorus", "2" = "M. sinensis", "3" = "M. x gig")[d$group3]

## --- Validation against the known 91-genotype subset ---
known <- read.csv("../data/markers/pca_validation_91genotypes.csv", stringsAsFactors = FALSE)[, c("sample.id", "group3")]
chk <- merge(d, known, by = "sample.id", suffixes = c(".new", ".known"))
stopifnot(nrow(chk) == 91)
agree <- mean(chk$group3.new == chk$group3.known)
cat("Validation against known 91-genotype labels: agreement =", agree, "\n")
stopifnot(agree == 1)
cat("PASSED: exact agreement with the previously established 91-genotype clustering.\n")

out <- d[, c("sample.id", "PopGroup")]
colnames(out)[1] <- "Genotype"
write.csv(out, "../data/popgroups_125.csv", row.names = FALSE)
cat("\nWrote ../data/popgroups_125.csv (", nrow(out), "genotypes )\n")
print(table(out$PopGroup))
