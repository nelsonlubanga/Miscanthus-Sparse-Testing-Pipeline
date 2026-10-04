## Heritability with population group fitted as random (per environment, ASReml):
## H2 = (Vpop + Vg)/(Vpop + Vg + Ve/r_bar) includes group differences;
## H2_within = Vg/(Vg + Ve/r_bar) excludes them and matches script 02 (group fixed).
## Run: Rscript scripts/02b_estimate_heritability_popgroup_random.R

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_file <- if (!is.null(source_file)) source_file else if (length(script_arg)) gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE) else NULL
root <- if (!is.null(script_file)) dirname(dirname(normalizePath(script_file))) else if (basename(getwd()) == "scripts") dirname(getwd()) else getwd()
outdir <- file.path(root, "results", "heritability", "popgroup_random")
suppressPackageStartupMessages(library(dplyr))

run_h2 <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) stop("ASReml-R required")
  suppressPackageStartupMessages(library(asreml))
d <- read.csv(file.path(root, "data/combined_phenotypes_long.csv"), stringsAsFactors = FALSE)
mapping <- unique(d[c("Genotype", "PopGroup")])
if (anyNA(mapping) || anyDuplicated(mapping$Genotype)) stop("Each genotype must map to exactly one non-missing PopGroup")
d$Genotype <- factor(d$Genotype)
d$PopGroup <- factor(d$PopGroup)
d$Block <- factor(d$Block)
d$Row <- factor(d$Row)
d$Col <- factor(d$Col)

traits <- unique(d$trait)
envs <- unique(d$Environment)

fit_h2_one <- function(sub, has_row_col) {
  needed <- c("value", "Genotype", "PopGroup", "Block",
              if (has_row_col) c("Row", "Col"))
  sub <- sub[complete.cases(sub[, needed]), , drop = FALSE]
  sub <- sub[is.finite(sub$value), , drop = FALSE]
  sub$Genotype <- droplevels(sub$Genotype)
  sub$PopGroup <- droplevels(sub$PopGroup)
  sub$Block <- droplevels(sub$Block)
  sub$Row <- droplevels(sub$Row)
  sub$Col <- droplevels(sub$Col)
  if (n_distinct(sub$Genotype) < 5 || nlevels(sub$PopGroup) < 2) return(NULL)

  if (all(table(sub$Genotype) == 1L)) {
    warning("No within-genotype replication: genotype and residual variance cannot be separated")
    return(NULL)
  }

  ## Aberystwyth: Block+Row+Col (Block independently recorded, not
  ## confounded with Row/Col). Braunschweig: Block only (its Block is
  ## reconstructed from Row, so adding Row/Col would duplicate the same
  ## spatial information -- see 01_estimate_BLUEs.R and
  ## 00_prepare_combined_phenotypes.R for the reconstruction).
  random_form <- if (has_row_col) ~ PopGroup + Genotype + Block + Row + Col else ~ PopGroup + Genotype + Block
  fit <- tryCatch(
    asreml(fixed = value ~ 1, random = random_form, data = sub,
           trace = FALSE, maxit = 30),
    error = function(e) { warning(conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(NULL)
  for (i in 1:8) {
    if (isTRUE(fit$converge)) break
    fit <- tryCatch(update.asreml(fit), error = function(e) fit)
  }

  if (!isTRUE(fit$converge)) {
    warning("Non-converged ASReml fit; no heritability reported")
    return(NULL)
  }
  vc <- summary(fit)$varcomp
  components <- data.frame(term = rownames(vc), vc, row.names = NULL)
  get_v <- function(term) {
    if (!term %in% rownames(vc)) stop("Missing variance component: ", term)
    as.numeric(vc[term, "component"])
  }
  sigma2_PopGroup <- get_v("PopGroup")
  if (!("Genotype" %in% rownames(vc))) return(NULL)
  sigma2_G <- get_v("Genotype")
  sigma2_e <- get_v("units!R")
  r_bar <- mean(table(sub$Genotype))
  ## Same approximate entry-mean convention as script 02. Spatial terms
  ## adjust the fit but are not added to the entry-mean denominator.
  ## H2 assumes between-group variance represents genetic differences;
  ## this does not establish their causal origin or additive inheritance.
  H2_within <- sigma2_G / (sigma2_G + sigma2_e / r_bar)
  denom <- sigma2_PopGroup + sigma2_G + sigma2_e / r_bar
  H2 <- (sigma2_PopGroup + sigma2_G) / denom
  genotype_fraction_total <- sigma2_G / denom

  result <- data.frame(sigma2_PopGroup = sigma2_PopGroup, H2_within = H2_within,
             genotype_fraction_total = genotype_fraction_total, sigma2_G = sigma2_G, sigma2_e = sigma2_e, r_bar = r_bar,
             H2 = H2, n_genotypes = n_distinct(sub$Genotype),
             converged = isTRUE(fit$converge))
  attr(result, "components") <- components
  result
}

results <- list()
vc_results <- list()
failed <- list()
for (env in envs) {
  has_row_col <- grepl("Aber-", env, fixed = TRUE)
  for (tr in traits) {
    sub <- d |> filter(Environment == env, trait == tr)
    if (nrow(sub) == 0) next
    out <- fit_h2_one(sub, has_row_col)
    if (is.null(out)) {
      failed[[paste(env, tr)]] <- data.frame(Environment = env, trait = tr,
          status = "Insufficient data, fit error or non-convergence; see warnings")
      next
    }
    comps <- attr(out, "components")
    comps$Environment <- env; comps$trait <- tr
    vc_results[[paste(env, tr)]] <- comps
    out$Environment <- env
    out$trait <- tr
    results[[paste(env, tr)]] <- out
    cat(env, tr, ": H2=", round(out$H2, 3), " r_bar=", round(out$r_bar, 2),
        " converged=", out$converged, "\n")
  }
}

if (length(failed)) write.csv(bind_rows(failed), file.path(outdir, "heritability_failures.csv"), row.names = FALSE)
if (!length(results)) stop("No successful heritability fits. Check ASReml licence and warnings.")
write.csv(bind_rows(vc_results), file.path(outdir, "heritability_variance_components.csv"), row.names = FALSE)
h2_tab <- bind_rows(results) |>
  select(Environment, trait, n_genotypes, r_bar, sigma2_PopGroup, sigma2_G, sigma2_e, H2, H2_within, genotype_fraction_total, converged)
write.csv(h2_tab, file.path(outdir, "heritability_by_environment.csv"), row.names = FALSE)
cat("\nSaved random-PopGroup heritability outputs in", outdir, "\n")
print(h2_tab, digits = 3)

}

run_h2()
