# Reviewer DM analysis

This separate pipeline addresses group differences and repeated years within a trial. Original analyses remain available. Start with RESULTS.md for findings and limitations.

## Phenotypes and group removal

The phenotype is one pooled DM BLUE per genotype per physical trial, estimated previously with ASReml using fixed Year + Genotype, random Genotype:Year + persistent Plot + Year:Block, with Year:Row and Year:Col additionally at ABR33, and year-specific residual variance. Genotype predictions equally average fixed years. See ../scripts/site_pooled_blue_revision.R and supporting_results/ for model checks and estimates. ABR33 has 2015/2016; JKI has 2014/2015/2016. The JKI-2014 residual variance is on its lower boundary. These are pooled BLUEs, not genotype BLUPs.

For each validation split, regress the training BLUEs on Trial × PopGroup (equivalently subtract each training trial/group mean), then apply those same offsets to validation BLUEs. Training residual means are zero; validation group means need not be zero. Full-data zero-mean residuals are saved only as descriptive checks, never as model inputs. This implements the reviewer's expressly suggested regression alternative. We do not assume random effects have realized group means exactly zero, or that a full fitted BLUP including group effects removes species differences.

## Prediction models

Both raw and group-adjusted branches omit a POP effect from BGLR. M1 contains intercept, random trial and independent random line; M2 adds a genomic effect with the supplied G matrix; M3 additionally adds genomic-by-trial covariance G multiplied elementwise by the trial incidence covariance. The genomic matrix is not recomputed within groups. Genomic structure may therefore still be present; within-group evaluation is essential.

Each branch uses identical splits and model-specific RNG seeds. Explicit priors use only raw training variance and are numerically identical between response branches. Variance prior df=5; each random term's prior mode corresponds to 0.1 times raw training variance, adjusted by mean kernel diagonal; residual prior mode is 0.5 times raw training variance. Main chains use 8,000 iterations, 2,000 burn-in, thinning 5. The source records exact settings and package versions.

CV1: five-fold genotype-wise masking at BOTH trials, ten repetitions. CV2: mask genotype/trial records while retaining the genotype at the OTHER physical trial; the genotype observed at only one trial is training-only. Sparse25: exactly 25 observed genotypes per trial, 0/5/10/15/20/25 shared training genotypes, ten repetitions. Unique training genotypes equal 50 minus overlap. All yearly records are pooled before prediction, so a masked trial phenotype has no separate same-trial year phenotype in the prediction training set.

There are 160 splits × 3 models × 2 response branches = 960 primary fits. M1 CV1 has constant predictions within a trial and cannot rank unseen genotypes: its PA is undefined, not zero or a Monte Carlo-generated correlation. Its RMSE is still meaningful. Independent line solutions for unseen genotypes are set to their exact posterior mean zero.

## Metrics and figures

fold_scores.csv contains fold-level target-scale PA and RMSE. CV1_paired_targets and CV2_paired_targets show those correlations. Raw and residual PA refer to different response targets and must not be interpreted as an exact variance decomposition.

per_repeat_metrics.csv pools out-of-fold predictions within each repetition. PA_within_pooled centers raw observations and reconstructed raw predictions by group FOR SCORING ONLY. This evaluates ranking beyond group means; its values differ from target-scale PA, which can retain validation group offsets. Reconstructed raw predictions add each training-derived group offset back to residual predictions. within_group_summary.csv and species_specific_prediction.png show each group separately. CV1_within_group and CV2_within_group summarize repeated pooled within-group correlations. Group-only RMSE provides a training-fitted baseline on the same raw scale.

Ten repeated splits reuse the same genotypes: boxplots describe partition sensitivity, not ten independent biological experiments. No naive significance tests across repetitions are used.

## Reproduction and validation

Run scripts in numerical order with Rscript using the original project directory layout. Script 01 refuses to reuse checkpoints if its input/signature changes. Script 02 independently validates masks, training-only means, paired targets, sparse budgets and all saved fold correlations/RMSE. Scripts 03 and 04 run representative independent chains; they do not certify all 960 fits. input_provenance.csv records first-stage sources and hashes. results/sessionInfo.txt records the main R session.

First-stage BLUEs were estimated from the full plot dataset before CV. This is second-stage adjusted-phenotype validation, not fully nested plot-level validation. BLUE uncertainty is not propagated or weighted in BGLR. Results concern two observed trials and do not establish prediction at arbitrary new sites. See RESULTS.md for numerical chain-check limitations.
