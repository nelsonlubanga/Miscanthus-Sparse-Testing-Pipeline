## POOLING SENSITIVITY: DM only, one adjusted BLUE per genotype per physical trial.
## Trials are ABR33 and JKI; a held-out trial record represents ALL its years.
## CV2 allows observations only at the other trial. CV1 masks both trials.
## PopGroup remains random as in the previous analysis. This isolates pooling
## years; it does not remove group/species means or resolve that separate concern.
## CV1/CV2 with random PopGroup; separate outputs, original scripts unchanged.
## Run: Rscript scripts/cv1_cv2_popgroup_random.R
## CV1: genotype held out in every environment.
## CV2: observed cells assigned within genotype to distinct folds, ensuring
## another observed environment remains in training. Singletons train only.
## Five folds, ten repeats, identical folds across M1/M2/M3.
## All models: intercept + random PopGroup + environment + line.
## M2 adds genomic; M3 adds genomic-by-environment. Stage-1 BLUEs are reused
## without precision weights; validation is at the BLUE/second-stage level.
## Relationship matrix/group membership use genotypes, not held-out phenotypes.
## Priors: df0=5; each shared term has prior mode 0.1 * training variance on
## mean marginal covariance scale; residual mode 0.5 * training variance.
## Unlike automatic BGLR allocation, shared priors do not change with model.
## 8000 iterations, 2000 burn-in, thin 5. Completion is not convergence proof.
## Fold-level correlations from repeated CV are dependent, not independent CIs.
suppressPackageStartupMessages({library(BGLR);library(dplyr);library(ggplot2)})
arg <- grep('^--file=', commandArgs(), value=TRUE)
sf <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if(is.null(sf) && length(arg)) sf <- gsub('~+~',' ',sub('^--file=','',arg[1]),fixed=TRUE)
root <- if(!is.null(sf)) dirname(dirname(normalizePath(sf))) else if(basename(getwd())=='scripts') dirname(getwd()) else getwd()
outdir <- file.path(root,'results','site_pooled_revision','cv')
dir.create(file.path(outdir,'checkpoints'),recursive=TRUE,showWarnings=FALSE)
N_REPS <- 10L; N_FOLDS <- 5L; N_CORES <- 4L
N_ITER <- 8000L; BURN_IN <- 2000L; THIN <- 5L
TRAITS <- 'DM'; ENVS <- c('ABR33','JKI')
b <- read.csv(file.path(root,'results/site_pooled_revision/pooled_DM_BLUEs.csv'))
pops <- read.csv(file.path(root,'data/popgroups_125.csv'))
ge <- new.env();load(file.path(root,'data/markers/Gmatrix_125.rda'),envir=ge);G <- ge$G
gids <- rownames(G)
stopifnot(identical(gids,colnames(G)),!anyDuplicated(gids),!anyDuplicated(pops$Genotype),
          !anyDuplicated(b[c('Genotype','Environment','trait')]),all(is.finite(G)),max(abs(G-t(G)))<1e-8,
          all(b$converged),all(gids %in% pops$Genotype))
obs <- expand.grid(Genotype=gids,Environment=ENVS,stringsAsFactors=FALSE)
obs$PopGroup <- pops$PopGroup[match(obs$Genotype,pops$Genotype)]
stopifnot(!anyNA(obs$PopGroup))
Zp <- model.matrix(~factor(obs$PopGroup)-1)
Zl <- model.matrix(~factor(obs$Genotype,levels=gids)-1)
Ze <- model.matrix(~factor(obs$Environment,levels=ENVS)-1)
Kg <- Zl %*% G %*% t(Zl); Kge <- Kg * tcrossprod(Ze)
## Precompute eigen decompositions once, preserving exactly the same kernels.
eigen_kernel <- function(k) {
 e <- eigen(k,symmetric=TRUE)
 stopifnot(min(e$values)> -1e-7)
 keep <- e$values>1e-10
 list(V=e$vectors[,keep,drop=FALSE],d=e$values[keep])
}
eg <- eigen_kernel(Kg); ege <- eigen_kernel(Kge)
ys <- setNames(lapply(TRAITS,function(tr) {
 z <- b[b$trait==tr,]; z$BLUE[match(paste(obs$Genotype,obs$Environment),paste(z$Genotype,z$Environment))]
}),TRAITS)
stopifnot(all(vapply(ys,length,1L)==nrow(obs)))
folds <- list(); fold_audit <- list()
for(tr in TRAITS) for(cv in c('CV1','CV2')) for(r in seq_len(N_REPS)) {
 set.seed(100000L*match(cv,c('CV1','CV2'))+1000L*r+match(tr,TRAITS))
 y <- ys[[tr]]; f <- integer(length(y))
 if(cv=='CV1') {
  fg <- sample(rep(seq_len(N_FOLDS),length.out=length(gids)))
  f <- fg[match(obs$Genotype,gids)]
 } else {
  for(g in gids) {
   ii <- which(obs$Genotype==g & is.finite(y))
   if(length(ii)>=2) f[ii] <- sample(seq_len(N_FOLDS),length(ii),replace=FALSE)
  }
 }
 for(k in seq_len(N_FOLDS)) {
  test <- f==k & is.finite(y); train <- f!=k & is.finite(y)
  overlap <- intersect(obs$Genotype[test],obs$Genotype[train])
  if(cv=='CV1') stopifnot(length(overlap)==0)
  if(cv=='CV2') stopifnot(all(obs$Genotype[test] %in% obs$Genotype[train]))
  stopifnot(!any(test & train))
 }
 key <- paste(tr,cv,r,sep='_'); folds[[key]] <- f
 fold_audit[[key]] <- data.frame(trait=tr,CV=cv,rep=r,obs,fold=f,observed=is.finite(y))
}
write.csv(bind_rows(fold_audit),file.path(outdir,'fold_assignments.csv'),row.names=FALSE)
saveRDS(folds,file.path(outdir,'fold_assignments.rds'))
tasks <- expand.grid(trait=TRAITS,model=c('M1','M2','M3'),CV=c('CV1','CV2'),rep=seq_len(N_REPS),fold=seq_len(N_FOLDS),stringsAsFactors=FALSE)
## Refuse stale checkpoints if inputs, code or settings change.
input_files <- c(file.path(root,'data/markers/Gmatrix_125.rda'),file.path(root,'data/popgroups_125.csv'),file.path(root,'results/site_pooled_revision/pooled_DM_BLUEs.csv'),if(!is.null(sf)) normalizePath(sf))
signature <- list(md5=tools::md5sum(input_files),iterations=c(N_ITER,BURN_IN,THIN),folds=folds)
sigfile <- file.path(outdir,'run_signature.rds')
if(file.exists(sigfile)) stopifnot(identical(readRDS(sigfile),signature)) else saveRDS(signature,sigfile)
fit_task <- function(i) {
 t <- tasks[i,]; key <- paste(t$trait,t$CV,t$rep,sep='_')
 y <- ys[[t$trait]]; masked <- folds[[key]]==t$fold
 yt <- y; yt[masked] <- NA_real_
 stopifnot(all(is.na(yt[masked])))
 v <- var(yt,na.rm=TRUE);stopifnot(is.finite(v),v>0)
 brr <- function(x) list(X=x,model='BRR',df0=5,S0=7*.1*v/mean(rowSums(x*x)))
 rk <- function(e,k) list(V=e$V,d=e$d,model='RKHS',df0=5,S0=7*.1*v/mean(diag(k)))
 ETA <- list(POP=brr(Zp),ENV=brr(Ze),LINE=brr(Zl))
 if(t$model!='M1') ETA$GENO <- rk(eg,Kg)
 if(t$model=='M3') ETA$GXE <- rk(ege,Kge)
 td <- tempfile('cv_random_');dir.create(td);on.exit(unlink(td,recursive=TRUE))
 set.seed(900000L+i)
 fit <- BGLR(y=yt,ETA=ETA,df0=5,S0=7*.5*v,nIter=N_ITER,burnIn=BURN_IN,thin=THIN,
             verbose=FALSE,saveAt=file.path(td,'bg_'))
 stopifnot(all(is.finite(fit$yHat)))
 test <- masked & is.finite(y)
 predictions <- data.frame(t,obs[test,],observed=y[test],predicted=fit$yHat[test],row.names=NULL)
 scores <- bind_rows(lapply(ENVS,function(e) {
  ix <- test & obs$Environment==e
  reason <- if(sum(ix)<5) 'fewer_than_5' else if(sd(y[ix])==0 || sd(fit$yHat[ix])==0) 'zero_variance' else 'ok'
  data.frame(t,Environment=e,n_test=sum(ix),predictive_ability=if(reason=='ok') cor(y[ix],fit$yHat[ix]) else NA_real_,status=reason,row.names=NULL)
 }))
 traces <- lapply(list.files(td,pattern='(varB|varU|varE)\\.dat$',full.names=TRUE),function(f) scan(f,quiet=TRUE)[-(seq_len(BURN_IN/THIN))])
 names(traces) <- basename(list.files(td,pattern='(varB|varU|varE)\\.dat$',full.names=TRUE))
 list(scores=scores,predictions=predictions,traces=traces)
}
cat('Running',nrow(tasks),'fits with random PopGroup\n')
for(start in seq(1,nrow(tasks),by=36)) {
 ids <- start:min(start+35,nrow(tasks))
 rr <- parallel::mclapply(ids,function(i) {
  path <- file.path(outdir,'checkpoints',sprintf('task_%04d.rds',i))
  if(file.exists(path)) return(TRUE)
  tryCatch({result <- fit_task(i);saveRDS(result,path);TRUE},error=function(e) conditionMessage(e))
 },mc.cores=N_CORES)
 if(!all(vapply(rr,isTRUE,TRUE))) stop('Task failure: ',paste(rr,collapse='; '))
 cat('Completed',max(ids),'/',nrow(tasks),'\n');flush.console()
}
all_results <- lapply(seq_len(nrow(tasks)),function(i)readRDS(file.path(outdir,'checkpoints',sprintf('task_%04d.rds',i))))
pa <- bind_rows(lapply(all_results,`[[`,'scores')); pred <- bind_rows(lapply(all_results,`[[`,'predictions'))
stopifnot(nrow(pa)==600,!anyDuplicated(pa[c('trait','model','CV','rep','fold','Environment')]))
write.csv(pa,file.path(outdir,'predictive_ability_cv1_cv2.csv'),row.names=FALSE)
write.csv(pred,file.path(outdir,'heldout_predictions.csv'),row.names=FALSE)
summary <- pa |> group_by(CV,trait,Environment,model) |> summarise(n_valid=sum(is.finite(predictive_ability)),mean_PA=mean(predictive_ability,na.rm=TRUE),median_PA=median(predictive_ability,na.rm=TRUE),sd_PA=sd(predictive_ability,na.rm=TRUE),.groups='drop')
write.csv(summary,file.path(outdir,'predictive_ability_summary.csv'),row.names=FALSE)
pa$Environment <- factor(pa$Environment,levels=ENVS,labels=ENVS)
pa$trait <- factor(pa$trait,levels=TRAITS);pa$model <- factor(pa$model,levels=c('M1','M2','M3'))
ylim <- range(c(0,pa$predictive_ability),na.rm=TRUE)+c(-.04,.04)
for(cv in c('CV1','CV2')) {
 p <- ggplot(pa[pa$CV==cv & is.finite(pa$predictive_ability),],aes(Environment,predictive_ability,fill=model))+
 geom_hline(yintercept=0,linewidth=.3,colour='grey60')+
 geom_boxplot(position=position_dodge(.8),width=.7,outlier.size=.6,linewidth=.3)+
 facet_wrap(~trait,nrow=1)+scale_fill_manual(values=c(M1='#F08A7E',M2='#4FAE62',M3='#5B9BD5'),name='Model')+
 coord_cartesian(ylim=ylim)+labs(x=NULL,y='Predictive ability',subtitle=paste(cv,'- pooled DM; random population group'))+
 theme_bw(base_size=11)+theme(strip.background=element_rect(fill='grey85'),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=45,hjust=1))
 for(ext in c('png','pdf')) ggsave(file.path(outdir,paste0('predictive_ability_',cv,'_pooled_DM.',ext)),p,width=5,height=4,dpi=300,bg='white')
}
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo.txt'))
writeLines(c('CV1: genotypes held out across all environments. CV2: per-genotype observed cells assigned to distinct folds; at least one other environment remains in training.',
 '10 repeats x 5 folds; splits shared by models. Fold-level correlations are dependent.',
 'Fixed intercept; random PopGroup, environment and line; M2 adds genomic; M3 adds genomic-by-environment.',
 'All fits use 8000 iterations, 2000 burn-in, thinning 5. Variance traces retained in task checkpoints; fit completion does not prove convergence.',
 'Explicit shared-component priors use training data only. See script header for specification.',
 'Existing stage-1 BLUEs used without precision weights; this is second-stage BLUE validation, not a full plot-level refit within each fold.',
 'These runs differ from the original fixed-PopGroup analysis in CV2 assignment, priors and chain length; they are not a controlled fixed-versus-random comparison.'),file.path(outdir,'analysis_notes.txt'))
cat('Saved CV1/CV2 plots and validation results in',outdir,'\n')
