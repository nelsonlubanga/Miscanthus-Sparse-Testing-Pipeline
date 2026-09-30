## Two physical trials; 25 observed training genotypes per trial, overlap 0:25 by 5.
## Unique training genotypes = 50 - overlap. Ten repeats, paired across models.
## Training budgets, rather than unique genotype counts, are held constant.
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
outdir <- file.path(root,'results','site_pooled_revision','sparse25')
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
common <- gids[vapply(gids,function(g) all(is.finite(ys[['DM']][obs$Genotype==g])),TRUE)]
for(o in seq(0,25,5)) for(r in seq_len(N_REPS)) {
 set.seed(800000L+100L*o+r)
 perm <- sample(common)
 shared <- if(o>0)perm[seq_len(o)] else character()
 a <- c(shared,perm[o+seq_len(25-o)])
 bb <- c(shared,perm[25+seq_len(25-o)])
 ## seq_len(0) contributes no entries at full overlap.
 train <- (obs$Environment=='ABR33' & obs$Genotype %in% a) | (obs$Environment=='JKI' & obs$Genotype %in% bb)
 stopifnot(length(intersect(a,bb))==o,sum(train & obs$Environment=='ABR33')==25,sum(train & obs$Environment=='JKI')==25,all(is.finite(ys[['DM']][train])))
 f <- ifelse(train,0L,1L)
 cv<-paste0('Overlap',o)
 key<-paste('DM',cv,r,sep='_');folds[[key]]<-f
 fold_audit[[key]]<-data.frame(trait='DM',CV=cv,rep=r,obs,fold=f,observed=is.finite(ys[['DM']]))
}
write.csv(bind_rows(fold_audit),file.path(outdir,'training_designs.csv'),row.names=FALSE)
saveRDS(folds,file.path(outdir,'training_designs.rds'))
tasks<-expand.grid(trait='DM',model=c('M1','M2','M3'),CV=paste0('Overlap',seq(0,25,5)),rep=seq_len(N_REPS),fold=1L,stringsAsFactors=FALSE)
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
stopifnot(nrow(pa)==360,!anyDuplicated(pa[c('trait','model','CV','rep','fold','Environment')]))
write.csv(pa,file.path(outdir,'predictive_ability_cv1_cv2.csv'),row.names=FALSE)
write.csv(pred,file.path(outdir,'heldout_predictions.csv'),row.names=FALSE)
summary <- pa |> group_by(CV,trait,Environment,model) |> summarise(n_valid=sum(is.finite(predictive_ability)),mean_PA=mean(predictive_ability,na.rm=TRUE),median_PA=median(predictive_ability,na.rm=TRUE),sd_PA=sd(predictive_ability,na.rm=TRUE),.groups='drop')
write.csv(summary,file.path(outdir,'predictive_ability_summary.csv'),row.names=FALSE)
pa$Environment<-factor(pa$Environment,levels=ENVS)
pa$overlap<-factor(sub('Overlap','',pa$CV),levels=as.character(seq(0,25,5)))
pa$model<-factor(pa$model,levels=c('M1','M2','M3'))
p<-ggplot(pa,aes(overlap,predictive_ability,fill=model))+
 geom_hline(yintercept=0,linewidth=.3,colour='grey60')+
 geom_boxplot(position=position_dodge(.8),width=.7,outlier.size=.6,linewidth=.3)+
 facet_wrap(~Environment,nrow=1)+scale_fill_manual(values=c(M1='#F08A7E',M2='#4FAE62',M3='#5B9BD5'),name='Model')+
 labs(x='Training genotypes shared between trials (25 per trial)',y='Predictive ability',subtitle='Pooled DM - random population group')+
 theme_bw(base_size=11)+theme(strip.background=element_rect(fill='grey85'),panel.grid.minor=element_blank())
for(ext in c('png','pdf'))ggsave(file.path(outdir,paste0('sparse25_pooled_DM.',ext)),p,width=9,height=3.8,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo.txt'))
cat('Saved sparse-25 predictions and plots in',outdir,'\n')
