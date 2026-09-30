## Reviewer DM analysis: paired raw/group-adjusted prediction on TWO pooled trials.
## Run with Rscript; original scripts and results are read-only inputs.
## Group adjustment is learned exclusively from TRAINING genotype/trial records.
## BGLR has NO POP term. Groups are used only for preprocessing/baseline/evaluation.
## Both target branches have identical partitions and numerical priors based on
## raw training variance. No test value sets a group mean, prior or model parameter.
## Stage-1 pooled BLUEs are pre-estimated; this is second-stage validation.
suppressPackageStartupMessages({library(BGLR);library(dplyr);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE)
sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
pipeline<-dirname(dirname(normalizePath(sf))); project<-dirname(pipeline)
out<-file.path(pipeline,'results');dir.create(file.path(out,'checkpoints'),recursive=TRUE,showWarnings=FALSE)
source_dir<-file.path(project,'results/site_pooled_revision')
N_ITER<-8000L;BURN<-2000L;THIN<-5L;CORES<-4L
b<-read.csv(file.path(source_dir,'pooled_DM_BLUEs.csv'))
pops<-read.csv(file.path(project,'data/popgroups_125.csv'))
checks<-read.csv(file.path(source_dir,'pooled_model_checks.csv'))
stopifnot(all(checks$converged),all(checks$stable),all(b$converged),all(b$trait=='DM'),!anyDuplicated(b[c('Genotype','Environment')]),!anyDuplicated(pops$Genotype))
ge<-new.env();load(file.path(project,'data/markers/Gmatrix_125.rda'),envir=ge);G<-ge$G
gids<-rownames(G);sites<-c('ABR33','JKI')
stopifnot(identical(gids,colnames(G)),!anyDuplicated(gids),all(is.finite(G)),max(abs(G-t(G)))<1e-8)
obs<-expand.grid(Genotype=gids,Environment=sites,stringsAsFactors=FALSE)
obs$PopGroup<-pops$PopGroup[match(obs$Genotype,pops$Genotype)]
stopifnot(!anyNA(obs$PopGroup))
key<-paste(obs$Genotype,obs$Environment)
y<-b$BLUE[match(key,paste(b$Genotype,b$Environment))]
Zl<-model.matrix(~factor(obs$Genotype,levels=gids)-1);Ze<-model.matrix(~factor(obs$Environment,levels=sites)-1)
Kg<-Zl%*%G%*%t(Zl);Kge<-Kg*tcrossprod(Ze)
eig<-function(k){e<-eigen(k,symmetric=TRUE);stopifnot(min(e$values)>-1e-7);j<-e$values>1e-10;list(V=e$vectors[,j,drop=FALSE],d=e$values[j])}
eg<-eig(Kg);ege<-eig(Kge)
## Reuse the exact site-pooled CV and sparse-25 assignments for comparability.
fcv<-read.csv(file.path(source_dir,'cv/fold_assignments.csv'))
fsp<-read.csv(file.path(source_dir,'sparse25/training_designs.csv'))
splits<-list()
for(cv in c('CV1','CV2'))for(r in 1:10)for(k in 1:5){z<-fcv[fcv$CV==cv & fcv$rep==r,];ff<-z$fold[match(key,paste(z$Genotype,z$Environment))];stopifnot(!anyNA(ff));splits[[length(splits)+1]]<-list(scheme=cv,rep=r,fold=k,overlap=NA_integer_,test=ff==k & is.finite(y),train=ff!=k & is.finite(y))}
for(o in seq(0,25,5))for(r in 1:10){z<-fsp[fsp$CV==paste0('Overlap',o)&fsp$rep==r,];ff<-z$fold[match(key,paste(z$Genotype,z$Environment))];stopifnot(!anyNA(ff));splits[[length(splits)+1]]<-list(scheme='Sparse25',rep=r,fold=1L,overlap=o,test=ff==1 & is.finite(y),train=ff==0 & is.finite(y))}
group_key<-paste(obs$Environment,obs$PopGroup,sep='::')
for(i in seq_along(splits)) {
 s<-splits[[i]];stopifnot(!any(s$train&s$test))
 if(s$scheme=='CV1')stopifnot(!any(obs$Genotype[s$test]%in%obs$Genotype[s$train]))
 if(s$scheme=='CV2')stopifnot(all(obs$Genotype[s$test]%in%obs$Genotype[s$train]))
 if(s$scheme=='Sparse25'){
  stopifnot(all(table(factor(obs$Environment[s$train],levels=sites))==25))
  stopifnot(length(intersect(obs$Genotype[s$train&obs$Environment=='ABR33'],obs$Genotype[s$train&obs$Environment=='JKI']))==s$overlap)
 }
 means<-tapply(y[s$train],group_key[s$train],mean)
 s$offset<-unname(means[group_key]);stopifnot(all(is.finite(s$offset[is.finite(y)])))
 s$adjusted<-y-s$offset
 stopifnot(max(abs(tapply(s$adjusted[s$train],group_key[s$train],mean)))<1e-10)
 ## Regression coefficients must not change if all held-out values are perturbed.
 yp<-y;yp[s$test]<-yp[s$test]+1000
 stopifnot(identical(means,tapply(yp[s$train],group_key[s$train],mean)))
 splits[[i]]<-s
}
assignment<-bind_rows(lapply(seq_along(splits),function(i){s<-splits[[i]];data.frame(split=i,scheme=s$scheme,rep=s$rep,fold=s$fold,overlap=s$overlap,obs,train=s$train,test=s$test,training_group_mean=s$offset,raw=y,adjusted=s$adjusted)}))
write.csv(assignment,file.path(out,'split_assignments_and_adjustments.csv'),row.names=FALSE)
## Descriptive full-data residuals ONLY: exact group-zero check, never model input.
descriptive<-data.frame(obs,BLUE=y) |> filter(is.finite(BLUE)) |> group_by(Environment,PopGroup) |> mutate(group_mean=mean(BLUE),full_data_residual=BLUE-group_mean) |> ungroup()
write.csv(descriptive,file.path(out,'descriptive_group_residuals_not_CV_inputs.csv'),row.names=FALSE)
group_summary<-descriptive |> group_by(Environment,PopGroup) |> summarise(n=n(),mean_BLUE=mean(BLUE),sd_BLUE=sd(BLUE),mean_residual=mean(full_data_residual),.groups='drop')
write.csv(group_summary,file.path(out,'group_means_before_after.csv'),row.names=FALSE)
r2<-bind_rows(lapply(sites,function(site){z<-descriptive[descriptive$Environment==site,];data.frame(Environment=site,n=nrow(z),group_R2=summary(lm(BLUE~PopGroup,z))$r.squared)}))
write.csv(r2,file.path(out,'descriptive_population_R2.csv'),row.names=FALSE)
tasks<-expand.grid(split=seq_along(splits),model=c('M1','M2','M3'),target=c('Raw','Group_adjusted'),stringsAsFactors=FALSE)
inputs<-c(file.path(source_dir,'pooled_DM_BLUEs.csv'),file.path(source_dir,'cv/fold_assignments.csv'),file.path(source_dir,'sparse25/training_designs.csv'),file.path(project,'data/markers/Gmatrix_125.rda'),file.path(project,'data/popgroups_125.csv'),normalizePath(sf))
sig<-list(md5=tools::md5sum(inputs),iterations=c(N_ITER,BURN,THIN),tasks=tasks)
sg<-file.path(out,'run_signature.rds');if(file.exists(sg))stopifnot(identical(readRDS(sg),sig))else saveRDS(sig,sg)
write.csv(tasks,file.path(out,'tasks.csv'),row.names=FALSE)
cor_safe<-function(a,b){if(length(a)<5||sd(a)<1e-10||sd(b)<1e-10)NA_real_ else cor(a,b)}
score<-function(z){data.frame(n=nrow(z),PA=cor_safe(z$observed,z$predicted),raw_scale_PA=cor_safe(z$observed_raw,z$predicted_raw),RMSE=sqrt(mean((z$observed-z$predicted)^2)),baseline_RMSE=sqrt(mean((z$observed_raw-z$group_baseline)^2)))}
run_task<-function(i){
 t<-tasks[i,];s<-splits[[t$split]];response<-if(t$target=='Raw')y else s$adjusted
 yt<-response;yt[!s$train]<-NA_real_;stopifnot(all(is.na(yt[s$test])))
 v<-var(y[s$train]);stopifnot(is.finite(v),v>0)
 brr<-function(x)list(X=x,model='BRR',df0=5,S0=7*.1*v/mean(rowSums(x*x)))
 rk<-function(e,k)list(V=e$V,d=e$d,model='RKHS',df0=5,S0=7*.1*v/mean(diag(k)))
 ETA<-list(ENV=brr(Ze),LINE=brr(Zl))
 if(t$model!='M1')ETA$GENO<-rk(eg,Kg)
 if(t$model=='M3')ETA$GXE<-rk(ege,Kge)
 td<-tempfile('reviewer_dm_');dir.create(td);on.exit(unlink(td,recursive=TRUE))
 ## Same RNG seed for the paired target branches of a given split/model.
 set.seed(550000L+10L*t$split+match(t$model,c('M1','M2','M3')))
 fit<-BGLR(y=yt,ETA=ETA,df0=5,S0=7*.5*v,nIter=N_ITER,burnIn=BURN,thin=THIN,verbose=FALSE,saveAt=file.path(td,'bg_'))
 line<-fit$ETA$LINE$b
 ## Independent line effects for wholly unobserved genotypes have exact mean 0;
 ## remove simulation noise rather than scoring a spurious CV1 M1 ranking.
 unseen<-!gids%in%obs$Genotype[s$train];line[unseen]<-0
 pred<-as.numeric(fit$mu+Ze%*%fit$ETA$ENV$b+Zl%*%line)
 if(t$model!='M1')pred<-pred+fit$ETA$GENO$u
 if(t$model=='M3')pred<-pred+fit$ETA$GXE$u
 stopifnot(all(is.finite(pred)))
 rawpred<-if(t$target=='Raw')pred else pred+s$offset
 pr<-data.frame(split=t$split,scheme=s$scheme,rep=s$rep,fold=s$fold,overlap=s$overlap,target=t$target,model=t$model,obs[s$test,],observed=response[s$test],predicted=pred[s$test],observed_raw=y[s$test],predicted_raw=rawpred[s$test],group_baseline=s$offset[s$test],row.names=NULL)
 sc<-bind_rows(lapply(sites,function(e){z<-pr[pr$Environment==e,];data.frame(split=t$split,scheme=s$scheme,rep=s$rep,fold=s$fold,overlap=s$overlap,target=t$target,model=t$model,Environment=e,score(z))}))
 files<-list.files(td,pattern='(varB|varU|varE)\\.dat$',full.names=TRUE)
 traces<-setNames(lapply(files,function(f)scan(f,quiet=TRUE)[-seq_len(BURN/THIN)]),basename(files))
 stopifnot(all(lengths(traces)==(N_ITER-BURN)/THIN))
 list(predictions=pr,scores=sc,traces=traces)
}
cat('Running',nrow(tasks),'paired DM fits on two physical trials\n')
for(start in seq(1,nrow(tasks),by=40)){
 ids<-start:min(start+39,nrow(tasks))
 ans<-parallel::mclapply(ids,function(i){path<-file.path(out,'checkpoints',sprintf('task_%04d.rds',i));if(file.exists(path))return(TRUE);tryCatch({saveRDS(run_task(i),path);TRUE},error=function(e)conditionMessage(e))},mc.cores=CORES)
 if(!all(vapply(ans,isTRUE,TRUE)))stop(paste(ans,collapse='; '))
 cat('Completed',max(ids),'/',nrow(tasks),'\n');flush.console()
}
res<-lapply(seq_len(nrow(tasks)),function(i)readRDS(file.path(out,'checkpoints',sprintf('task_%04d.rds',i))))
pr<-bind_rows(lapply(res,`[[`,'predictions'));sc<-bind_rows(lapply(res,`[[`,'scores'))
write.csv(pr,file.path(out,'heldout_predictions.csv'),row.names=FALSE)
write.csv(sc,file.path(out,'fold_scores.csv'),row.names=FALSE)
## Group-only predictor fitted on training data; raw-scale PA and RMSE baseline.
base<-pr |> filter(target=='Raw',model=='M1') |> group_by(scheme,rep,fold,overlap,Environment) |> summarise(n=n(),PA=cor_safe(observed_raw,group_baseline),RMSE=sqrt(mean((observed_raw-group_baseline)^2)),.groups='drop')
write.csv(base,file.path(out,'group_only_baseline.csv'),row.names=FALSE)
## Within-group metrics on RAW held-out observations/predictions, so target
## offsets do not change the outcome being evaluated. Group centering here is
## for SCORING ONLY, not for fitting a model or tuning a hyperparameter.
repmetrics<-pr |> group_by(scheme,rep,overlap,target,model,Environment) |> group_modify(function(z,key){
 ac<-z$observed_raw-ave(z$observed_raw,z$PopGroup,FUN=mean)
 pc<-z$predicted_raw-ave(z$predicted_raw,z$PopGroup,FUN=mean)
 data.frame(n=nrow(z),PA_target=cor_safe(z$observed,z$predicted),PA_raw=cor_safe(z$observed_raw,z$predicted_raw),PA_within_pooled=cor_safe(ac,pc),RMSE_raw=sqrt(mean((z$observed_raw-z$predicted_raw)^2)),RMSE_group_only=sqrt(mean((z$observed_raw-z$group_baseline)^2)))
}) |> ungroup()
bygroup<-pr |> group_by(scheme,rep,overlap,target,model,Environment,PopGroup) |> summarise(n=n(),PA=cor_safe(observed_raw,predicted_raw),RMSE=sqrt(mean((observed_raw-predicted_raw)^2)),.groups='drop')
## Do not give the no-genomics/no-group CV1 model a ranking based on fold offsets.
repmetrics[repmetrics$scheme=='CV1'&repmetrics$model=='M1',c('PA_target','PA_raw','PA_within_pooled')]<-NA_real_
bygroup$PA[bygroup$scheme=='CV1'&bygroup$model=='M1']<-NA_real_
write.csv(repmetrics,file.path(out,'per_repeat_metrics.csv'),row.names=FALSE)
write.csv(bygroup,file.path(out,'within_group_per_repeat.csv'),row.names=FALSE)
mean_safe<-function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE)
summary<-repmetrics |> group_by(scheme,overlap,target,model,Environment) |> summarise(across(c(PA_target,PA_raw,PA_within_pooled,RMSE_raw,RMSE_group_only),mean_safe),.groups='drop')
write.csv(summary,file.path(out,'prediction_summary.csv'),row.names=FALSE)
write.csv(bygroup |> group_by(scheme,overlap,target,model,Environment,PopGroup) |> summarise(n_min=min(n),n_repeats=sum(is.finite(PA)),mean_PA=mean_safe(PA),mean_RMSE=mean(RMSE),.groups='drop'),file.path(out,'within_group_summary.csv'),row.names=FALSE)
cols<-c(M1='#F08A7E',M2='#4FAE62',M3='#5B9BD5')
for(cv in c('CV1','CV2')) {
 d<-sc[sc$scheme==cv & is.finite(sc$PA),];d$target<-factor(d$target,levels=c('Raw','Group_adjusted'),labels=c('Group means retained','Training-only group adjustment'))
 p<-ggplot(d,aes(Environment,PA,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.7)+facet_wrap(~target,nrow=1)+scale_fill_manual(values=cols,name='Model')+labs(x=NULL,y='Predictive ability',subtitle=paste(cv,'- pooled DM; no POP term in BGLR'),caption=if(cv=='CV1')'M1 cannot rank unseen genotypes; its correlation is undefined.' else NULL)+theme_bw(base_size=11)
 for(ext in c('png','pdf'))ggsave(file.path(out,paste0(cv,'_paired_targets.',ext)),p,width=9,height=4,dpi=300,bg='white')
 d<-repmetrics[repmetrics$scheme==cv & is.finite(repmetrics$PA_within_pooled),]
 p<-ggplot(d,aes(Environment,PA_within_pooled,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.7)+facet_wrap(~target,nrow=1)+scale_fill_manual(values=cols,name='Model')+labs(x=NULL,y='Within-group pooled correlation',subtitle=paste(cv,'- pooled DM; 10 repetitions'),caption='Both observations and predictions centered by group for scoring only.')+theme_bw(base_size=11)
 for(ext in c('png','pdf'))ggsave(file.path(out,paste0(cv,'_within_group.',ext)),p,width=9,height=4,dpi=300,bg='white')
}
d<-bygroup[bygroup$scheme%in%c('CV1','CV2') & bygroup$target=='Group_adjusted' & is.finite(bygroup$PA),]
p<-ggplot(d,aes(Environment,PA,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.5)+facet_grid(scheme~PopGroup)+scale_fill_manual(values=cols,name='Model')+labs(x=NULL,y='Within-group predictive ability',subtitle='DM - training-only group adjustment; 10 repetitions')+theme_bw(base_size=10)
for(ext in c('png','pdf'))ggsave(file.path(out,paste0('species_specific_prediction.',ext)),p,width=11,height=6,dpi=300,bg='white')
d<-repmetrics[repmetrics$scheme=='Sparse25'&repmetrics$target=='Group_adjusted',];d$overlap<-factor(d$overlap,levels=seq(0,25,5))
p<-ggplot(d,aes(overlap,PA_within_pooled,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.5)+facet_wrap(~Environment,nrow=1)+scale_fill_manual(values=cols,name='Model')+labs(x='Shared training genotypes (25 per trial)',y='Within-group pooled correlation',subtitle='DM - training-only group adjustment')+theme_bw(base_size=11)
for(ext in c('png','pdf'))ggsave(file.path(out,paste0('sparse25_within_group.',ext)),p,width=9,height=4,dpi=300,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Saved paired target, within-group, and sparse-25 outputs in',out,'\n')
