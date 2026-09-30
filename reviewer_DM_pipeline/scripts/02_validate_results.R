suppressPackageStartupMessages(library(dplyr))
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE);p<-file.path(dirname(dirname(normalizePath(sf))),'results')
f<-read.csv(file.path(p,'split_assignments_and_adjustments.csv'))
pr<-read.csv(file.path(p,'heldout_predictions.csv'));sc<-read.csv(file.path(p,'fold_scores.csv'))
stopifnot(nrow(sc)==1920,!anyDuplicated(pr[c('split','target','model','Genotype','Environment')]),all(is.finite(pr$predicted)),all(is.finite(pr$observed)))
cor_safe<-function(a,b)if(length(a)<5||sd(a)<1e-10||sd(b)<1e-10)NA_real_ else cor(a,b)
for(i in unique(f$split)){
 a<-f[f$split==i,];train<-a[a$train,];test<-a[a$test,]
 means<-tapply(train$raw,paste(train$Environment,train$PopGroup),mean)
 off<-unname(means[paste(a$Environment,a$PopGroup)])
 stopifnot(max(abs(off-a$training_group_mean),na.rm=TRUE)<1e-10,max(abs(tapply(train$adjusted,paste(train$Environment,train$PopGroup),mean)))<1e-10)
 if(a$scheme[1]=='CV1')stopifnot(!any(test$Genotype%in%train$Genotype))
 if(a$scheme[1]=='CV2')stopifnot(all(test$Genotype%in%train$Genotype))
 if(a$scheme[1]=='Sparse25')stopifnot(all(table(train$Environment)==25),length(intersect(train$Genotype[train$Environment=='ABR33'],train$Genotype[train$Environment=='JKI']))==a$overlap[1])
 for(target in c('Raw','Group_adjusted'))for(model in c('M1','M2','M3')){
  z<-pr[pr$split==i & pr$target==target & pr$model==model,];j<-match(paste(z$Genotype,z$Environment),paste(test$Genotype,test$Environment));stopifnot(nrow(z)==nrow(test),!anyNA(j))
  expected<-if(target=='Raw')test$raw[j] else test$raw[j]-test$training_group_mean[j]
  stopifnot(max(abs(z$observed-expected))<1e-10)
 }
}
recalc<-pr |> group_by(split,target,model,Environment) |> summarise(r=cor_safe(observed,predicted),rmse=sqrt(mean((observed-predicted)^2)),n=n(),.groups='drop')
z<-left_join(sc,recalc,by=c('split','target','model','Environment'))
stopifnot(identical(is.na(z$PA),is.na(z$r)),all(abs(z$PA-z$r)<1e-10,na.rm=TRUE),all(abs(z$RMSE-z$rmse)<1e-10),all(z$n.x==z$n.y))
stopifnot(all(is.na(sc$PA[sc$scheme=='CV1'&sc$model=='M1'])))
cat('PASS: 960 fits, 1920 trial scores; target pairing, independent training-only adjustment, group means, physical-trial masks, sparse-25 budgets, held-out targets, correlations and RMSE verified.\n')
cat('Undefined correlations from constant predictions:',sum(is.na(sc$PA)),'\n')
writeLines(c('All checks passed.',paste('Undefined constant-prediction correlations:',sum(is.na(sc$PA)))) ,file.path(p,'validation_status.txt'))
