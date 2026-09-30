## Independent sampling chains for representative M3 fits, same data and priors.
## Reports conventional and rank-normalized split-Rhat and prediction stability;
## these checks cover selected cases, not all 960 models.
suppressPackageStartupMessages({library(BGLR);library(dplyr);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE);self<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
main<-file.path(dirname(normalizePath(self)),'01_run_DM_reviewer_analysis.R')
e<-new.env(parent=globalenv());e$sf<-main;e$pipeline<-dirname(dirname(main));e$project<-dirname(e$pipeline);e$out<-file.path(e$pipeline,'results')
# Evaluate setup/functions only, skipping writes, path detection and run loop.
for(x in parse(main)){
 txt<-paste(deparse(x),collapse=' ')
 if(grepl('^for \\(start',txt))break
 if(grepl('^cat\\(',txt)||grepl('^write\\.',txt)||grepl('^dir.create',txt))next
 if(is.call(x)&&identical(x[[1]],as.name('<-'))&&is.symbol(x[[2]])&&as.character(x[[2]])%in%c('a','sf','pipeline','project','out','inputs','sig','sg'))next
 if(grepl('^if \\(file.exists\\(sg',txt))next
 eval(x,e)
}
# CV1, CV2, disjoint/full-overlap sparse25, both target branches, model M3.
ids<-which(e$tasks$model=='M3' & e$tasks$split %in% c(1,51,101,151))
dest<-file.path(e$out,'chain_checks');dir.create(dest,showWarnings=FALSE)
orig_body<-body(e$run_task)
for(chain in 2:3){
 body(e$run_task)<-parse(text=gsub('550000L',paste0(550000L+chain*1000000L,'L'),paste(deparse(orig_body),collapse='\n')))[[1]]
 results<-parallel::mclapply(ids,function(i){path<-file.path(dest,paste0('task_',i,'_chain_',chain,'.rds'));if(!file.exists(path))saveRDS(e$run_task(i),path);TRUE},mc.cores=4)
 stopifnot(all(vapply(results,isTRUE,TRUE)));cat('Finished independent chain',chain,'for',length(ids),'fits\n')
}
split_rhat<-function(mat){n<-nrow(mat);h<-floor(n/2);s<-cbind(mat[seq_len(h),,drop=FALSE],mat[(n-h+1):n,,drop=FALSE]);W<-mean(apply(s,2,var));B<-h*var(colMeans(s));sqrt(((h-1)/h*W+B/h)/W)}
rank_rhat<-function(mat){v<-as.vector(mat);z<-matrix(qnorm((rank(v,ties.method='average')-3/8)/(length(v)+1/4)),nrow=nrow(mat));v2<-abs(v-median(v));zf<-matrix(qnorm((rank(v2,ties.method='average')-3/8)/(length(v2)+1/4)),nrow=nrow(mat));max(split_rhat(z),split_rhat(zf))}
diag<-list();predcheck<-list()
for(i in ids){
 r<-list(readRDS(file.path(e$out,'checkpoints',sprintf('task_%04d.rds',i))),readRDS(file.path(dest,paste0('task_',i,'_chain_2.rds'))),readRDS(file.path(dest,paste0('task_',i,'_chain_3.rds'))))
 for(term in names(r[[1]]$traces)){
  mat<-sapply(r,function(z)z$traces[[term]])
  diag[[length(diag)+1]]<-data.frame(task=i,e$tasks[i,],term=term,split_Rhat=split_rhat(mat),rank_split_Rhat=rank_rhat(mat))
 }
 for(site in c('ABR33','JKI')){
  rr<-sapply(r,function(z)z$scores$PA[z$scores$Environment==site])
  pp<-sapply(r,function(z)z$predictions$predicted[z$predictions$Environment==site])
  predcheck[[length(predcheck)+1]]<-data.frame(task=i,e$tasks[i,],Environment=site,PA_chain1=rr[1],PA_chain2=rr[2],PA_chain3=rr[3],PA_range=diff(range(rr)),minimum_prediction_correlation=min(cor(pp)[lower.tri(cor(pp))]))
 }
}
d<-bind_rows(diag);pc<-bind_rows(predcheck)
write.csv(d,file.path(dest,'variance_chain_diagnostics.csv'),row.names=FALSE);write.csv(pc,file.path(dest,'prediction_stability.csv'),row.names=FALSE)
cat('Maximum rank-normalized split Rhat:',max(d$rank_split_Rhat),'\n');cat('Maximum PA range across chains:',max(pc$PA_range),'\n');print(pc)
