## Convergence check for FW or MC, mirroring the DM checks (reviewer_DM_pipeline
## scripts 03 and 04): the same eight representative M3 fits (CV1, CV2, sparse-25
## with 0 and 25 shared genotypes; both phenotype sets), three independent chains
## at the main setting (8,000 iterations; chain 1 = main analysis) and at
## 80,000 iterations (burn-in 20,000, thin 10). Same R-hat functions as for DM.
## Usage: Rscript 08_representative_chain_check.R FW|MC
suppressPackageStartupMessages({library(BGLR);library(dplyr)})
trait<-commandArgs(trailingOnly=TRUE)[1];stopifnot(trait%in%c('FW','MC'))
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
here<-dirname(dirname(normalizePath(sf)))
main_out<-file.path(here,'results',trait)
txt<-readLines(file.path(here,'scripts/01_run_trait_reviewer_analysis.R'));txt<-txt[seq_len(grep("^cat\\('Running'",txt)[1]-1)]
txt<-gsub("file.path(pipeline,'results',trait)","OUTTMP",txt,fixed=TRUE)
txt<-gsub("set.seed(550000L+","set.seed(CHAIN_OFFSET+550000L+",txt,fixed=TRUE)
stopifnot(sum(grepl('OUTTMP',txt))==1,sum(grepl('CHAIN_OFFSET',txt))==1)
OUTTMP<-tempfile('chaincheck_');CHAIN_OFFSET<-0L
eval(parse(text=txt))
ids<-which(tasks$model=='M3' & tasks$split %in% c(1,51,101,151));stopifnot(length(ids)==8)
split_rhat<-function(mat){n<-nrow(mat);h<-floor(n/2);s<-cbind(mat[seq_len(h),,drop=FALSE],mat[(n-h+1):n,,drop=FALSE]);W<-mean(apply(s,2,var));B<-h*var(colMeans(s));sqrt(((h-1)/h*W+B/h)/W)}
rank_rhat<-function(mat){v<-as.vector(mat);z<-matrix(qnorm((rank(v,ties.method='average')-3/8)/(length(v)+1/4)),nrow=nrow(mat));v2<-abs(v-median(v));zf<-matrix(qnorm((rank(v2,ties.method='average')-3/8)/(length(v2)+1/4)),nrow=nrow(mat));max(split_rhat(z),split_rhat(zf))}
run_set<-function(label,chains,iter,burn,thin,chain1_main){
 dest<-file.path(here,'results/chain_checks',trait,label);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 N_ITER<<-iter;BURN<<-burn;THIN<<-thin
 for(chain in chains){
  CHAIN_OFFSET<<-as.integer(chain*1000000L+if(label=='long')500000L else 0L)
  ans<-parallel::mclapply(ids,function(i){p<-file.path(dest,sprintf('task_%04d_chain_%d.rds',i,chain));if(!file.exists(p))saveRDS(run_task(i),p);TRUE},mc.cores=8L)
  stopifnot(all(vapply(ans,isTRUE,TRUE)));cat(trait,label,'chain',chain,'done\n')
 }
 diag<-list();pc<-list()
 for(i in ids){
  r<-lapply(1:3,function(ch)if(ch==1&&chain1_main)readRDS(file.path(main_out,'checkpoints',sprintf('task_%04d.rds',i))) else readRDS(file.path(dest,sprintf('task_%04d_chain_%d.rds',i,ch))))
  for(term in names(r[[1]]$traces)){mat<-sapply(r,function(z)z$traces[[term]]);diag[[length(diag)+1]]<-data.frame(trait=trait,setting=label,task=i,tasks[i,],term=term,split_Rhat=split_rhat(mat),rank_split_Rhat=rank_rhat(mat))}
  for(site in c('ABR33','JKI')){rr<-sapply(r,function(z)z$scores$PA[z$scores$Environment==site]);pp<-sapply(r,function(z)z$predictions$predicted[z$predictions$Environment==site])
   pc[[length(pc)+1]]<-data.frame(trait=trait,setting=label,task=i,tasks[i,],Environment=site,PA_range=diff(range(rr)),minimum_prediction_correlation=min(cor(pp)[lower.tri(cor(pp))]))}
 }
 d<-bind_rows(diag);p<-bind_rows(pc)
 write.csv(d,file.path(dest,'variance_chain_diagnostics.csv'),row.names=FALSE);write.csv(p,file.path(dest,'prediction_stability.csv'),row.names=FALSE)
 cat(trait,label,': max rank-normalised split Rhat =',round(max(d$rank_split_Rhat),4),'; max PA range =',round(max(p$PA_range,na.rm=TRUE),3),'; min prediction correlation =',round(min(p$minimum_prediction_correlation,na.rm=TRUE),3),'\n')
}
run_set('main',2:3,8000L,2000L,5L,TRUE)
run_set('long',1:3,80000L,20000L,10L,FALSE)
