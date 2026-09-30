## Monte Carlo variability at the level reported: pooled within-group predictive
## ability per repetition and its mean over 10 repetitions, for M3 under CV1/CV2.
## Chain 1 = the main-analysis checkpoints; chains 2 and 3 rerun the same tasks
## with different seeds and otherwise identical data, priors and settings.
## Reuses the setup and run_task() of the main scripts (01) without modifying them;
## setup outputs go to a temporary directory. Usage: Rscript 06_chain_variability.R DM|FW|MC [M3|M2]
suppressPackageStartupMessages({library(BGLR);library(dplyr)})
trait<-commandArgs(trailingOnly=TRUE)[1];stopifnot(trait%in%c('DM','FW','MC'))
MODEL<-commandArgs(trailingOnly=TRUE)[2];if(is.na(MODEL))MODEL<-'M3';stopifnot(MODEL%in%c('M2','M3'))
suffix<-if(MODEL=='M3')'' else paste0('_',MODEL)
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
here<-dirname(dirname(normalizePath(sf)));project<-dirname(here)
main_src<-if(trait=='DM')file.path(project,'reviewer_DM_pipeline/scripts/01_run_DM_reviewer_analysis.R') else file.path(here,'scripts/01_run_trait_reviewer_analysis.R')
main_out<-if(trait=='DM')file.path(project,'reviewer_DM_pipeline/results') else file.path(here,'results',trait)
outdir<-file.path(here,'results/chain_variability');dir.create(file.path(outdir,trait),recursive=TRUE,showWarnings=FALSE)
txt<-readLines(main_src);txt<-txt[seq_len(grep("^cat\\('Running'",txt)[1]-1)]
txt<-gsub("file.path(pipeline,'results',trait)","OUTTMP",txt,fixed=TRUE)
txt<-gsub("file.path(pipeline,'results')","OUTTMP",txt,fixed=TRUE)
txt<-gsub("set.seed(550000L+","set.seed(CHAIN_OFFSET+550000L+",txt,fixed=TRUE)
stopifnot(sum(grepl('OUTTMP',txt))==1,sum(grepl('CHAIN_OFFSET',txt))==1)
OUTTMP<-tempfile('chainvar_');CHAIN_OFFSET<-0L
eval(parse(text=txt))
ids<-which(tasks$model==MODEL & vapply(tasks$split,function(i)splits[[i]]$scheme,'')%in%c('CV1','CV2'))
stopifnot(length(ids)==200)
for(chain in 2:3) {
 CHAIN_OFFSET<-as.integer(chain*1000000L)
 todo<-ids[!file.exists(file.path(outdir,trait,sprintf('task_%04d_chain_%d.rds',ids,chain)))]
 ans<-parallel::mclapply(todo,function(i){tryCatch({saveRDS(run_task(i)$predictions,file.path(outdir,trait,sprintf('task_%04d_chain_%d.rds',i,chain)));TRUE},error=function(e)conditionMessage(e))},mc.cores=10L)
 if(!all(vapply(ans,isTRUE,TRUE)))stop(paste(ans,collapse='; '))
 cat(trait,MODEL,'chain',chain,'done\n')
}
pr<-bind_rows(lapply(ids,function(i)data.frame(chain=1L,readRDS(file.path(main_out,'checkpoints',sprintf('task_%04d.rds',i)))$predictions)),
 lapply(2:3,function(ch)bind_rows(lapply(ids,function(i)data.frame(chain=ch,readRDS(file.path(outdir,trait,sprintf('task_%04d_chain_%d.rds',i,ch))))))))
rep_pa<-pr |> group_by(chain,scheme,rep,target,Environment) |> summarise(PA=cor(observed_raw-ave(observed_raw,PopGroup),predicted_raw-ave(predicted_raw,PopGroup)),.groups='drop')
## Check chain 1 reproduces the reported per-repetition values exactly.
rep_main<-read.csv(file.path(main_out,'per_repeat_metrics.csv')) |> filter(model==MODEL,scheme%in%c('CV1','CV2'))
z<-inner_join(rep_pa[rep_pa$chain==1,],rep_main,by=c('scheme','rep','target','Environment'));stopifnot(nrow(z)==80,max(abs(z$PA-z$PA_within_pooled))<1e-10)
per_rep<-rep_pa |> group_by(scheme,rep,target,Environment) |> summarise(range=diff(range(PA)),.groups='drop')
means<-rep_pa |> group_by(chain,scheme,target,Environment) |> summarise(mean_PA=mean(PA),.groups='drop')
mean_range<-means |> group_by(scheme,target,Environment) |> summarise(min=min(mean_PA),max=max(mean_PA),range=max-min,.groups='drop')
paired<-rep_pa |> tidyr::pivot_wider(names_from=target,values_from=PA) |> mutate(diff=Raw-Group_adjusted) |> group_by(chain,scheme,Environment) |> summarise(mean_diff=mean(diff),.groups='drop') |> group_by(scheme,Environment) |> summarise(min_diff=min(mean_diff),max_diff=max(mean_diff),range=max_diff-min_diff,.groups='drop')
write.csv(rep_pa,file.path(outdir,paste0(trait,suffix,'_per_repeat_by_chain.csv')),row.names=FALSE)
write.csv(mean_range,file.path(outdir,paste0(trait,suffix,'_mean_over_repeats_by_chain.csv')),row.names=FALSE)
write.csv(paired,file.path(outdir,paste0(trait,suffix,'_adjustment_difference_by_chain.csv')),row.names=FALSE)
cat(trait,MODEL,': max per-repetition range across chains =',round(max(per_rep$range),3),'; max range of 10-repetition mean =',round(max(mean_range$range),3),'; max range of adjusted-unadjusted difference =',round(max(paired$range),3),'\n')
