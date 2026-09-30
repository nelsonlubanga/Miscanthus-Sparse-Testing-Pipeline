## Impact of the corrected G matrix on objective (iii): M2 and M3 under CV1/CV2,
## both phenotype sets, identical splits, priors and seeds to the reported results.
## Only the G file differs (Gmatrix_125_corrected.rda). Main results are not modified.
## Usage: Rscript 10_corrected_G_impact.R DM|FW|MC
suppressPackageStartupMessages({library(BGLR);library(dplyr)})
trait<-commandArgs(trailingOnly=TRUE)[1];stopifnot(trait%in%c('DM','FW','MC'))
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
here<-dirname(dirname(normalizePath(sf)));project<-dirname(here)
main_src<-if(trait=='DM')file.path(project,'reviewer_DM_pipeline/scripts/01_run_DM_reviewer_analysis.R') else file.path(here,'scripts/01_run_trait_reviewer_analysis.R')
main_out<-if(trait=='DM')file.path(project,'reviewer_DM_pipeline/results') else file.path(here,'results',trait)
outdir<-file.path(here,'results/corrected_G',trait);dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
txt<-readLines(main_src);txt<-txt[seq_len(grep("^cat\\('Running'",txt)[1]-1)]
txt<-gsub("file.path(pipeline,'results',trait)","OUTTMP",txt,fixed=TRUE)
txt<-gsub("file.path(pipeline,'results')","OUTTMP",txt,fixed=TRUE)
txt<-gsub("data/markers/Gmatrix_125.rda","data/markers/Gmatrix_125_corrected.rda",txt,fixed=TRUE)
stopifnot(sum(grepl('OUTTMP',txt))==1,sum(grepl('Gmatrix_125_corrected',txt))>=1)
OUTTMP<-tempfile('corrG_');eval(parse(text=txt))
ids<-which(tasks$model%in%c('M2','M3') & vapply(tasks$split,function(i)splits[[i]]$scheme,'')%in%c('CV1','CV2'))
stopifnot(length(ids)==400)
todo<-ids[!file.exists(file.path(outdir,sprintf('task_%04d.rds',ids)))]
ans<-parallel::mclapply(todo,function(i){tryCatch({saveRDS(run_task(i)$predictions,file.path(outdir,sprintf('task_%04d.rds',i)));TRUE},error=function(e)conditionMessage(e))},mc.cores=10L)
if(!all(vapply(ans,isTRUE,TRUE)))stop(paste(ans,collapse='; '))
pr_new<-bind_rows(lapply(ids,function(i)readRDS(file.path(outdir,sprintf('task_%04d.rds',i)))))
pr_old<-bind_rows(lapply(ids,function(i)readRDS(file.path(main_out,'checkpoints',sprintf('task_%04d.rds',i)))$predictions))
metr<-function(pr)pr |> group_by(scheme,rep,target,model,Environment) |> summarise(PA_within=cor(observed_raw-ave(observed_raw,PopGroup),predicted_raw-ave(predicted_raw,PopGroup)),PA_across=cor(observed_raw,predicted_raw),.groups='drop') |> group_by(scheme,target,model,Environment) |> summarise(across(c(PA_within,PA_across),mean),.groups='drop')
cmp<-inner_join(metr(pr_old),metr(pr_new),by=c('scheme','target','model','Environment'),suffix=c('_old','_new')) |> mutate(d_within=PA_within_new-PA_within_old,d_across=PA_across_new-PA_across_old)
write.csv(cmp,file.path(here,'results/corrected_G',paste0(trait,'_old_vs_corrected_G.csv')),row.names=FALSE)
cat(trait,': max |change| within-group PA =',round(max(abs(cmp$d_within)),3),'; across-group PA =',round(max(abs(cmp$d_across)),3),'; mean change within =',round(mean(cmp$d_within),3),'\n')
