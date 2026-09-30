## CV1/CV2 within-group figures in the original DM layout (Group_adjusted | Raw),
## with FW and MC added in the same row. DM read from the DM pipeline, unchanged.
## Usage: Rscript 05_plot_within_group_one_row.R <figure output folder>
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
pipeline<-dirname(dirname(normalizePath(sf)));project<-dirname(pipeline)
figdir<-commandArgs(trailingOnly=TRUE)[1];stopifnot(!is.na(figdir));dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
srcs<-c(DM=file.path(project,'reviewer_DM_pipeline/results'),FW=file.path(pipeline,'results/FW'),MC=file.path(pipeline,'results/MC'))
rm<-bind_rows(lapply(names(srcs),function(t)data.frame(trait=t,read.csv(file.path(srcs[[t]],'per_repeat_metrics.csv')))))
rm$Environment<-factor(rm$Environment,levels=c('ABR33','JKI'),labels=c('Aber','Brau'))
rm$trait<-factor(rm$trait,levels=names(srcs));rm$model<-factor(rm$model,levels=c('M1','M2','M3'));rm$target<-factor(rm$target,levels=c('Raw','Group_adjusted'),labels=c('Population group unadjusted','Population group adjusted'))
cols<-c(M1='#F08A7E',M2='#4FAE62',M3='#5B9BD5')
for(cv in c('CV1','CV2')) {
 d<-rm[rm$scheme==cv & is.finite(rm$PA_within_pooled),]
 p<-ggplot(d,aes(Environment,PA_within_pooled,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.7,show.legend=TRUE)+facet_grid(.~trait+target)+scale_fill_manual(values=cols,name='Model',drop=FALSE)+labs(x=NULL,y='Predictive ability',title=paste('Pooled within-group predictive ability under',cv))+theme_bw(base_size=11)+theme(plot.title=element_text(size=12,hjust=.5))
 for(ext in c('png','pdf'))ggsave(file.path(figdir,paste0(cv,'_within_group_DM_FW_MC.',ext)),p,width=16,height=4.5,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else NULL)
}
cat('Saved one-row CV1/CV2 within-group figures in',figdir,'\n')
