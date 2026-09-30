## Figures for FW and MC in the same layout as the DM figures sent to Gancho
## (reply_to_gancho_v2/figures_Aber_Brau), with trials labelled Aber and Brau.
## Also writes a DM/FW/MC summary table of the main M3 comparison.
## Usage: Rscript 03_plot_figures.R <figure output folder>
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
pipeline<-dirname(dirname(normalizePath(sf)));project<-dirname(pipeline)
figdir<-commandArgs(trailingOnly=TRUE)[1];stopifnot(!is.na(figdir));dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
relabel<-function(d){d$Environment<-factor(d$Environment,levels=c('ABR33','JKI'),labels=c('Aber','Brau'));if('model'%in%names(d))d$model<-factor(d$model,levels=c('M1','M2','M3'));if('PopGroup'%in%names(d))d$PopGroup<-factor(d$PopGroup,levels=c('M. sacchariflorus','M. sinensis','M. x gig'),labels=c('M. sacchariflorus','M. sinensis','M. × giganteus'));d}
cols<-c(M1='#F08A7E',M2='#4FAE62',M3='#5B9BD5')
save_both<-function(p,stem,w,h)for(ext in c('png','pdf'))ggsave(file.path(figdir,paste0(stem,'.',ext)),p,width=w,height=h,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else NULL)
for(trait in c('FW','MC')) {
 res<-file.path(pipeline,'results',trait)
 rm<-relabel(read.csv(file.path(res,'per_repeat_metrics.csv')));rm$target_label<-factor(rm$target,levels=c('Raw','Group_adjusted'),labels=c('Population group unadjusted','Population group adjusted'));bg<-relabel(read.csv(file.path(res,'within_group_per_repeat.csv')))
 for(cv in c('CV1','CV2')) {
  d<-rm[rm$scheme==cv & is.finite(rm$PA_within_pooled),]
  p<-ggplot(d,aes(Environment,PA_within_pooled,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.7,show.legend=TRUE)+facet_wrap(~target_label,nrow=1)+scale_fill_manual(values=cols,name='Model',drop=FALSE)+labs(x=NULL,y='Predictive ability',title=paste('Pooled within-group predictive ability for',trait,'under',cv))+theme_bw(base_size=11)+theme(plot.title=element_text(size=12,hjust=.5))
  save_both(p,paste0(trait,'_',cv,'_within_group'),9,4)
 }
 d<-bg[bg$scheme%in%c('CV1','CV2') & bg$target=='Group_adjusted' & is.finite(bg$PA),]
 p<-ggplot(d,aes(Environment,PA,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.5,show.legend=TRUE)+facet_grid(scheme~PopGroup)+scale_fill_manual(values=cols,name='Model',drop=FALSE)+labs(x=NULL,y='Predictive ability',title=paste('Predictive ability within each population group for',trait))+theme_bw(base_size=10)+theme(plot.title=element_text(size=11,hjust=.5),strip.text.x=element_text(face='italic'))
 save_both(p,paste0(trait,'_species_specific_prediction'),11,6)
 d<-rm[rm$scheme=='Sparse25'&rm$target=='Group_adjusted',];d$overlap<-factor(d$overlap,levels=seq(0,25,5))
 p<-ggplot(d,aes(overlap,PA_within_pooled,fill=model))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.5,show.legend=TRUE)+facet_wrap(~Environment,nrow=1)+scale_fill_manual(values=cols,name='Model',drop=FALSE)+labs(x='Shared training genotypes (25 per trial)',y='Predictive ability',title=paste('Pooled within-group predictive ability for',trait,'under sparse testing'))+theme_bw(base_size=11)+theme(plot.title=element_text(size=12,hjust=.5))
 save_both(p,paste0(trait,'_sparse25_within_group'),9,4)
}
## Summary across traits: DM is read from the existing DM pipeline, unchanged.
srcs<-c(DM=file.path(project,'reviewer_DM_pipeline/results'),FW=file.path(pipeline,'results/FW'),MC=file.path(pipeline,'results/MC'))
s<-bind_rows(lapply(names(srcs),function(t)data.frame(trait=t,read.csv(file.path(srcs[[t]],'prediction_summary.csv')))))
r2<-bind_rows(lapply(names(srcs),function(t)data.frame(trait=t,read.csv(file.path(srcs[[t]],'descriptive_population_R2.csv')))))
s<-relabel(s);r2<-relabel(r2)
main<-s |> filter(scheme%in%c('CV1','CV2'),model=='M3') |> group_by(trait,scheme,Environment) |> summarise(
 overall_PA_unadjusted=PA_raw[target=='Raw'],within_group_PA_unadjusted=PA_within_pooled[target=='Raw'],within_group_PA_group_adjusted=PA_within_pooled[target=='Group_adjusted'],
 RMSE_M3_adjusted=RMSE_raw[target=='Group_adjusted'],RMSE_group_only=RMSE_group_only[target=='Group_adjusted'],.groups='drop')
write.csv(main,file.path(figdir,'summary_M3_DM_FW_MC.csv'),row.names=FALSE)
models<-s |> filter(target=='Group_adjusted',scheme%in%c('CV1','CV2')) |> select(trait,scheme,Environment,model,PA_within_pooled) |> tidyr::pivot_wider(names_from=model,values_from=PA_within_pooled)
write.csv(models,file.path(figdir,'summary_adjusted_within_group_M1_M3.csv'),row.names=FALSE)
sp<-s |> filter(scheme=='Sparse25',target=='Group_adjusted') |> select(trait,Environment,overlap,model,PA_within_pooled) |> tidyr::pivot_wider(names_from=model,values_from=PA_within_pooled)
write.csv(sp,file.path(figdir,'summary_sparse25_adjusted_within_group.csv'),row.names=FALSE)
write.csv(r2,file.path(figdir,'summary_group_R2.csv'),row.names=FALSE)
## Comparison figure: M3 overall vs within-group (adjusted), all three traits.
d<-main |> tidyr::pivot_longer(c(overall_PA_unadjusted,within_group_PA_unadjusted,within_group_PA_group_adjusted),names_to='metric',values_to='r')
d$metric<-factor(d$metric,levels=c('overall_PA_unadjusted','within_group_PA_unadjusted','within_group_PA_group_adjusted'),labels=c('Overall','Within population groups, unadjusted','Within population groups, adjusted'))
p<-ggplot(d,aes(Environment,r,fill=metric))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_col(width=.75,position=position_dodge(.8))+geom_text(aes(label=sprintf('%.2f',r)),position=position_dodge(.8),vjust=ifelse(d$r>=0,-.3,1.2),size=2.6)+facet_grid(scheme~trait)+scale_fill_manual(values=c('#9E9E9E','#8FB8DE','#2E6DA4'),name=NULL)+labs(x=NULL,y='Predictive ability (M3, mean of 10 repetitions)')+theme_bw(base_size=11)+theme(plot.title=element_text(size=12,hjust=.5))+theme(legend.position='bottom')
save_both(p,'M3_overall_vs_within_group_DM_FW_MC',10,6)
print(as.data.frame(main),digits=3);print(as.data.frame(models),digits=3);print(as.data.frame(r2),digits=3)
