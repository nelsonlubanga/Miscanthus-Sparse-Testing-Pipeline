## Combined DM/FW/MC figures, traits as facets, trials labelled Aber and Brau.
## DM is read from the existing DM pipeline results, unchanged.
## Usage: Rscript 04_plot_all_traits.R <figure output folder>
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
pipeline<-dirname(dirname(normalizePath(sf)));project<-dirname(pipeline)
figdir<-commandArgs(trailingOnly=TRUE)[1];stopifnot(!is.na(figdir));dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
srcs<-c(DM=file.path(project,'reviewer_DM_pipeline/results'),FW=file.path(pipeline,'results/FW'),MC=file.path(pipeline,'results/MC'))
load_all<-function(f)bind_rows(lapply(names(srcs),function(t)data.frame(trait=t,read.csv(file.path(srcs[[t]],f)))))
tidy<-function(d){
 d$Environment<-factor(d$Environment,levels=c('ABR33','JKI'),labels=c('Aber','Brau'))
 d$trait<-factor(d$trait,levels=names(srcs));if('PopGroup'%in%names(d))d$PopGroup<-factor(d$PopGroup,levels=c('M. sacchariflorus','M. sinensis','M. x gig'),labels=c('M. sacchariflorus','M. sinensis','M. × giganteus'));d$model<-factor(d$model,levels=c('M1','M2','M3'))
 d$target<-factor(d$target,levels=c('Raw','Group_adjusted'),labels=c('Population group unadjusted','Population group adjusted'))
 d}
rm<-tidy(load_all('per_repeat_metrics.csv'));bg<-tidy(load_all('within_group_per_repeat.csv'))
cols<-c(M1='#F08A7E',M2='#4FAE62',M3='#5B9BD5')
box<-function(p)p+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+geom_boxplot(width=.7,position=position_dodge(.8),outlier.size=.5,show.legend=TRUE)+scale_fill_manual(values=cols,name='Model',drop=FALSE)+theme_bw(base_size=11)+theme(plot.title=element_text(size=12,hjust=.5))
save_both<-function(p,stem,w,h)for(ext in c('png','pdf'))ggsave(file.path(figdir,paste0(stem,'.',ext)),p,width=w,height=h,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else NULL)
for(cv in c('CV1','CV2')) {
 d<-rm[rm$scheme==cv & is.finite(rm$PA_within_pooled),]
 p<-box(ggplot(d,aes(Environment,PA_within_pooled,fill=model)))+facet_grid(target~trait)+
  labs(x=NULL,y='Predictive ability',title=paste('Pooled within-group predictive ability under',cv))
 save_both(p,paste0('ALL_TRAITS_',cv,'_within_group'),10,6)
 d<-bg[bg$scheme==cv & bg$target=='Population group adjusted' & is.finite(bg$PA),]
 p<-box(ggplot(d,aes(Environment,PA,fill=model)))+facet_grid(trait~PopGroup)+theme(strip.text.x=element_text(face='italic'))+
  labs(x=NULL,y='Predictive ability',title=paste('Predictive ability within each population group under',cv))
 save_both(p,paste0('ALL_TRAITS_',cv,'_species_specific'),10,8)
}
d<-rm[rm$scheme=='Sparse25' & rm$target=='Population group adjusted' & is.finite(rm$PA_within_pooled),];d$overlap<-factor(d$overlap,levels=seq(0,25,5))
p<-box(ggplot(d,aes(overlap,PA_within_pooled,fill=model)))+facet_grid(trait~Environment)+
 labs(x='Shared training genotypes (25 per trial)',y='Predictive ability',title='Pooled within-group predictive ability under sparse testing')
save_both(p,'ALL_TRAITS_sparse25_within_group',10,8)
cat('Saved combined DM/FW/MC figures in',figdir,'\n')
