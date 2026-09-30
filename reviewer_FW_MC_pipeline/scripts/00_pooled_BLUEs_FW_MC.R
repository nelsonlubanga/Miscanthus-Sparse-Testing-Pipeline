## One pooled BLUE per genotype per physical trial for FW and MC.
## Identical phenotypic model to the DM pooled BLUEs (../../scripts/site_pooled_blue_revision.R):
## fixed Year + Genotype, random Genotype:Year + persistent Plot + Year:Block,
## plus Year:Row and Year:Col at ABR33; year-specific residual variance.
## Predictions equally average fixed years. BLUEs retain group means.
suppressPackageStartupMessages({library(asreml);library(dplyr)})
a<-grep('^--file=',commandArgs(),value=TRUE)
sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
pipeline<-dirname(dirname(normalizePath(sf)));project<-dirname(pipeline)
outdir<-file.path(pipeline,'results','pooled_BLUEs');dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d0<-read.csv(file.path(project,'data/combined_phenotypes_long.csv'))
d0$Trial<-ifelse(grepl('^Aber-',d0$Environment),'ABR33','JKI')
for(trait in c('FW','MC')) {
 d<-d0[d0$trait==trait,]
 blue_out<-list();vcs<-list();fitchecks<-list()
 for(site in c('ABR33','JKI')) {
  sub<-d[d$Trial==site & is.finite(d$value),]
  stopifnot(!anyNA(sub[c('Genotype','Year','Block','Row','Col')]))
  sub$Plot<-interaction(sub$Block,sub$Row,sub$Col,drop=TRUE)
  stopifnot(!anyDuplicated(sub[c('Year','Plot')]))
  map<-unique(sub[c('Plot','Genotype')]);stopifnot(!anyDuplicated(map$Plot))
  for(n in c('Genotype','Year','Block','Row','Col','Plot'))sub[[n]]<-factor(sub[[n]])
  sub<-sub[order(sub$Year,sub$Plot),]
  random_form<-if(site=='ABR33')~Genotype:Year+Plot+Year:Block+Year:Row+Year:Col else ~Genotype:Year+Plot+Year:Block
  fit<-asreml(fixed=value~Year+Genotype,random=random_form,residual=~dsum(~units|Year),data=sub,maxit=60,trace=FALSE)
  for(i in 1:5) {
   fit<-update(fit,maxit=60)
   vc<-summary(fit)$varcomp
   change_col<-grep('ch',colnames(vc),ignore.case=TRUE,value=TRUE)
   stable<-!length(change_col)||all(abs(vc[[change_col[1]]][vc$bound=='P'])<1,na.rm=TRUE)
   if(isTRUE(fit$converge)&&stable)break
  }
  if(!isTRUE(fit$converge))stop('Pooled model did not converge: ',trait,' ',site)
  av<-list(Year=rep(1/nlevels(sub$Year),nlevels(sub$Year)))
  pv<-predict(fit,classify='Genotype',average=av)$pvals
  stopifnot(all(pv$status=='Estimable'),all(is.finite(pv$predicted.value)))
  blue_out[[site]]<-data.frame(Genotype=as.character(pv$Genotype),Environment=site,trait=trait,BLUE=pv$predicted.value,SE=pv$std.error,converged=fit$converge)
  vcs[[site]]<-data.frame(trait=trait,Trial=site,term=rownames(vc),vc,row.names=NULL)
  fitchecks[[site]]<-data.frame(trait=trait,Trial=site,years=paste(levels(sub$Year),collapse=','),n_records=nrow(sub),n_plots=nlevels(sub$Plot),n_genotypes=nlevels(sub$Genotype),converged=fit$converge,stable=stable)
  saveRDS(fit,file.path(outdir,paste0('pooled_fit_',trait,'_',site,'.rds')))
  cat(trait,site,'pooled BLUEs:',nrow(pv),'converged:',fit$converge,'stable:',stable,'\n')
 }
 write.csv(bind_rows(blue_out),file.path(outdir,paste0('pooled_',trait,'_BLUEs.csv')),row.names=FALSE)
 write.csv(bind_rows(vcs),file.path(outdir,paste0('pooled_variance_components_',trait,'.csv')),row.names=FALSE)
 write.csv(bind_rows(fitchecks),file.path(outdir,paste0('pooled_model_checks_',trait,'.csv')),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(outdir,'phenotypic_sessionInfo.txt'))
