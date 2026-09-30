## Restore site-year correlations and derive one pooled DM BLUE per trial.
## Default uses every available year (ABR33 2015/2016; JKI 2014/2015/2016).
## Correlations also supplied for the common 2015/2016 four-cell subset.
## Trial phenotypic model: fixed Year + Genotype, random Genotype:Year,
## persistent Plot, Year:Block, plus Year:Row and Year:Col at ABR33.
## Residual variances differ by Year. Predictions equally average fixed years.
## Genotype BLUEs span the same fitted fixed-effect space as the original
## nested fixed-population/genotype parameterization. They retain group means.
## This script addresses pooling years; it does not remove population structure.
suppressPackageStartupMessages({library(asreml);library(dplyr);library(tidyr);library(ggplot2)})
arg <- grep('^--file=',commandArgs(),value=TRUE)
sf <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if(is.null(sf)&&length(arg)) sf<-gsub('~+~',' ',sub('^--file=','',arg[1]),fixed=TRUE)
root<-if(!is.null(sf))dirname(dirname(normalizePath(sf))) else if(basename(getwd())=='scripts')dirname(getwd()) else getwd()
outdir<-file.path(root,'results','site_pooled_revision');dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
bl<-read.csv(file.path(root,'results/BLUEs/BLUEs_all_environments.csv'))
pg<-read.csv(file.path(root,'data/popgroups_125.csv'))
bl<-left_join(bl,pg,by='Genotype')
stopifnot(!anyNA(bl$PopGroup),!anyDuplicated(bl[c('Genotype','Environment','trait')]))
bl$Environment<-sub('^Aber-','ABR33-',sub('^Brau-','JKI-',bl$Environment))
lev<-c('ABR33-2015','ABR33-2016','JKI-2014','JKI-2015','JKI-2016')
bl<-bl |> group_by(trait,Environment,PopGroup) |> mutate(group_residual=BLUE-mean(BLUE)) |> ungroup()
cc<-list()
for(tr in unique(bl$trait)) for(a in seq_along(lev)) for(b in a:length(lev)) {
 x<-bl[bl$trait==tr & bl$Environment==lev[a],]; y<-bl[bl$trait==tr & bl$Environment==lev[b],]
 z<-merge(x,y,by='Genotype',suffixes=c('.x','.y'))
 cc[[length(cc)+1L]]<-data.frame(trait=tr,environment1=lev[a],environment2=lev[b],n=nrow(z),
 Pearson=cor(z$BLUE.x,z$BLUE.y),Spearman=cor(z$BLUE.x,z$BLUE.y,method='spearman'),
 group_adjusted_Pearson=cor(z$group_residual.x,z$group_residual.y))
}
cc<-bind_rows(cc);write.csv(cc,file.path(outdir,'BLUE_correlations_all_five.csv'),row.names=FALSE)
write.csv(cc[cc$environment1!='JKI-2014' & cc$environment2!='JKI-2014',],file.path(outdir,'BLUE_correlations_four_2015_2016.csv'),row.names=FALSE)
for(tr in unique(cc$trait)) {
 x<-cc[cc$trait==tr,];mat<-matrix(NA_real_,5,5,dimnames=list(lev,lev));nn<-mat
 for(i in seq_len(nrow(x))) {j<-match(x$environment1[i],lev);k<-match(x$environment2[i],lev);mat[j,k]<-mat[k,j]<-x$Pearson[i];nn[j,k]<-nn[k,j]<-x$n[i]}
 write.csv(mat,file.path(outdir,paste0('BLUE_correlation_matrix_',tr,'.csv')))
 write.csv(nn,file.path(outdir,paste0('BLUE_pairwise_sample_sizes_',tr,'.csv')))
}
plotdata<-bind_rows(cc,cc[cc$environment1!=cc$environment2,] |> mutate(tmp=environment1,environment1=environment2,environment2=tmp) |> select(-tmp))
plotdata$environment1<-factor(plotdata$environment1,levels=lev);plotdata$environment2<-factor(plotdata$environment2,levels=rev(lev))
p<-ggplot(plotdata,aes(environment1,environment2,fill=Pearson))+geom_tile(colour='white')+
 geom_text(aes(label=sprintf('%.2f\n(n=%d)',Pearson,n)),size=2.7)+facet_wrap(~trait,nrow=1)+
 scale_fill_gradient2(low='#B2182B',mid='white',high='#2166AC',midpoint=0,limits=c(-1,1),name='Pearson r')+
 labs(x=NULL,y=NULL)+theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
for(ext in c('png','pdf'))ggsave(file.path(outdir,paste0('BLUE_correlations.',ext)),p,width=12,height=4,dpi=300,bg='white')
d<-read.csv(file.path(root,'data/combined_phenotypes_long.csv'));d<-d[d$trait=='DM',]
d$Trial<-ifelse(grepl('^Aber-',d$Environment),'ABR33','JKI')
blue_out<-list();vcs<-list();fitchecks<-list()
for(site in c('ABR33','JKI')) {
 sub<-d[d$Trial==site & is.finite(d$value),]
 stopifnot(!anyNA(sub[c('Genotype','Year','Block','Row','Col')]))
 ## Coordinates checked for uniqueness within Year; persistent plot identity
 ## uses physical position, avoiding treating repeated years as new plots.
 sub$Plot<-interaction(sub$Block,sub$Row,sub$Col,drop=TRUE)
 stopifnot(!anyDuplicated(sub[c('Year','Plot')]))
 map<-unique(sub[c('Plot','Genotype')]);stopifnot(!anyDuplicated(map$Plot))
 for(n in c('Genotype','Year','Block','Row','Col','Plot'))sub[[n]]<-factor(sub[[n]])
 sub<-sub[order(sub$Year,sub$Plot),]
 random_form<-if(site=='ABR33')~Genotype:Year+Plot+Year:Block+Year:Row+Year:Col else ~Genotype:Year+Plot+Year:Block
 fit<-asreml(fixed=value~Year+Genotype,random=random_form,residual=~dsum(~units|Year),data=sub,maxit=60,trace=FALSE)
 ## Force an update to check stability, even if the first convergence flag is TRUE.
 for(i in 1:5) {
  fit<-update(fit,maxit=60)
  vc<-summary(fit)$varcomp
  change_col<-grep('ch',colnames(vc),ignore.case=TRUE,value=TRUE)
  stable<-!length(change_col)||all(abs(vc[[change_col[1]]][vc$bound=='P'])<1,na.rm=TRUE)
  if(isTRUE(fit$converge)&&stable)break
 }
 if(!isTRUE(fit$converge))stop('Pooled model did not converge: ',site)
 av<-list(Year=rep(1/nlevels(sub$Year),nlevels(sub$Year)))
 pv<-predict(fit,classify='Genotype',average=av)$pvals
 stopifnot(all(pv$status=='Estimable'),all(is.finite(pv$predicted.value)))
 blue_out[[site]]<-data.frame(Genotype=as.character(pv$Genotype),Environment=site,trait='DM',BLUE=pv$predicted.value,SE=pv$std.error,converged=fit$converge)
 vcs[[site]]<-data.frame(Trial=site,term=rownames(vc),vc,row.names=NULL)
 fitchecks[[site]]<-data.frame(Trial=site,years=paste(levels(sub$Year),collapse=','),n_records=nrow(sub),n_plots=nlevels(sub$Plot),n_genotypes=nlevels(sub$Genotype),converged=fit$converge,stable=stable)
 saveRDS(fit,file.path(outdir,paste0('pooled_fit_',site,'.rds')))
 cat(site,'pooled BLUEs:',nrow(pv),'stable:',stable,'\n')
}
pooled<-bind_rows(blue_out)
write.csv(pooled,file.path(outdir,'pooled_DM_BLUEs.csv'),row.names=FALSE)
write.csv(bind_rows(vcs),file.path(outdir,'pooled_variance_components.csv'),row.names=FALSE)
write.csv(bind_rows(fitchecks),file.path(outdir,'pooled_model_checks.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'phenotypic_sessionInfo.txt'))
