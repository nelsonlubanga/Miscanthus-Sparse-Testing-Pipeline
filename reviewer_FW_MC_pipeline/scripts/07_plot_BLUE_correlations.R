## Correlation matrix of site-year BLUEs for DM, FW and MC, redrawn for legibility
## (large numbers, white text on dark tiles), as requested by Gancho.
## Input: results/site_pooled_revision/BLUE_correlations_all_five.csv (unchanged).
## Usage: Rscript 07_plot_BLUE_correlations.R <figure output folder>
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
project<-dirname(dirname(dirname(normalizePath(sf))))
figdir<-commandArgs(trailingOnly=TRUE)[1];stopifnot(!is.na(figdir))
x<-read.csv(file.path(project,'results/site_pooled_revision/BLUE_correlations_all_five.csv'))
lab<-function(e)sub('ABR33-20','Aber',sub('JKI-20','Brau',e))
lev<-c('Aber15','Aber16','Brau14','Brau15','Brau16')
x$e1<-factor(lab(x$environment1),levels=lev);x$e2<-factor(lab(x$environment2),levels=lev)
d<-x[x$environment1!=x$environment2,]
## lower triangle: row = later site-year, column = earlier
d<-transform(d,row=factor(e2,levels=rev(lev)),col=factor(e1,levels=lev))
d$label<-sprintf('%.2f',d$Pearson)
p<-ggplot(d,aes(col,row,fill=Pearson))+geom_tile(colour='white',linewidth=1)+
 geom_text(aes(label=label,colour=Pearson>=0.6),size=6.5,fontface='bold',show.legend=FALSE)+
 scale_colour_manual(values=c(`TRUE`='white',`FALSE`='black'))+
 scale_fill_gradientn(colours=c('#F7FBFF','#DEEBF7','#9ECAE1','#4292C6','#08519C','#08306B'),values=c(0,.2,.45,.62,.8,1),limits=c(0,1),name='Pearson r')+
 facet_wrap(~trait,nrow=1)+labs(x=NULL,y=NULL)+coord_equal()+
 theme_minimal(base_size=14)+theme(panel.grid=element_blank(),axis.text=element_text(size=12,colour='black'),axis.text.x=element_text(angle=45,hjust=1),strip.text=element_text(size=15,face='bold'),legend.title=element_text(size=12))
for(ext in c('png','pdf'))ggsave(file.path(figdir,paste0('BLUE_correlations_site_years.',ext)),p,width=15,height=5.5,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else NULL)
cat('Saved legible site-year BLUE correlation matrix\n')
