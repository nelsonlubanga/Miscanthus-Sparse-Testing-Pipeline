## Impact of the corrected G on the ORIGINAL site-year CV1/CV2 analysis
## (scripts/07_fit_cv1_cv2_models.R): same data, PopGroup fixed, same fold
## assignments and BGLR settings. Each task is fitted twice, with the original
## and the corrected G, using an identical MCMC seed, so differences reflect G only.
## The original script sets no MCMC seed, so its published values are not reused.
## Usage: Rscript 11_corrected_G_siteyear_check.R [trait=DM]
suppressMessages({library(dplyr);library(BGLR);library(parallel)})
TRAIT<-commandArgs(trailingOnly=TRUE)[1];if(is.na(TRAIT))TRAIT<-'DM'
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
here<-dirname(dirname(normalizePath(sf)));project<-dirname(here)
N_REPS<-10;N_FOLDS<-5;ENVS<-c('Aber-2015','Aber-2016','Brau-2014','Brau-2015','Brau-2016')
blues<-read.csv(file.path(project,'results/BLUEs/BLUEs_all_environments.csv'),stringsAsFactors=FALSE)
popgroups<-read.csv(file.path(project,'data/popgroups_125.csv'),stringsAsFactors=FALSE)
gload<-function(f){e<-new.env();load(file.path(project,'data/markers',f),envir=e);e$G}
Gs<-list(original=gload('Gmatrix_125.rda'),corrected=gload('Gmatrix_125_corrected.rda'))
stopifnot(identical(rownames(Gs$original),rownames(Gs$corrected)))
GENOTYPES<-rownames(Gs$original)
obs<-expand.grid(Genotype=GENOTYPES,Environment=ENVS,stringsAsFactors=FALSE)
obs$PopGroup<-setNames(popgroups$PopGroup,popgroups$Genotype)[obs$Genotype]
Zg<-model.matrix(~factor(obs$Genotype,levels=GENOTYPES)-1);Ze<-model.matrix(~factor(obs$Environment,levels=ENVS)-1)
Xpop<-model.matrix(~factor(obs$PopGroup))[,-1,drop=FALSE]
K<-lapply(Gs,function(G){Ko<-Zg%*%G%*%t(Zg);list(obs=Ko,gxe=Ko*(Ze%*%t(Ze)))})
sub<-blues[blues$trait==TRAIT,c('Genotype','Environment','BLUE')]
y_full<-left_join(obs[c('Genotype','Environment')],sub,by=c('Genotype','Environment'))$BLUE
## fold assignment copied from scripts/07_fit_cv1_cv2_models.R
cv1_mask<-function(r,k){set.seed(1000*r+nchar(TRAIT));g<-sample(GENOTYPES);f<-setNames(rep(1:N_FOLDS,length.out=length(g)),g);obs$Genotype%in%names(f)[f==k]}
cv2_mask<-function(r,k){set.seed(2000*r+nchar(TRAIT));ri<-which(!is.na(y_full));p<-sample(ri);f<-setNames(rep(1:N_FOLDS,length.out=length(p)),p);m<-rep(FALSE,length(y_full));m[as.integer(names(f)[f==k])]<-TRUE;m}
fit<-function(y,model,Kk,seed){ETA<-list(POP=list(X=Xpop,model='FIXED'),ENV=list(X=Ze,model='BRR'),LINE=list(X=Zg,model='BRR'),GENO=list(K=Kk$obs,model='RKHS'))
 if(model=='M3')ETA$GXE<-list(K=Kk$gxe,model='RKHS')
 d<-tempfile('bg_');dir.create(d);on.exit(unlink(d,recursive=TRUE));set.seed(seed)
 BGLR(y=y,ETA=ETA,nIter=1500,burnIn=300,thin=3,verbose=FALSE,saveAt=file.path(d,'bg_'))$yHat}
tasks<-expand.grid(model=c('M2','M3'),CV=c('CV1','CV2'),rep=1:N_REPS,fold=1:N_FOLDS,stringsAsFactors=FALSE)
res<-mclapply(seq_len(nrow(tasks)),function(i){t<-tasks[i,];m<-if(t$CV=='CV1')cv1_mask(t$rep,t$fold) else cv2_mask(t$rep,t$fold)
 y<-ifelse(m,NA_real_,y_full);seed<-7000000L+i
 bind_rows(lapply(names(K),function(g){yh<-fit(y,t$model,K[[g]],seed)
  data.frame(G=g,t,Environment=ENVS,PA=sapply(ENVS,function(e){ix<-m&obs$Environment==e&!is.na(y_full);if(sum(ix)<5)NA_real_ else cor(yh[ix],y_full[ix])}))}))},mc.cores=10L)
out<-bind_rows(res)
dir.create(file.path(here,'results/corrected_G'),showWarnings=FALSE)
write.csv(out,file.path(here,'results/corrected_G',paste0(TRAIT,'_siteyear_fold_PA_original_vs_corrected_G.csv')),row.names=FALSE)
s<-out|>group_by(G,model,CV,Environment)|>summarise(PA=mean(PA,na.rm=TRUE),.groups='drop')|>tidyr::pivot_wider(names_from=G,values_from=PA)|>mutate(diff=corrected-original)
write.csv(s,file.path(here,'results/corrected_G',paste0(TRAIT,'_siteyear_mean_PA_original_vs_corrected_G.csv')),row.names=FALSE)
print(as.data.frame(s),digits=3);cat(TRAIT,'site-year: max |change| in mean PA =',round(max(abs(s$diff)),3),'\n')
