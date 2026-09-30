## Rebuild G with the marker coding rrBLUP::A.mat expects (-1/0/1).
## The stored data/markers/Gmatrix_125.rda was built from 0/1/2 codes, which made
## A.mat mis-estimate allele frequencies and silently drop 1,786 of 13,687 SNPs.
## Output is a NEW file; the original G is left untouched.
suppressMessages(library(rrBLUP))
a<-grep('^--file=',commandArgs(),value=TRUE);sf<-gsub('~+~',' ',sub('^--file=','',a[1]),fixed=TRUE)
project<-dirname(dirname(dirname(normalizePath(sf))))
e<-new.env();load(file.path(project,'data/markers/SNPs_125genotypes.rda'),envir=e);X<-e$SNPs
stopifnot(all(X%in%c(0,1,2,NA)))
imp<-A.mat(X-1,return.imputed=TRUE);G<-imp$A
stopifnot(ncol(imp$imputed)==ncol(X))
rownames(G)<-colnames(G)<-rownames(X)
ev<-eigen(G,symmetric=TRUE,only.values=TRUE)$values
if(min(ev)<=1e-8){cat('adding small ridge; min eigenvalue',min(ev),'\n');G<-G+diag(1e-6,nrow(G))}
old<-new.env();load(file.path(project,'data/markers/Gmatrix_125.rda'),envir=old)
stopifnot(identical(rownames(old$G),rownames(G)))
save(G,file=file.path(project,'data/markers/Gmatrix_125_corrected.rda'))
cat('SNPs used:',ncol(imp$imputed),'| cor off-diagonal with old G:',round(cor(G[upper.tri(G)],old$G[upper.tri(old$G)]),4),'\n')
