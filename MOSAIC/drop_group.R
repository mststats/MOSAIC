# function to drop groups, indexed by rem.group
if (!exists("updateparas")) updateparas=F
Mu<-Mu[-rem.group,];Mu<-t(t(Mu)/colSums(Mu)) # rescale
panels<-panels[-rem.group]
keep<-!(label[KNOWN]%in%rem.group)
for (ch in 1:nchrno)
{
  f=function(g) {ans=d.w[[ch]]$w[[g]][keep];ans=ans[!is.na(ans)];ans}
  d.w[[ch]]$w<-lapply(1:G[ch], f)
}
mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
KNOWN<-KNOWN[!(label%in%rem.group)]
full.NUMP=NUMP
NUMP<-NUMP-sum(NL[rem.group]);NN<-NN-sum(NL[rem.group]);NL<-NL[-rem.group]
max.donors<-min(full.max.donors,NUMP) # keep as high as possible
tmp.d=-log(1-rho)/full.NUMP;rho<-1-exp(-tmp.d*NUMP)# this will decrease; base prediction on effective genetic distance estimate
tmp.phi=theta*full.NUMP/(1-theta);theta=tmp.phi/(tmp.phi+NUMP)# this will increase; base prediction on effective error estimate
label<-label[!(label%in%rem.group)]
# now relabel all to go from 1 to LL
relabel<-label
for (k in 1:kLL)
  relabel[label==unique(label)[k]]<-k
label<-relabel
LL<-LL-length(rem.group);kLL<-kLL-length(rem.group)
for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,Q[[ind]],Mu,rho,NL)
source("initProb.R")
LOG=F
if (prethin)
  source("pre_all_donates.R")
source("all_donates.R")
if (updateparas)
{
  o.dotheta=dotheta;o.dorho=dorho;o.doQ=doQ;o.doMu=doMu
  doMu=T;dotheta=F;dorho=F;doQ=F
  cloglike=NaN;total=50;source("mosaic.R") # could try increasing rho and there first as they will be expected to increase with fewer donors available
  dotheta=o.dotheta;dorho=o.dorho;doQ=o.doQ;doMu=o.doMu
  source("all_donates.R") 
}
LOG=o.LOG
