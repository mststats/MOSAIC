#rm(list=ls());
writelog<-function(alg) # single consistent function to write to EMlogfile
  write(file=EMlogfile,c(alg,signif(diff.time,4),signif(t(Mu),4),signif(rho,4),c(sapply(Q, function(x) signif(t(x),4))),
			 sapply(alpha, function(x) signif(x,4)),sapply(lambda,function(x) round(x,4)),signif(theta,4),round(cloglike,4)),ncol=length(lognames),append=T)
NUMA=2;ANC=T
require(doParallel)
MC=as.integer(detectCores()/2);if (is.na(MC)) {MC=2;warning("using 2 cores as detectCores() has failed",immediate.=T)} # use 2 if can't use detectCores() 
registerDoParallel(cores=MC)
HPC=2;ffpath="/dev/shm/" #cores and whether to use ff()
subG<-1;subNL=100
initonly=F;
FLAT=F
optlevel=3
USEHAPMIX=F
PHASE=F;RPE=0
o.lambda=1000
get_switches=F
OUTHAPMIX=F
tol=1e-8
verbose=T; 
shargs<-commandArgs(trailingOnly=TRUE)
datasource="chromopainter/"
panels=strsplit(shargs[1], ",")[[1]]
#panels=c("French","Yoruba")
chrnos=22:1;
GpcM=60;nl<-100;
prop.don=0.99
NUMA=2
L=1
LOG=T
total=200
nchrno=length(chrnos)
EM=T
dotheta=T
doMu=F
dorho=T
doQ=F
PLOT=F

chunkcounts=NULL
load("panels/hgdp_95_panels.22.RData") # should find a faster way...
lpanels=rep(NaN,length(panels)); for (np in 1:length(panels)) lpanels[np]=nrow(multipanels[[which(names(multipanels)==panels[np])]])/2
for (np in 1:length(panels))
{
  target=panels[np]
  for (firstind in 1:lpanels[np])
  {
    cat("Looking at individual",firstind,"in", target,"panel\n")
    source("setup.R") # includes introducing phase errors if RPE>0
    source("noanc.R") # always need to run noanc.R b/c need good paras for init_Mu
    switches=list();
    for (ch in 1:nchrno)
    {
      switches[[ch]]=list()
      tmp<-foreach(ind=1:NUMI) %dopar%
      create_donates(T,ch,ind,umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],flips[[ind]][[ch]],kLL,Mu,rho,theta,HPC) 
      for (ind in 1:NUMI)
      {
	if (NUMA==1)
	  switches[[ch]][[ind]]<-tmp[[ind]]$switches[[1]] 
	if (NUMA==2)
	  switches[[ch]][[ind]]<-tmp[[ind]]$switches[[1]] +tmp[[ind]]$switches[[2]] 
      }
      rm(tmp)
    } 
    tmpchunkcounts=rowSums(sapply(switches, function(x) rowSums(sapply(x,colSums))))/2
    tmpchunkcounts=tmpchunkcounts[seq(1,NUMP,2)]+tmpchunkcounts[seq(2,NUMP,2)] # sum over both haps in each donor ind
    chunkcounts=rbind(chunkcounts,tmpchunkcounts)
    #boxplot(tmpchunkcounts~label[1:NUMP],names=pops$panels)
  }
}
nind=nrow(chunkcounts)
tmp=matrix(0,nind,nind);for (i in 1:nind) tmp[i,-i]=chunkcounts[i,]
chunkcounts=tmp
tmp=rep(NaN,sum(lpanels));l=0;for (np in 1:length(panels)) for (i in 1:lpanels[np]) {l=l+1;tmp[l]=paste(panels[np],i,sep=":")}
rownames(chunkcounts)=colnames(chunkcounts)=tmp
#boxplot(chunkcounts[i,-i]~tapply((1:nind)[-i],(1:nind)[-i], function(j) strsplit(colnames(chunkcounts)[j],":")[[1]][1]),main=colnames(chunkcounts)[i])
#image(1:nind,1:nind,chunkcounts)
save(chunkcounts,file=paste0(paste(panels,collapse="_"),"_chromochunks.RData"))
