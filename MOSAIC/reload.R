require(mosaicpackage)
rseed=1
######################## first set some options ###############################
# important things to set
#shargs<-commandArgs(trailingOnly=TRUE) # read in arguments from the command line; 
#filename=shargs[1]  # e.g. simulated_2way_1-1_13-15_200_60_0.99_100.RData
#pathin=shargs[2]  # e.g. RESULTS/
#pathout=shargs[3]  # e.g. RESULTS/
#pathin="RESULTS/";pathout="RESULTS/";filename="simulated_2way_1-4_15-22_362_60_0.99_100.RData"
datasource="HGDP/";ANC=NULL
doMu=dotheta=dorho=doQ=T
load(paste0(pathin,filename));a.Mu=Mu;a.Q=Q;a.lambda=lambda;a.alpha=alpha;a.rho=rho;a.theta=theta;
load(paste0(pathin,"localanc_",filename));
if (target=="simulated") ANC=T
if ((sum(kLL)!=94) & (target!="simulated"))
{
  panels<-read.table(paste(datasource,"sample.names",sep=""), header=F);panels<-as.character(unique(panels[,1]))
  tmp=match(c(target,rownames(Mu)),panels);tmp=tmp[!is.na(tmp)]
  mask=panels[-tmp]
}
a.max.donors=max.donors;a.prop.don=prop.don
samp_chrnos=chrnos;subNUMA=NUMA;subNL=1000 # use them all
MC=NaN;verbose=T;USEHAPMIX=F;OUTHAPMIX=F;PHASE=T;FLAT=F;nl=100;prethin=F;EM=F
HPC=2 # whether to use ff() chromosome-by-chromosome (HPC=1) or chromosomeXind-by-chromsomeXind(HPC=2) or not at all (HPC=F);
ffpath="/dev/shm/" # location of fast-files
ffcleanup=T

if (target=="simulated") {o.lambda=20;set.seed(rseed)}
if (target!="simulated") o.lambda=mean(unlist(lambda))

require(doParallel)
if (is.na(MC)) {
  MC=as.integer(detectCores()/2)
  if (is.na(MC)) {MC=2;warning("using 2 cores as detectCores() has failed",immediate.=T)} # use 2 if can't use detectCores() 
}
if (verbose) cat("using", MC, "cores\n")
registerDoParallel(cores=MC)
#chrnos=21:22;nchrno=length(chrnos)
source("setup.R") # includes introducing phase errors if RPE>0
Mu=a.Mu;Q=a.Q;lambda=a.lambda;alpha=a.alpha;rho=a.rho;theta=a.theta; # reset to loaded versions as setup may have changed these
# need to re-create transitions and emissions 
mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,Q[[ind]],Mu,rho,NL)
if (exists("final.flips")) 
  flips=final.flips

