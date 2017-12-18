shargs<-commandArgs(trailingOnly=TRUE) # read in arguments from the command line; 
filename=shargs[1]
resultsdir=pathin="RESULTS/"
GpcM=60;chrnos=1:22
L=3
min.Mu.ratio=0;o.drop.threshold=-0.005;o.M=M=0

s.total=10;o.total=100;REPS=3
require(doParallel)
MC=16;registerDoParallel(cores=MC)
cat("Looking at ", filename, "\n")
target=strsplit(filename,"_")[[1]][1]
L=as.integer(strsplit(strsplit(filename,"_")[[1]][2],"way")[[1]][1])
firstind=as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]][1])
NUMI=diff(as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]]))+1
tmp=as.integer(strsplit(strsplit(filename,"_")[[1]][4],"-")[[1]]);chrnos=tmp[1]:tmp[2]
NN=as.integer(strsplit(filename,"_")[[1]][5])
GpcM=as.integer(strsplit(filename,"_")[[1]][6])
prop.don=as.numeric(strsplit(filename,"_")[[1]][7])
max.donors=as.integer(strsplit(strsplit(filename,"_")[[1]][8],"R")[[1]][1])
ANC=NULL
source("reload.R")
get_switches=F;LOG=T;PLOT=F;eps.lower=log(2);min.bg=0.1;max.bg=1.0;require(boot);
cloglike=NaN
source("initProb.R")
EM=TRUE
source("intermediate_calcs.R")
# run additional rounds of thin/phase/EM
total<-s.total # number of EM steps with repeated loop
eps=log(1.01) # i.e. a 1% increase in relative likelihood
old.runtime<-as.numeric(Sys.time())
writelog<-function(alg) # single consistent function to write to EMlogfile
  write(file=EMlogfile,c(alg,signif(diff.time,4),signif(t(Mu),4),signif(rho,4),c(sapply(Q, function(x) signif(t(x),4))),
			 sapply(alpha, function(x) signif(x,4)),sapply(lambda,function(x) round(x,4)),signif(theta,4),round(cloglike,4)),ncol=length(lognames),append=T)
source("create_logfile.R")
source("all_donates.R")
if (EM) 
{
  for (reps in 1:REPS) 
  {
    cat("######################## round ", reps, "of ",  REPS, "#######################\n")
    if (reps==(REPS-1)) # drop on second last pass so that EM and phasing are done once afterwards 
    {
      drop.threshold=o.drop.threshold
      #save.image(paste0(resultsdir,"drop_",target,"_",L,"_",firstind,"_",NUMA,".RData"));stop("check group dropping performance")
    } 
    if (reps!=(REPS-1)) drop.threshold=NaN # don't drop on any other pass 
    if (reps==REPS) 
      M=o.M # on last rep do more MCMC phasing
    source("mosaic.R") # location of this an issue. If above thin&phase, low lambda. If below then first rep has lowered log-like
    old.kLL=kLL
    if (min.Mu.ratio>0) # removal of unused groups
      source("rem_groups.R")
    # if we've dropped any groups then we don't need to re-run thinning as drop_group.R has done this already; otherwise call all_donates.R
    if (max.donors<NUMP & (kLL==old.kLL))
      source("all_donates.R") # decide on donor set using updated parameters
    if (PHASE) # & reps>1) # phasing early not the cause of low lambda. 
    {
      #stop(1)
      source("phase_hunt.R") 
      if (M>0) 
	source("phase_mcmc.R") # now do MCMC re-phasing w/o output to console (ugly due to txtProgressBar)
    }
  }
  total=o.total # do longer run EM on last rep
  if (verbose)
    cat("run one final round of EM\n")
  source("mosaic.R") 
}
final.flips=flips
source("ancunaware.R") # functions for getting localanc and gfbs that are ancestry unaware
EM=F;getnoancgfbs=T;eps=log(1.01);o.LOG=F;PLOT=F;get_switches=F;
a.Mu=Mu;a.rho=rho;a.theta=theta;a.Q=Q;a.alpha=alpha;a.lambda=lambda
a.o.Mu=o.Mu;a.o.rho=o.rho;a.o.theta=o.theta;a.o.Q=o.Q;a.o.alpha=o.alpha;a.o.lambda=o.lambda
samp_chrnos=chrnos;subNUMA=NUMA;subNL=max(NL) # use them all
source("coancestry.R")

######### fully Mosaic curves with Mosaic phasing ############
gfbs<-get_gfbs();source("localanc.R") 
if (verbose) cat("saving localanc results to file\n")
if (target!="simulated")
  save(file=paste0(resultsdir,"localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), 
       localanc, final.flips, g.loc)
if (target=="simulated")
  save(file=paste0(resultsdir,"localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), 
       localanc, g.true_anc, final.flips, g.loc)
save(file=paste0(resultsdir,"gfbs_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), gfbs)
source("plot_funcs.R");flocalanc=phase_localanc(localanc,final.flips)
if (verbose) cat("calculating ancestry aware re-phased coancestry curves\n"); acoancs=create_coancs(localanc,dr,"DIP");
#stop("done")
######### GlobeTrotter style curves Mosaic phasing  ###########
#source("noanc.R")
#Mu=a.Mu;rho=a.rho;theta=a.theta;Q=a.Q;alpha=a.alpha;lambda=a.lambda
#o.Mu=a.o.Mu;o.rho=a.o.rho;o.theta=a.o.theta;o.Q=a.o.Q;o.alpha=a.o.alpha;o.lambda=a.o.lambda
#noanc_rephased_localanc=get_ancunaware_localanc() # works off noanc_gfbs
#if (verbose) cat("calculating ancestry unaware re-phased coancestry curves\n"); pcoancs=create_coancs(noanc_rephased_localanc,dr,"DIP")

######## GlobeTrotter style curves original phasing ##########
for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # undo phase flips
source("noanc.R")
source("cleanup.R")
Mu=a.Mu;rho=a.rho;theta=a.theta;Q=a.Q;alpha=a.alpha;lambda=a.lambda
o.Mu=a.o.Mu;o.rho=a.o.rho;o.theta=a.o.theta;o.Q=a.o.Q;o.alpha=a.o.alpha;o.lambda=a.o.lambda
noanc_unphased_localanc=get_ancunaware_localanc() # works off noanc_gfbs
save(file=paste0(resultsdir,"noanc_unphased_localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), 
     noanc_unphased_localanc, flips, g.loc)
if (verbose) cat("calculating ancestry unaware input phasing coancestry curves\n"); coancs=create_coancs(noanc_unphased_localanc,dr,"DIP")

if (verbose) cat("saving final results to file\n")
save(file=paste0(resultsdir,"",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), 
     target, phase.error.locs, o.Mu, o.lambda, o.theta, o.alpha, o.Q, o.rho, Mu, lambda, theta, alpha, Q, rho, L, NUMA, nchrno, chrnos, tol, dr, NL, kLL, RPE, acoancs, coancs)#, pcoancs)

