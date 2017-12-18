#rm(list=ls());set.seed(1);mean.sim.alpha=c(1,1,8)
require(mosaicpackage)
######################## first set some options ###############################
# important things to set
shargs<-commandArgs(trailingOnly=TRUE) # read in arguments from the command line; 
target=shargs[1]  # e.g. SpainPopn_1
datasource=shargs[2]; # i.e. spanish/ ; note the / at the end is crucial
L=as.integer(shargs[3]) # #ways of admixture
firstind<-as.integer(shargs[4]); # which target individual to start from. If NUMA=2 then only this ind is run
NUMA=as.integer(shargs[5]) # number of target admixed haplotypes
MC=as.integer(shargs[6]) # number of cores to use
chrnos=1:22 # which chromosomes to run on
#tmp=scan("regional.txt",what="character",quiet=T,skip=as.integer(target)-1,nlines=1);target=tmp[2];mask=tmp[-(1:2)] # read regional.txt for inputs; target is which line to run
#chrnos=1:22;firstind=1;NUMA=440;nl=200;L=2;datasource="HGDP/";target="NorthAfrican";ANC=T # use this line for testing
#chrnos=10:15;prop.missing=0.0;firstind=1;target="simulated";NUMA=8;L=3;datasource="HGDP/";RPE=2.0;ANC=T;#ANC=c("Ireland","Yoruba"); # use this line for testing
chrnos=10:10;prop.missing=0.0;firstind=1;target="simulated";NUMA=2;L=2;datasource="HGDP/";RPE=0.0;ANC=T;
#chrnos=10:22;firstind=1;NUMA=16;nl=200;L=5;datasource="HGDP/";target="SanKhomani";ANC=NULL # use this line for testing
#chrnos=10:15;firstind=1;NUMA=8;nl=200;L=2;datasource="HGDP/";target="Hazara";ANC=T # use this line for testing
#chrnos=21:22;firstind=1;NUMA=4;L=2;datasource="spanish/";target="SpainPopn_2";mask=c("Spain","Portugal") # use this line for testing
nchrno=length(chrnos) # number of chromosomes for these target haplotypes
HPC=2 # whether to use ff() chromosome-by-chromosome (HPC=1) or chromosomeXind-by-chromsomeXind(HPC=2) or not at all (HPC=F);
ffpath="/dev/shm/" # location of fast-files
# somewhat less important things to set
set.seed(1)
if (!exists("PHASE")) PHASE=T
nl<-1000; # maximum number of haps per population 
max.donors<-100;prop.don=0.99; # reasonable defaults. 
total=200;s.total=10; # number of EM iterations at end and w/in reps respectively
Q.total=10; # EM iterations for Q only at start respectively; useful before re-phasing
REPS=2*L+1 # maximum number of iterations through thin/phase/EM cycle
# s.M is w/in each iteration of thin / phase / EM and M is after convergence. 
s.M<-0.00;M=0.00 # these are now multiplied by the #gridpoints in each chromosome to determine how many MCMC iterations are run. No longer really needed
eps.lower=log(2) # threshold ratio of new likelihood to old required to phase flip
min.bg=0.1;max.bg=1.0 # default cM length of buffer around phase hunting spikes; hunter will start at min.bg and ramp up to max.bg
resultsdir="RESULTS/"
###############################################################################

# the rest are mostly used in debugging, comparing answers w/ Hapmix, etc
if (!exists("ANC")) ANC=NULL; # no a-priori knowledge of which panels to use for which ancestries
min.donors=as.integer(10); # take at least this number of haps at each gridpoint
if (!exists("RPE")) RPE=0; # Rate of Phase Errors: used to purposefully introduce phasing errors that change ancestry in simulated example (when RPE>0)
min.Mu.ratio<-0.4 # drop groups of haps that are each used less than min.Mu.ratio/NUMP at max across ancs; i.e. min.Mu.ratio times the prior
FLAT=F # set to FALSE to use the recombination rate map. If set to TRUE then map is flattened and one gridpoint per obs is used (this is for debugging purposes). 
optlevel=3 # used to compile some functions to speed up
USEHAPMIX=F # don't use Hapmix EM values to initialise
tol=1e-8
if (!exists("OUTHAPMIX")) OUTHAPMIX=F # output parameters to be used by Hapmix?
verbose=T # print certain statements of progress as algorithm runs?
PLOT=F # create plots as the code runs? 
EM=T # run EM algorithm?
doMu=T # update copying matrix parameters?
doQ=T # update ancestry switching parameters parameters?
dorho=T # update recombination w/in same ancestry parameters? 
dotheta=T # update error / mutation parameters?
initonly=F; # if TRUE just run until ready to find Mu
earlydrop=T # drop some groups immediately after no-ancestry model fit?
if (nchrno==22) samp_chrnos=c(1,3,7,10,15,17) # indices of chromosomes used in no-ancestry initial fit; swap for contiguous 5Mb blocks of all chromosomes?
if (nchrno!=22) samp_chrnos=chrnos[1:5] # just use first 5
subNUMA=NUMA # =NUMA=>use all; number of target haps used in no-ancestry initial fit; don't use less than min(2,NUMA)
subNL=100 # #individuals from each panel in no-ancestry initial fit
dropfast=T # drop multiple panels at once due to unclear signal? 
drop.threshold=-0.005 # drop all panels below this expected change in r^2; use NaN to turn this off
drop.threshold=NaN;min.Mu.ratio=0; # no group dropping at all
maxdrops=5 # maximum number of panels to drop at any one time
o.LOG=T;# o.LOG turns on and off reporting of log-like after each thin and each phase (EM always reports as always needed to check convergence)
mcmcprog=F # whether to plot a progress bar for the MCMC phasing; makes for ugly log files!
ffcleanup=T # whether to remove all ff files at the end
if (target!="simulated") o.lambda=20 else o.lambda=50 # this is less important now for phasing steps as we get o.lambda from init_Mu 
require(doParallel)
if (is.na(MC)) {
  MC=as.integer(detectCores()/2)
  if (is.na(MC)) {MC=2;warning("using 2 cores as detectCores() has failed",immediate.=T)} # use 2 if can't use detectCores() 
}
if (verbose) cat("using", MC, "cores\n")
registerDoParallel(cores=MC)

if (!exists("prethin")) prethin=F # don't use this by default
source("setup.R") # includes introducing phase errors if RPE>0
if (max.donors==NUMP & prop.don<1)
{
  warning("can't use prop.don<1 and all donors: setting prop.don to 1", immediate.=T)
  prop.don=1
}
old.runtime<-as.numeric(Sys.time())
o.total=total
writelog<-function(alg) # single consistent function to write to EMlogfile
  write(file=EMlogfile,c(alg,signif(diff.time,4),signif(t(Mu),4),signif(rho,4),c(sapply(Q, function(x) signif(t(x),4))),
			 sapply(alpha, function(x) signif(x,4)),sapply(lambda,function(x) round(x,4)),signif(theta,4),round(cloglike,4)),ncol=length(lognames),append=T)

o.drop.threshold=drop.threshold;drop.threshold=NaN;
eps=log(1.01) # i.e. a 1% increase in relative likelihood
get_switches=F
if (prethin)
  source("pre_all_donates.R")
Mu<-matrix(rep(1/kLL,L*kLL),kLL);for (ind in 1:NUMI) alpha[[ind]]=rep(1/L,L) # flatten out w.r.t. ancestry
source("noanc.R") # always need to run noanc.R b/c need good paras for init_Mu
source("initProb.R")
runtime<-as.numeric(Sys.time())
if (kLL>L) # otherwise can't cluster kLL things into L clusters
{
  source("init_Mu.R")
  get_switches=T;LOG=F;source("all_donates.R") # use this to get #switches in noanc model w/o writelog
  windowed_copying<-window_chunks(nswitches=noanc_gswitches,ww=0.5,verbose=verbose) # similar in the same windows
  if (initonly)
  {
    save(file=paste0(resultsdir,"init_",target,"_",L, "way_", firstind, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,".RData"), 
     target, panels, Mu, rho, theta, alpha, lambda, Q, windowed_copying, L, NUMA, nchrno, chrnos, g.loc, tol, dr, NL, kLL, RPE)
    source("cleanup.R")
    stop("saving initialisation and quitting")
  }
  tmp<-cluster_windows(windowed_copying,PLOT=F,t.L=L,verbose=verbose)
  Mu<-tmp$Mu
  alpha<-tmp$alpha
  Q<-tmp$Q
  all.o.lambda=tmp$lambda
  #lambda=as.list(sapply(tmp$lambda,mean,na.rm=T)) # doesn't work well!
  lambda=o.lambda # re-use o.lambda from above
  rm(windowed_copying,tmp)
  o.Mu<-Mu;o.alpha<-alpha;o.lambda=lambda;o.Q=Q # these are the official starting ancestry related parameters now
}
get_switches=F
mutmat<-fmutmat(theta, L, maxmiss, maxmatch); for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,Q[[ind]],Mu,rho,NL)
rm(noanc_gswitches)
if (OUTHAPMIX)
{
  for (ch in 1:nchrno)
    output_hapmix_initial(ch)
}
o.M<-M;M<-s.M
if (EM) source("create_logfile.R")
LOG=F
if (earlydrop & min.Mu.ratio>0) 
  source("rem_groups.R") # this will do all_donates.R
LOG=o.LOG;cloglike=NaN 
if (!(earlydrop & min.Mu.ratio>0))
  source("all_donates.R") # decide on donor set using initial parameters

############ a few Q only updates first; very useful to do before first re-phasing
if (Q.total>0)
{
  o.doMu=doMu;o.dotheta=dotheta;o.dorho=dorho;o.doQ=doQ;doQ=T;dorho=dotheta=doMu=F;
  if (verbose) cat("Inferring ancestry switching rates holding other parameters fixed\n");
  total=Q.total;source("mosaic.R")
  if (!absorbrho | !commonrho | !commontheta) source("all_donates.R") 
  doMu=o.doMu;dorho=o.dorho;dotheta=o.dotheta;doQ=o.doQ
}

total<-s.total # number of EM steps with repeated loop
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

##### FLAG ##### 
#if (verbose) cat("saving final results to file\n")
#save(file=paste0(resultsdir,"",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), 
#     target, phase.error.locs, o.Mu, o.lambda, o.theta, o.alpha, o.Q, o.rho, Mu, lambda, theta, alpha, Q, rho, L, NUMA, nchrno, chrnos, tol, dr, NL, kLL, RPE)#, acoancs, coancs)#, pcoancs)
################ 

#source("plot_funcs.R");flocalanc=phase_localanc(localanc,final.flips)
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

