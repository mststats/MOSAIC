# script that sets default parameters, creates some required objects, and creates functions based on choice of parallelisation strategy (HPC)
require(parallel)
require(Rcpp)
if (HPC==1)
{
  require(ff)
  getdonates<-function(t.donates,NUMI) # note that this is used for donates, donatesl, and donatesr
  {
    donates_chr=list()
    for(ind in 1:NUMI) 
    {
      open(t.donates[[ind]])
      donates_chr[[ind]]=t.donates[[ind]][] # open all donors on this chromosome
      close(t.donates[[ind]])
    }
    donates_chr
  }
}
if (HPC==2)
{
  require(ff)
  getdonates_ind<-function(t.donates) # note that this is used for donates, donatesl, and donatesr
  {
    open(t.donates)
    donates_chr_ind=t.donates[] # open all donors on this chromosome
    close(t.donates)
    donates_chr_ind
  }
}
G<-NULL
if (NUMA==1 & PHASE)
{
  warning("can't do re-phasing on a single haplotype: turning off phasing", immediate.=T)
  PHASE=F
}
gobs<-list()
# some default values
if (!exists("nl")) nl=1000; # maximum number of haps per population 
max.donors<-100;prop.don=0.99; # reasonable defaults. 
total=200;s.total=10; # number of EM iterations at end and w/in reps respectively
PI.total=10; # EM iterations for PI only at start respectively; useful before re-phasing
if (!exists("REPS")) REPS=2*L+1 # maximum number of iterations through thin/phase/EM cycle
# s.M is w/in each iteration of thin / phase / EM and M is after convergence. 
s.M<-0.00;M=0.00 # these are now multiplied by the #gridpoints in each chromosome to determine how many MCMC iterations are run. No longer really needed
eps.lower=log(2) # threshold ratio of new likelihood to old required to phase flip
min.bg=0.1;max.bg=1.0 # default cM length of buffer around phase hunting spikes; hunter will start at min.bg and ramp up to max.bg
resultsdir="RESULTS/"
###############################################################################
PLOT=F # create plots as the code runs? 
initonly=F; # if TRUE just run until ready to find Mu
min.donors=as.integer(10); # take at least this number of haps at each gridpoint
optlevel=3 # used to compile some functions to speed up
tol=1e-8
subNUMA=NUMA # =NUMA=>use all; number of target haps used in no-ancestry initial fit; don't use less than min(2,NUMA)
subNL=100 # #individuals from each panel in no-ancestry initial fit
o.LOG=T;# o.LOG turns on and off reporting of log-like after each thin and each phase (EM always reports as always needed to check convergence)
mcmcprog=F # whether to plot a progress bar for the MCMC phasing; makes for ugly log files!
ffcleanup=T # whether to remove all ff files at the end
if (nchrno==22) samp_chrnos=c(1,3,7,10,15,17) # indices of chromosomes used in no-ancestry initial fit; swap for contiguous 5Mb blocks of all chromosomes?
if (nchrno!=22) samp_chrnos=chrnos[1:5] # just use first 5
if (!exists("GpcM")) GpcM=60 # number of gridpoints per centiMorgan
if (!exists("Ne")) Ne=9e4 # effective population size
if (!exists("S")) S<-rep(NaN,nchrno) # if no limit is set
dr<-1/(GpcM*100) # GpcM is #gridpoints per centiMorgan cM
if (target!="simulated") o.lambda=20 else o.lambda=50 # this is less important now for phasing steps as we get o.lambda from init_Mu 
require(doParallel)
if (is.na(MC)) {
  MC=as.integer(detectCores()/2)
  if (is.na(MC)) {MC=2;warning("using 2 cores as detectCores() has failed",immediate.=T)} # use 2 if can't use detectCores() 
}
if (verbose) cat("using", MC, "cores\n")
registerDoParallel(cores=MC)
g.loc<-list()
maxmiss=maxmatch=0 # these get set in grid.R which is called by read_panels.R
if (!exists("singlePI")) singlePI=F
FLAT=F # set to FALSE to use the recombination rate map. If set to TRUE then map is flattened and one gridpoint per obs is used (this is for debugging purposes). 
source("read_panels.R")
if (!exists("min.donors")) min.donors=10
if (max.donors==NUMP & prop.don<1)
{
  warning("can't use prop.don<1 and all donors: setting prop.don to 1", immediate.=T)
  prop.don=1
}
if (max.donors>NUMP) max.donors<-NUMP # try using less than NUMP
if (min.donors>NUMP) min.donors<-NUMP 
phi.theta<-0.2
if (!exists("absorbrho")) absorbrho=T # should ancestry self-switches be included in rho or diag of PI?
if (!exists("commonrho")) commonrho=T # needs to be false if it includes PI[i,i]
if (!exists("commontheta")) commontheta=T
if (!exists("prethin")) prethin=F
#theta=o.theta<-rep(phi.theta/(phi.theta+NUMP/L), L) # as per Hapmix
theta=o.theta<-rep(phi.theta/(phi.theta+max.donors/L), L) # as per Hapmix
#invsum=1/sum(1/(1:NUMP));o.theta<-rep(0.5*invsum/(NUMP+invsum), L) # Watterson's estimator
rho=o.rho=rep(1-exp(-Ne/(NUMP/L)*dr),L) # similar to HapMix choice but transformed; 1/L as this will include anc self-switches
source("donates.R") # find which haps are useful donors at which gridpoints to which admixed recipients.
source("create_logfile.R")
source("intermediate_calcs.R")
source("mix_hmm.R")
source("phase_funcs.R")
source("initProb.R")
source("EM_updates.R")
if (LL<L) {stop("Can't fit more latent ancs than panels");}
##################################
# now use principal components to set sensible starting values for the copying matrix
alpha=lambda<-list()
for (ind in 1:NUMI)
{
  alpha[[ind]]=rep(1/L,L)
  lambda[[ind]]<-o.lambda
}
rm(Y) # leave in if planning to re-grid
fmutmat<-function(theta, L, maxmiss, maxmatch)
{
  mutmat<-array(NaN, c(L,maxmiss+1,maxmatch+1)) # +1 to include 0
  ltheta<-log(theta)
  lmtheta<-log(1-theta)
  for(l in 1:L) 
    for (i in 1:(maxmiss+1)) # +1 to include 0
      for (j in 1:(maxmatch+1)) # +1 to include 0
	mutmat[l,i,j]=exp(ltheta[l]*(i-1)+lmtheta[l]*(j-1)) # theta^y + (1-theta)^(1-y)
  mutmat[is.nan(mutmat)]=0 # in case of empty ancestries
  mutmat
}
mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
PI<-create_PI(alpha,lambda,L,dr,NUMI)
klabel<-unique(label[KNOWN])
##################################
Mu=matrix(1/kLL,kLL,L)
transitions<-list()
for (ind in 1:NUMI)
  transitions[[ind]]<-s_trans(L,kLL,PI[[ind]],Mu,rho,NL)
if (!EM) total=1
# change next two lines perhaps
o.PI<-PI;o.alpha<-alpha;o.rho<-rho
o.theta<-theta;o.phi.theta<-phi.theta
flips<-list()
for (ind in 1:NUMI) 
{
  flips[[ind]]<-list()
  for (ch in 1:nchrno) 
    flips[[ind]][[ch]]<-rep(F,G[ch]) 
}
phase.error.locs<-list()
