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
if (!exists("GpcM")) GpcM=60 # number of gridpoints per centiMorgan
if (!exists("Ne")) Ne=9e4 # effective population size
if (!exists("S")) S<-rep(NaN,nchrno) # if no limit is set
dr<-1/(GpcM*100) # GpcM is #gridpoints per centiMorgan cM
g.loc<-list()
maxmiss=maxmatch=0 # these get set in grid.R which is called by read_panels.R
#sourceCpp("grid.cpp")
if (!exists("singleQ")) singleQ=F
source("read_panels.R")
if (!exists("prop.don")) prop.don<-1
if (!exists("max.donors")) max.donors<-NUMP # try using less than NUMP
if (!exists("min.donors")) min.donors<-2 
if (max.donors>NUMP) max.donors<-NUMP # try using less than NUMP
if (min.donors>NUMP) min.donors<-NUMP 
phi.theta<-0.2
if (!exists("absorbrho")) absorbrho=T # should ancestry self-switches be included in rho or diag of Q?
if (!exists("commonrho")) commonrho=T # needs to be false if it includes Q[i,i]
if (!exists("commontheta")) commontheta=T
#theta=o.theta<-rep(phi.theta/(phi.theta+NUMP/L), L) # as per Hapmix
theta=o.theta<-rep(phi.theta/(phi.theta+max.donors/L), L) # as per Hapmix
#invsum=1/sum(1/(1:NUMP));o.theta<-rep(0.5*invsum/(NUMP+invsum), L) # Watterson's estimator
rho=o.rho=rep(1-exp(-Ne/(NUMP/L)*dr),L) # similar to HapMix choice but transformed; 1/L as this will include anc self-switches
source("donates.R") # find which haps are useful donors at which gridpoints to which admixed recipients.
source("intermediate_calcs.R")
source("mix_hmm.R")
source("phase_funcs.R")
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
Q<-create_Q(alpha,lambda,L,dr,NUMI)
klabel<-unique(label[KNOWN])
if (USEHAPMIX) source("input_hapmix.R")
##################################
Mu=matrix(1/kLL,kLL,L)
transitions<-list()
for (ind in 1:NUMI)
  transitions[[ind]]<-s_trans(L,kLL,Q[[ind]],Mu,rho,NL)
if (!EM) total=1
# change next two lines perhaps
o.Q<-Q;o.alpha<-alpha;o.rho<-rho
o.theta<-theta;o.phi.theta<-phi.theta
flips<-list()
for (ind in 1:NUMI) 
{
  flips[[ind]]<-list()
  for (ch in 1:nchrno) 
    flips[[ind]][[ch]]<-rep(F,G[ch]) 
}
phase.error.locs<-list()
