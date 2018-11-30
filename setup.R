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
# function that sets default parameters, creates some required objects, and creates functions based on choice of parallelisation strategy (HPC)
setup_data_etc=function(t.NUMA,t.target,t.chrnos,t.pops,L,datasource,EM,gens,ratios,MC,
  # some default values
  verbose=T,
  HPC=2, # whether to use ff() chromosome-by-chromosome (HPC=1) or chromosomeXind-by-chromsomeXind(HPC=2) or not at all (HPC=F);
  PHASE=TRUE,
  nl=1000, # maximum number of haps per population 
  max.donors=100, # maximum number of donors to consider per t.target at each gridpoint
  min.donors=10L, # take at least this number of haps at each gridpoint
  prop.don=0.99, # stop taking donors once this degree of probability is reached
  s.total=10, # maximum number of EM iterations in each round
  total=200, # maximum number of EM iterations in final round
  PI.total=0, # EM iterations for PI only at start; useful before re-phasing
  REPS=0, # maximum number of iterations through thin/phase/EM cycle
  s.M=0.00, # w/in each iteration of thin / phase / EM 
  M=0.00, # MCMC phasing after convergence  
  # above two are multiplied by the #gridpoints in each chromosome to determine how many MCMC iterations are run. No longer really needed
  eps.lower=log(2), # threshold ratio of new likelihood to old required to phase flip
  min.bg=0.1, 
  max.bg=1.0, # default cM length of buffer around phase hunting spikes; hunter will start at min.bg and ramp up to max.bg
  LOG=TRUE, # LOG turns on and off reporting of log-like after each thin and each phase (EM always reports as always needed to check convergence)
  mcmcprog=F, # whether to plot a progress bar for the MCMC phasing; makes for ugly log files!
  GpcM, # number of gridpoints per centiMorgan
  Ne=9e4, # effective population size
  mask=NULL,
  singlePI=FALSE,
  absorbrho=TRUE, # should ancestry self-switches be included in rho or diag of PI?
  commonrho=TRUE, # needs to be false if it includes PI[i,i]
  commontheta=TRUE,
  prethin=FALSE,
  resultsdir="MOSAIC_RESULTS/") # where to store results files
{
  if (!file.exists(resultsdir))
    dir.create(file.path(getwd(), resultsdir))
  if (REPS==0) REPS=2*L+1 # maximum number of iterations through thin/phase/EM cycle
  REPS=ifelse(EM, REPS, 1) # no need for more than 1 if no EM parameter changes
  t.nchrno=length(t.chrnos)
  ans=list() # build a list to store resulting data, parameters, etc
  S<-rep(NaN,t.nchrno) # if no limit is set
  if (t.NUMA==1 & PHASE)
  {
    warning("can't do re-phasing on a single haplotype: turning off phasing", immediate.=T)
    PHASE=FALSE
  }
  # some defaults will get returned also
  ans$resultsdir=resultsdir
  ans$PHASE=PHASE
  ans$HPC=HPC
  ans$GpcM=GpcM
  ans$LOG=LOG
  ans$mcmcprog=mcmcprog
  ans$absorbrho=absorbrho
  ans$commonrho=commonrho
  ans$commontheta=commontheta
  ans$prethin=prethin
  ans$s.M=s.M
  ans$M=M
  ans$PI.total=PI.total
  ans$s.total=s.total
  ans$REPS=REPS
  ans$eps.lower=eps.lower
  ans$min.bg=min.bg
  ans$max.bg=max.bg
  if (t.nchrno==22) ans$samp_chrnos=c(1,3,7,10,15,17) # indices of chromosomes used in no-ancestry initial fit; swap for contiguous 5Mb blocks of all chromosomes?
  if (t.nchrno!=22) ans$samp_chrnos=t.chrnos[1:5] # just use first 5
  if (length(ans$samp_chrnos)>t.nchrno) ans$samp_chrnos=t.chrnos # use all if try to use too many
  ans$dr<-1/(ans$GpcM*100) # GpcM is #gridpoints per centiMorgan cM
  if (gens==0) { # if not provided
    if (t.target!="simulated") gens=10 else gens=50 # this is less important now for phasing steps as we get lambda from init_Mu 
  } # note that gens is a single date (not necessarily true for multiway L>2 admixture; EM will correct this quickly if on)
  if (is.null(ratios)){  # if not provided
    ratios=rep(1/L,L)
  } else {
    if (length(ratios)==L) {ratios=ratios/sum(ratios)} else {
      cat(ratios, " supplied as vector of mixing group ratios\n")
      ratios=rep(1/L,L)
      if (!EM)
        warning("length of ancestry ratios must be equal to number of mixing groups: using 1/", L, " for each group", immediate.=T)
    }
  }
  if (!is.null(t.pops)) {
    if (t.target=="simulated" & length(t.pops)>L & length(t.pops)<(2*L)) {
      warning("Provide at least (2 X #ancestries) named groups; first #ancestries to simulate from and rest to use as donor groups.
  Therefore using all other available groups as donors\n", immediate.=T)
      t.pops=t.pops[1:L]
    }
    if (t.target!="simulated" & length(t.pops)<L) {
      warning("Need at least " , L, " named groups to use as donors; using all available donor groups\n", immediate.=T)
      t.pops=NULL
    }
  }
  if (!EM & L>2)
    warning("Turning off EM and specifying a single mixing date is not advised", immediate.=T)
  if (MC==0) {
    MC=as.integer(detectCores()/2)
    if (is.na(MC)) {MC=2;warning("using 2 cores as detectCores() has failed",immediate.=T)} # use 2 if can't use detectCores() 
  }
  if (verbose) cat("using", MC, "cores\n")
  registerDoParallel(cores=MC)
  ans$FLAT=FALSE # FALSE to use the recombination rate map. If set to TRUE then map is flattened and one gridpoint per obs is used (this is for debugging purposes). 
  tmp=read_panels(datasource, t.target, t.chrnos, t.NUMA, L, t.pops, nl, ans$FLAT, ans$dr, gens, resultsdir, mask=mask, ratios=ratios) 
  if (verbose) cat("\nFitting model to ", tmp$NUMI, " ", t.target, " ", L, "-way admixed target individuals using ", tmp$kLL, " panels\n", sep="")
  if (verbose) cat("EM inference is ", ifelse(EM, "on", "off"), " and re-phasing is ", ifelse(PHASE, "on", "off"), "\n")
  ans$maxmatch=tmp$maxmatch;ans$maxmiss=tmp$maxmiss;ans$umatch=tmp$umatch;ans$d.w=tmp$d.w;ans$t.w=tmp$t.w;ans$g.loc=tmp$g.loc;ans$gobs=tmp$gobs
  ans$NUMP=tmp$NUMP;LL=tmp$LL;ans$NUMA=tmp$NUMA;ans$NUMI=tmp$NUMI;ans$label=tmp$label;ans$KNOWN=tmp$KNOWN;ans$kLL=tmp$kLL;ans$NL=tmp$NL;ans$G=tmp$G;ans$NN=tmp$NN;
  ans$maxmatchsize=tmp$maxmatchsize;ans$panels=tmp$panels
  if (t.target=="simulated")
    ans$g.true_anc=tmp$g.true_anc
  if (max.donors==ans$NUMP & prop.don<1)
  {
    warning("can't use prop.don<1 and all donors: setting prop.don to 1", immediate.=T)
    prop.don=1
  }
  if (max.donors>ans$NUMP) max.donors<-ans$NUMP # try using less than NUMP
  if (min.donors>ans$NUMP) min.donors<-ans$NUMP 
  ans$min.donors=min.donors
  phi.theta<-0.2
  ans$theta=o.theta<-rep(phi.theta/(phi.theta+max.donors/L), L) # as per Hapmix
  #invsum=1/sum(1/(1:ans$NUMP));o.theta<-rep(0.5*invsum/(ans$NUMP+invsum), L) # Watterson's estimator
  ans$rho=o.rho=rep(1-exp(-Ne/(ans$NUMP/L)*ans$dr),L) # similar to HapMix choice but transformed; 1/L as this will include anc self-switches
  if (LL<L) {stop("Can't fit more latent ancs than panels");}
  ##################################
  # now use principal components to set sensible starting values for the copying matrix
  ans$alpha=ans$lambda<-list()
  for (ind in 1:ans$NUMI)
  {
    ans$alpha[[ind]]=ratios
    ans$lambda[[ind]]=gens
  }
  mutmat<-fmutmat(ans$theta, L, ans$maxmiss, ans$maxmatch)
  ans$PI<-create_PI(ans$alpha,ans$lambda,L,ans$dr,ans$NUMI)
  ##################################
  ans$Mu=matrix(1/ans$kLL,ans$kLL,L)
  ans$transitions<-list()
  for (ind in 1:ans$NUMI)
    ans$transitions[[ind]]<-s_trans(L,ans$kLL,ans$PI[[ind]],ans$Mu,ans$rho,ans$NL)
  if (!EM) total=1
  ans$total=total
  # change next two lines perhaps
  ans$o.PI<-ans$PI;ans$o.alpha<-ans$alpha;ans$o.rho<-ans$rho;ans$o.lambda=ans$lambda
  ans$o.theta<-ans$theta;ans$o.phi.theta<-phi.theta
  ans$flips<-list()
  for (ind in 1:ans$NUMI) 
  {
    ans$flips[[ind]]<-list()
    for (ch in 1:t.nchrno) 
      ans$flips[[ind]][[ch]]<-rep(F,ans$G[ch]) 
  }
  ans$prop.don=prop.don
  ans$max.donors=max.donors
  return(ans)
}
