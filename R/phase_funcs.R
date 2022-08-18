# functions used in re-phasing target haplotypes based on current MOSAIC fit
flip.ll<-function(g, t.ch, t.A, ndonorsg, t.NUMP, NNL2, t.G, t.scalefactor, t.scalefactorb, t.fors, t.backs)
{
  g=g-1
  tmpvec=1:(ndonorsg*t.A)
  ans=0
  for (h in 1:2)
  {
    ans=ans-sum(log(t.scalefactor[[h]][1:g]))
    ans=ans-sum(log(t.scalefactorb[[h]][g:t.G[t.ch]]))
    tmp=sum(t.fors[[h]][(g-1)*NNL2+tmpvec]*t.backs[[3-h]][(g-1)*NNL2+tmpvec])
    ans=ans+log(tmp)-log(t.NUMP)-log(t.A)
  }
  ans
}

all.flip.ll<-function(g, t.ch, t.A, ndonorsg, t.NUMP, NNL2, t.G, cumsum.log.fors, rev.cumsum.log.backs, t.fors, t.backs)
{
  tmpvec=1:(ndonorsg*t.A)
  ans=0
  for (h in 1:2)
  {
    ans=ans+cumsum.log.fors[[h]][g-1]
    ans=ans+rev.cumsum.log.backs[[h]][t.G[t.ch]-g+2]
    tmp=sum(t.fors[[h]][(g-2)*NNL2+tmpvec]*t.backs[[3-h]][(g-2)*NNL2+tmpvec])
    ans=ans+log(tmp)-log(t.NUMP)-log(t.A)
  }
  ans
}


r.create.proposal<-function(t.ch,t.G,t.A,t.NUMP,t.max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll) 
{
  NNL2=t.max.donors*t.A
  cumsum.log.fors=list(-cumsum(log(t.scalefactor[[1]])),-cumsum(log(t.scalefactor[[2]])))
  rev.cumsum.log.backs=list(-cumsum(log(rev(t.scalefactorb[[1]]))),-cumsum(log(rev(t.scalefactorb[[2]]))))
  g.flip.ll<-function(g) all.flip.ll(g, t.ch, t.A, t.ndonors[g], t.NUMP, NNL2, t.G, cumsum.log.fors, rev.cumsum.log.backs, t.fors, t.backs)
  ll=c(ind.c.ll,vapply(2:t.G[t.ch], g.flip.ll,0))
  ll=ll-ind.c.ll
  ll[is.nan(ll)|is.infinite(ll)]=-1 # can happen if no path through
  ll
}
create.proposal<-cmpfun(r.create.proposal,list(optimize=3)) # 
#create.proposal<-r.create.proposal
# H1 gets: fors[1] to g, backs[2] from g and new fors from g, new backs to g
# H2 gets: fors[2] to g, backs[1] from g and new fors from g, new backs to g
rephaser<-function(t.ch, t.G, t.A, hap, g, t.NUMP, t.NUMA, t.kLL, t.transitions, t.umatch, t.maxmatchsize, t.dw, t.tw, t.gobs, t.fors, t.sumfors, t.backs, t.scalefactor, t.scalefactorb, t.flips, t.initProb, 
		   t.label, t.ndonors, t.donates, t.donatesl, t.donatesr) # returns log-likelihood of data when g is flipped for jth individual
{
  c.fors<-c.sumfors<-c.backs<-c.scalefactor<-c.scalefactorb<-list()
  for (h in 1:2) 
  {
    c.fors[[h]]<-clone(t.fors[[h]])
    c.sumfors[[h]]<-clone(t.sumfors[[h]])
    c.backs[[h]]<-clone(t.backs[[h]])
    c.scalefactor[[h]]<-clone(t.scalefactor[[h]])
    c.scalefactorb[[h]]<-clone(t.scalefactorb[[h]])
  }
  tmp.backs<-c.backs
  tmp.scalefactorb<-c.scalefactorb
  # swap the backward probabilities and scaling factors
  c.backs[[1]]<-tmp.backs[[2]];c.backs[[2]]<-tmp.backs[[1]]
  c.scalefactorb[[1]]<-tmp.scalefactorb[[2]];c.scalefactorb[[2]]<-tmp.scalefactorb[[1]]
  for (h in 1:2)
  {
    cppforward(hap[h],t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,g-1,t.G[t.ch],t.G[t.ch],t.transitions, t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.initProb[hap[h],],t.label,t.ndonors,t.donates,t.donatesl,t.flips,c.fors[[h]],c.sumfors[[h]],c.scalefactor[[h]]) # g-1
    cppbackward(hap[h],t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,g-1,t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.label,t.ndonors,t.donates,t.donatesr,t.flips,c.backs[[h]],c.scalefactorb[[h]]) # g-1
  }
  return(list(fors=c.fors, sumfors=c.sumfors, backs=c.backs, scalefactor=c.scalefactor, scalefactorb=c.scalefactorb))
}

r_phase_hunt<-function(t.eps.lower, t.ch, t.G, t.A, t.GpcM, t.ind, t.flips, verbose, t.ndonors, t.NUMP, t.NUMA, t.kLL, t.max.donors, t.donates, t.donatesl, t.donatesr, 
		       t.transitions, t.umatch, t.maxmatchsize, t.dw, t.tw, t.gobs, t.mutmat, t.maxmiss, t.initProb, t.label, lim=10, minbg=0.1, maxbg=1, mult=1.5)
{
  hap<-c(t.ind*2-1,t.ind*2) # t.ind indexes over genotypes, hap the two haplotypes
  # fb calcs moved to here to avoid storing all fors, backs, etc
  t.fors<-list();t.sumfors<-list();t.scalefactor<-list()
  t.backs<-list();t.scalefactorb<-list()
  THIN=ifelse(t.max.donors==t.NUMP,FALSE,TRUE)
  for (h in 1:2)
  {
    k=hap[h]
    t.fors[[h]]<-rep(0,t.G[t.ch]*t.max.donors*t.A);t.sumfors[[h]]<-matrix(0,t.G[t.ch],t.A);t.scalefactor[[h]]<-rep(0,t.G[t.ch])
    cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.initProb[k,],t.label,
	       t.ndonors,t.donates,t.donatesl,t.flips,t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
    t.backs[[h]]<-rep(0,t.G[t.ch]*t.max.donors*t.A);t.scalefactorb[[h]]<-rep(0,t.G[t.ch])
    cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.label,
		t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
  }
  ind.orig.ll<-ind.max.ll<-ind.c.ll<- -sum(log(t.scalefactor[[1]]))-sum(log(t.scalefactor[[2]])) # LL for a single individual
  ind.c.v<-create.proposal(t.ch,t.G,t.A,t.NUMP,t.max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll)
  c.probs<-1/(1+exp(-ind.c.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
  site.c.probs<-c.probs/sum(c.probs)
  nflips<-1
  t.nflips<-0
  iters=1
  end.eps.lower=t.eps.lower
  if (max(ind.c.v)>t.eps.lower)
    while (nflips>0 & iters<lim) {
      #t.eps.lower=end.eps.lower+0.5*log(lim/iters) # fewer flips at the start in phase hunter, down to 2 as lower threshold for a phase flip
      bg=((2*plogis(mult*iters-mult))-1)*(maxbg-minbg)+minbg # smaller band at the start, increasing quickly
      BG=as.integer(bg*t.GpcM) # bg is number of centiMorgans to each side of a proposed phase flip around which don't try more flips
      if (max(ind.c.v)<t.eps.lower) 
	break;
      old.ll<-ind.c.ll
      cand.l<-(ind.c.v>t.eps.lower) # TRUE for all those that would have log-like change greater than t.eps.lower if flipped
      cand.u<-which(cand.l) # which are candidates
      cand<-NULL # take the max spike, block off around it, repeat until none left
      while (max(ind.c.v[cand.u])>t.eps.lower & sum(cand.l)>0) { # while some are still worth flipping
	tmp<-which.max(ind.c.v[cand.u]) # identify maximum
	cand<-c(cand,cand.u[tmp]) # add to the vector of candidates
	cand.l[max(1,cand.u[tmp]-BG):min(cand.u[tmp]+BG,t.G[t.ch])]<-FALSE # remove block of BG gridpoints to each side
	cand.u<-which(cand.l)
	if (length(cand.u)==0)
	  break
      }
      cand<-sort(cand,method="quick") # required?
      nflips<-length(cand)
      t.nflips<-t.nflips+nflips
      if (nflips==0)
      {
	break;
      }
      for (g in cand) 
	t.flips[g:t.G[t.ch]]<-!t.flips[g:t.G[t.ch]]
      # full refit required here; potentially lots of re-phasing
      for (h in 1:2) {
	k<-hap[h]
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors,t.donates,t.donatesl,t.flips,
		   t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.label,t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
      }
      ind.c.ll<- -sum(sum(log(t.scalefactor[[1]])))-sum(sum(log(t.scalefactor[[2]])))
      if (ind.c.ll==old.ll) # just break
	break;
      if (ind.c.ll<old.ll) { # only flip max one
	t.nflips<-t.nflips-nflips+1
	nflips=1
	for (g in cand) 
	  t.flips[g:t.G[t.ch]]<-!t.flips[g:t.G[t.ch]] # flip all candidates back
	g=which.max(ind.c.v)
	t.flips[g:t.G[t.ch]]<-!t.flips[g:t.G[t.ch]] # flip max candidate only
	# full refit required again here; can't simply re-use half of forwards and backwards as these have changed in attempting multiple flips
	for (h in 1:2) {
	  k<-hap[h]
	  t.ind<-as.integer((k+1)/2)
	  cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors,t.donates,t.donatesl,t.flips,
		     t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
	  cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.label,t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
	}
	ind.c.ll<- -sum(sum(log(t.scalefactor[[1]])))-sum(sum(log(t.scalefactor[[2]])))
	if (is.nan(ind.c.ll)|is.infinite(ind.c.ll)) ind.c.ll=old.ll+ind.c.v[g];
	#break; # actually we often find more again after a single flip
      } 
      # if pass above checks then loglikelihood has improved and we try again based on a new proposal
      ind.c.v<-create.proposal(t.ch,t.G,t.A,t.NUMP,t.max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll)
      if (max(ind.c.v)<=0)
	break;
      if (verbose) cat("iter:", iters, "BG:", BG, "flips:", nflips, ind.c.ll, "\n")
      iters=iters+1
    }
  ind.max.ll<-ind.c.ll
  c.probs<-1/(1+exp(-ind.c.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
  site.c.probs<-c.probs/sum(c.probs) 
  if (verbose) cat(t.nflips, "phase flips made: log-likelihood ", ind.orig.ll, "->", ind.c.ll, " on Chr ", chrnos[t.ch], " for ind ", t.ind, 
		   "after ", iters, "hunting passes", "\n")
  # best so far so save to a file
  return(list(ind.max.ll=ind.max.ll,ind.orig.ll=ind.orig.ll,ind.max.flips=t.flips,nflips=t.nflips,niters=iters))
}

#phase_hunt<-cmpfun(r_phase_hunt,list(optimize=3)) # 
phase_hunt<-r_phase_hunt
r_phase_mcmc<-function(t.ch, t.G, t.A, t.ind, M, t.NUMP, t.NUMA, t.kLL, t.max.donors, t.initProb, t.label, t.flips, verbose, t.ndonors, t.donates, t.donatesl, t.donatesr, 
		       t.transitions, t.umatch, t.maxmatchsize, t.dw, t.tw, t.gobs, t.mutmat, t.maxmiss, mcmcprog, mcmchill=TRUE) 
{
  M=as.integer(M*t.G[t.ch])
  NNL2=t.max.donors*t.A
  if (mcmcprog) pb<-txtProgressBar(min=1,max=M,style=3)
  # fb calcs moved to here to avoid storing all fors, backs, etc
  t.fors<-list();t.sumfors<-list();t.scalefactor<-list()
  t.backs<-list();t.scalefactorb<-list()
  THIN=ifelse(t.max.donors==t.NUMP,FALSE,TRUE)
  for (h in 1:2)
  {
    k=hap[h]
    t.fors[[h]]<-rep(0,t.G[t.ch]*t.max.donors*t.A);t.sumfors[[h]]<-matrix(0,t.G[t.ch],t.A);t.scalefactor[[h]]<-rep(0,t.G[t.ch])
    cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.initProb[k,],t.label,
	       t.ndonors,t.donates,t.donatesl,t.flips,t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
    t.backs[[h]]<-rep(0,t.G[t.ch]*t.max.donors*t.A);t.scalefactorb[[h]]<-rep(0,t.G[t.ch])
    cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[t.ch],t.G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.mutmat,t.maxmiss,t.label,
		t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
  }
  ind.orig.ll<-ind.max.ll<-ind.c.ll<- -sum(log(t.scalefactor[[1]]))-sum(log(t.scalefactor[[2]])) # LL for a single individual
  ind.max.flips<-t.flips
  ind.mcmc.ll<-ind.c.ll # running record of log-likelihood for this individual
  hap<-c(t.ind*2-1,t.ind*2) # t.ind indexes over genotypes, hap the two haplotypes
  ind.mcmc.acc<-0
  ind.c.v<-create.proposal(t.ch,t.G,t.A,t.NUMP,t.max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll)
  c.probs<-1/(1+exp(-ind.c.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
  if (all(c.probs==0)) {cat("none left to propose at this hill climbing rate\n");return();}
  site.c.probs<-c.probs/sum(c.probs)
  if (verbose) cat("running MCMC chain...\n")
  p.flips<-t.flips
  for (m in 1:M) {
    #if (any(t.flips)) cat(m,1,"\n")
    ind.mcmc.acc<-c(ind.mcmc.acc,0) # default to zero, changes below if accepted
    #if (max(c.probs)==0) cat(range(ind.c.v), range(p.v), "\n") 
    g<-sample(2:t.G[t.ch],1,prob=c.probs[-1]) 
    # given j, choose to propose to flip or not flip with probability c.probs[j]
    p.flip<-rbinom(1,prob=c.probs[g],size=1) # this will be low. flip==1 => flip it
    #cat(ind.c.ll, mean(c.probs), c.probs[g], p.flip, "\n")
    if (p.flip==1) {
      # now update proposal distribution
      ##########################################################################################
      p.flips<-t.flips
      p.flips[g:t.G[t.ch]]<-!t.flips[g:t.G[t.ch]] # could have been flipped before in which case will get flipped back
      rephased<-rephaser(t.ch,t.G,t.A,hap,g,t.NUMP,t.NUMA,t.kLL,t.transitions,t.umatch,t.dw,t.tw,t.gobs,(t.fors),(t.sumfors),(t.backs),(t.scalefactor),(t.scalefactorb),
			 (p.flips),t.initProb,t.ndonors,t.donates,t.donatesl,t.donatesr)
      #gc()
      n.ll<-flip.ll(g, t.ch, t.A, t.ndonors[g], t.NUMP, NNL2, t.G, t.scalefactor, t.scalefactorb, t.fors, t.backs) # LL for a single individual
      p.v<-create.proposal(t.ch,t.G,t.A,t.NUMP,t.max.donors,rephased$fors,rephased$sumfors,rephased$backs,rephased$scalefactor,rephased$scalefactorb,t.ndonors,n.ll)
      # below is super fast but incorrect; leads to bad proposals in the long run
      #p.v<-ind.c.v;p.v[g]=flip.ll(g, t.ch, t.A, t.ndonors[g], t.NUMP, NNL2, t.G, rephased$scalefactor, rephased$scalefactorb, rephased$fors, rephased$backs)-n.ll
      p.probs<-1/(1+exp(-p.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
      if (all(p.probs==0)) {cat("none left to propose at this hill climbing rate\n");break();}
      site.p.probs<-p.probs/sum(p.probs) 
      ##########################################################################################
      acc<-exp(min(0,n.ll-ind.c.ll+log(p.probs[g])-log(c.probs[g])+log(site.p.probs[g])-log(site.c.probs[g]))) 
      tmp.ll<- -sum(log(rephased$scalefactor[[1]]))-sum(log(rephased$scalefactor[[2]])) # LL for a single individual
      rp<-runif(1,0,1) # random uniform prob unless hill climbing (then only take improvements)
      if (mcmchill & n.ll<ind.c.ll)
      {
	acc=0 # set acceptance prob to zero for lower log-likelihood
	c.probs[g]=0 # don't propose this flip again
      }
      if (acc>=rp) {
	ind.mcmc.acc[m+1]<-1
	t.flips[g:t.G[t.ch]]<-!t.flips[g:t.G[t.ch]] 
	for (h in 1:2) {
	  t.fors[[h]]<-clone(rephased$fors[[h]])
	  t.sumfors[[h]]<-clone(rephased$sumfors[[h]])
	  t.backs[[h]]<-clone(rephased$backs[[h]])
	  t.scalefactor[[h]]<-clone(rephased$scalefactor[[h]])
	  t.scalefactorb[[h]]<-clone(rephased$scalefactorb[[h]])
	}
	ind.c.ll<-n.ll # set next loglikelihood to new loglikelihood
	ind.c.v<-p.v
	c.probs<-p.probs
	site.c.probs<-site.p.probs
	if (ind.c.ll>ind.max.ll) { # new maximum likelihood value
	  ind.max.ll<-ind.c.ll
	  ind.max.flips<-t.flips
	}
      }
    }
    ind.mcmc.ll<-c(ind.mcmc.ll,ind.c.ll) # running record of log-likelihood for this individual
    if (mcmcprog) setTxtProgressBar(pb, m)
  }
  if (mcmcprog)
    close(pb)
  if (verbose) 
    cat(100*mean(ind.mcmc.acc), "% accepted: log-likelihood ", ind.orig.ll, " -> ", ind.max.ll, " on Chr ", chrnos[t.ch], " for ind ", t.ind, "\n", sep="")
  return(list(ind.max.flips=ind.max.flips,ind.mcmc.ll=ind.mcmc.ll,ind.max.ll=ind.max.ll,ind.mcmc.acc=ind.mcmc.acc))
}
#phase_mcmc<-cmpfun(r_phase_mcmc,list(optimize=3)) # 
phase_mcmc<-r_phase_mcmc


# function to run the fast phase-hunting algorithm to re-phase target genomes using MOSAIC fit
phase_hunt_all=function(t.donates, t.donatesl, t.donatesr, t.ndonors, t.NUMP, t.G, t.A, t.GpcM, t.max.donors, t.nchrno, t.NUMA, t.NUMI, t.kLL, t.flips, t.eps.lower, 
			t.transitions, t.umatch, t.maxmatchsize, t.d.w, t.t.w, t.gobs, t.mutmat, t.maxmiss, 
			t.initProb, t.label, t.min.bg, t.max.bg,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,
			t.old.runtime, t.logfile, t.HPC=2, verbose=TRUE, t.LOG=TRUE) {
  orig.ll<-max.ll<-c.ll<-list()
  nflips=0
  niters=0
  if (verbose) cat("re-phasing... ")
  if (t.HPC!=2)
  {
    for (ch in 1:t.nchrno) 
    {
      if (t.HPC==1)
      {
	donates_chr=getdonates(t.donates[[ch]],t.NUMI)
	donatesl_chr=getdonates(t.donatesl[[ch]],t.NUMI)
	donatesr_chr=getdonates(t.donatesr[[ch]],t.NUMI)
	orig.ll[[ch]]<-max.ll[[ch]]<-c.ll[[ch]]<-list()
	tmp<-foreach(ind=1:t.NUMI) %dopar%
	  phase_hunt(t.eps.lower,ch,t.G,t.A,t.GpcM,ind,t.flips[[ind]][[ch]], FALSE, t.ndonors[[ch]][[ind]], t.NUMP, t.NUMA, t.kLL, t.max.donors, donates_chr[[ind]], donatesl_chr[[ind]], donatesr_chr[[ind]], 
		     t.transitions[[ind]], t.umatch[[ch]], t.maxmatchsize[ch], t.d.w[[ch]], t.t.w[[ch]], t.gobs[[ch]][[ind]], t.mutmat, t.maxmiss, 
		     t.initProb, t.label, minbg=t.min.bg, maxbg=t.max.bg)
      }
      if (!t.HPC)
      {
	orig.ll[[ch]]<-max.ll[[ch]]<-c.ll[[ch]]<-list()
	tmp<-foreach(ind=1:t.NUMI) %dopar%
	  phase_hunt(t.eps.lower,ch,t.G,t.A,t.GpcM,ind,t.flips[[ind]][[ch]], FALSE, t.ndonors[[ch]][[ind]], t.NUMP, t.NUMA, t.kLL, t.max.donors, t.donates[[ch]][[ind]], t.donatesl[[ch]][[ind]], 
		     t.donatesr[[ch]][[ind]], t.transitions[[ind]], t.umatch[[ch]], t.maxmatchsize[ch], t.d.w[[ch]], t.t.w[[ch]], t.gobs[[ch]][[ind]], t.mutmat, t.maxmiss, 
		     t.initProb, t.label, minbg=t.min.bg, maxbg=t.max.bg)
      }
      for (ind in 1:(t.NUMI)) 
      {
	t.flips[[ind]][[ch]]<-tmp[[ind]]$ind.max.flips
	max.ll[[ch]][[ind]]<-c.ll[[ch]][[ind]]<-tmp[[ind]]$ind.max.ll
	orig.ll[[ch]][[ind]]<-tmp[[ind]]$ind.orig.ll
	nflips=nflips+tmp[[ind]]$nflips
	niters=niters+tmp[[ind]]$niters
      }
    } 
  }
  if (t.HPC==2)
  {
    tmp<-foreach(ch_ind=1:(t.nchrno*t.NUMI)) %dopar%
    {
      ch=as.integer((ch_ind-0.5)/t.NUMI)+1
      ind=(ch_ind-1)%%t.NUMI+1
      donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
      donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
      donatesr_chr_ind=getdonates_ind(t.donatesr[[ch]][[ind]])
      ans=phase_hunt(t.eps.lower,ch,t.G,t.A,t.GpcM,ind,t.flips[[ind]][[ch]], FALSE, t.ndonors[[ch]][[ind]], t.NUMP, t.NUMA, t.kLL, t.max.donors, donates_chr_ind, donatesl_chr_ind, donatesr_chr_ind, 
		     t.transitions[[ind]], t.umatch[[ch]], t.maxmatchsize[ch], t.d.w[[ch]], t.t.w[[ch]], t.gobs[[ch]][[ind]], t.mutmat, t.maxmiss, 
		     t.initProb, t.label, minbg=t.min.bg, maxbg=t.max.bg)
      ans
    }
      for (ch in 1:t.nchrno)
      {
	orig.ll[[ch]]<-max.ll[[ch]]<-c.ll[[ch]]<-list()
	for (ind in 1:t.NUMI)
	{
	  ch_ind=(ch-1)*t.NUMI+ind
	  t.flips[[ind]][[ch]]<-tmp[[ch_ind]]$ind.max.flips
	  max.ll[[ch]][[ind]]<-c.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.max.ll
	  orig.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.orig.ll
	  nflips=nflips+tmp[[ch_ind]]$nflips
	  niters=niters+tmp[[ch_ind]]$niters
	}
      }
  }
  rm(tmp)
  cloglike<-sum(unlist(max.ll))
  if (verbose)
    cat(nflips, " phase flips made after an average of ", niters/t.NUMI/t.nchrno, "hunts/ind/chromosome: log-likelihood", sum(unlist(orig.ll)), "-> ", cloglike, "\n")
  runtime<-as.numeric(Sys.time());diff.time<-runtime-t.old.runtime;t.old.runtime<-runtime;
  if (t.LOG) 
    writelog(t.logfile,"phasehunt",diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,cloglike) 
  ans=list()
  ans$flips=t.flips
  ans$runtime=runtime
  ans$cloglike=cloglike
  return(ans)
}


# function to run iterations of the MCMC algorithm for re-phasing target genomes based on MOSAIC fit. 
# phase-hunter is far more efficient and finds approximate best phasing in practice
phase_mcmc_all=function(t.donates, t.donatesl, t.donatesr, t.ndonors, t.NUMP, t.NUMA, t.G, t.A, t.max.donors, t.nchrno, t.NUMI, t.kLL, t.flips, 
			t.transitions, t.umatch, t.maxmatchsize, t.d.w, t.t.w, t.gobs, t.mutmat, t.maxmiss, 
			t.initProb, t.label,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,t.old.runtime, 
			t.logfile, t.HPC=2, verbose=TRUE, t.LOG=TRUE) {
  max.ll<-mcmc.ll<-mcmc.acc<-list()
  for (ch in 1:t.nchrno) 
  {
    max.ll[[ch]]<-mcmc.ll[[ch]]<-mcmc.acc[[ch]]<-list()
  }
  if (t.HPC!=2)
  {
    for (ch in 1:t.nchrno)
    {
      if (t.HPC==1)
      {  
	donates_chr=getdonates(t.donates[[ch]],t.NUMI)
	donatesl_chr=getdonates(t.donatesl[[ch]],t.NUMI)
	donatesr_chr=getdonates(t.donatesr[[ch]],t.NUMI)
	tmp<-foreach(ind=1:t.NUMI) %dopar%
	  phase_mcmc(ch,t.G,t.A,ind,M,t.max.donors,t.initProb,t.label, t.flips[[ind]][[ch]], verbose,t.ndonors[[ch]][[ind]], t.NUMP, t.NUMA, t.kLL, t.max.donors, 
		     donates_chr[[ind]], donatesl_chr[[ind]], donatesr_chr[[ind]], t.transitions[[ind]], t.umatch[[ch]], t.maxmatchsize[ch], t.d.w[[ch]], t.t.w[[ch]], t.gobs[[ch]][[ind]], t.mutmat, t.maxmiss, mcmcprog)
      }
      if (!t.HPC)
      {
	tmp<-foreach(ind=1:t.NUMI) %dopar%
	  phase_mcmc(ch,t.G,t.A,ind,M,t.max.donors,t.initProb,t.label, t.flips[[ind]][[ch]], verbose,t.ndonors[[ch]][[ind]], t.NUMP, t.NUMA, t.kLL, t.max.donors, 
		     t.donates[[ch]][[ind]], t.donatesl[[ch]][[ind]], t.donatesr[[ch]][[ind]], t.transitions[[ind]], t.umatch[[ch]], t.maxmatchsize[ch], t.d.w[[ch]], t.t.w[[ch]], t.gobs[[ch]][[ind]], t.mutmat, t.maxmiss,  mcmcprog)
      }
      for (ind in 1:t.NUMI)
      {
	if (length(tmp[[ind]]$ind.mcmc.acc)>0) 
	{
	  t.flips[[ind]][[ch]]<-tmp[[ind]]$ind.max.flips
	  mcmc.ll[[ch]][[ind]]<-tmp[[ind]]$ind.mcmc.ll
	  mcmc.acc[[ch]][[ind]]<-tmp[[ind]]$ind.mcmc.acc
	  max.ll[[ch]][[ind]]<-tmp[[ind]]$ind.max.ll
	}
      }
    }
  }
  if (t.HPC==2)
  {  
    tmp<-foreach(ch_ind=1:(t.nchrno*t.NUMI)) %dopar%
    {
      ch=as.integer((ch_ind-0.5)/t.NUMI)+1
      ind=(ch_ind-1)%%t.NUMI+1
      donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
      donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
      donatesr_chr_ind=getdonates_ind(t.donatesr[[ch]][[ind]])
      ans=phase_mcmc(ch,t.G,t.A,ind,M,t.max.donors,t.initProb,t.label,t.flips[[ind]][[ch]], verbose,t.ndonors[[ch]][[ind]], t.NUMP, t.NUMA, t.kLL, t.max.donors, 
		     donates_chr_ind, donatesl_chr_ind, donatesr_chr_ind, t.transitions[[ind]], t.umatch[[ch]], t.maxmatchsize[ch], t.d.w[[ch]], t.t.w[[ch]], t.gobs[[ch]][[ind]], t.mutmat, t.maxmiss, mcmcprog)
      ans
    }
      for (ch in 1:t.nchrno)
      {
	for (ind in 1:t.NUMI)
	{
	  ch_ind=(ch-1)*t.NUMI+ind
	  if (length(tmp[[ch_ind]]$ind.mcmc.acc)>0) 
	  {
	    t.flips[[ind]][[ch]]<-tmp[[ch_ind]]$ind.max.flips
	    mcmc.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.mcmc.ll
	    mcmc.acc[[ch]][[ind]]<-tmp[[ch_ind]]$ind.mcmc.acc
	    max.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.max.ll
	  }
	}
      }
  }
  rm(tmp)
  cloglike<-sum(unlist((max.ll)))
  runtime<-as.numeric(Sys.time());diff.time<-runtime-t.old.runtime;t.old.runtime<-runtime;
  if (t.LOG) 
    writelog(t.logfile,"phasemcmc",diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,cloglike) 
  ans=list()
  ans$flips=t.flips
  ans$runtime=runtime
  ans$cloglike=cloglike
  return(ans)
}
