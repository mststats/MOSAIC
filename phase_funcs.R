# functions used in re-phasing target haplotypes based on current MOSAIC fit
require(compiler)
require(parallel)
require(bit) 
flip.ll<-function(g, t.ch, ndonorsg, NNL2, t.scalefactor, t.scalefactorb, t.fors, t.backs)
{
  g=g-1
  tmpvec=1:(ndonorsg*L)
  ans=0
  for (h in 1:2)
  {
    ans=ans-sum(log(t.scalefactor[[h]][1:g]))
    ans=ans-sum(log(t.scalefactorb[[h]][g:G[t.ch]]))
    tmp=sum(t.fors[[h]][(g-1)*NNL2+tmpvec]*t.backs[[3-h]][(g-1)*NNL2+tmpvec])
    ans=ans+log(tmp)-log(NUMP)-log(L)
  }
  ans
}

r.create.proposal<-function(t.ch,max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll) 
{
  NNL2=max.donors*L
  g.flip.ll<-function(g) flip.ll(g, t.ch, t.ndonors[g], NNL2, t.scalefactor, t.scalefactorb, t.fors, t.backs)
  ll=c(ind.c.ll,vapply(2:G[t.ch], g.flip.ll,0))
  ll=ll-ind.c.ll
  ll[is.nan(ll)|is.infinite(ll)]=-1 # can happen if no path through
  ll
}
#create.proposal<-cmpfun(r.create.proposal,list(optimize=optlevel)) # 
create.proposal<-r.create.proposal
# H1 gets: fors[1] to g, backs[2] from g and new fors from g, new backs to g
# H2 gets: fors[2] to g, backs[1] from g and new fors from g, new backs to g
rephaser<-function(t.ch, hap, g, t.transitions, t.umatch, t.maxmatchsize, t.dw, t.tw, t.gobs, t.fors, t.sumfors, t.backs, t.scalefactor, t.scalefactorb, t.flips, t.initProb, t.ndonors, t.donates, t.donatesl, t.donatesr) # returns log-likelihood of data when g is flipped for jth individual
{
  c.fors<-c.sumfors<-c.backs<-c.scalefactor<-c.scalefactorb<-list()
  for (h in 1:2) 
  {
    c.fors[[h]]<-bit::clone(t.fors[[h]])
    c.sumfors[[h]]<-bit::clone(t.sumfors[[h]])
    c.backs[[h]]<-bit::clone(t.backs[[h]])
    c.scalefactor[[h]]<-bit::clone(t.scalefactor[[h]])
    c.scalefactorb[[h]]<-bit::clone(t.scalefactorb[[h]])
  }
  tmp.backs<-c.backs
  tmp.scalefactorb<-c.scalefactorb
  # swap the backward probabilities and scaling factors
  c.backs[[1]]<-tmp.backs[[2]];c.backs[[2]]<-tmp.backs[[1]]
  c.scalefactorb[[1]]<-tmp.scalefactorb[[2]];c.scalefactorb[[2]]<-tmp.scalefactorb[[1]]
  for (h in 1:2)
  {
    cppforward(hap[h],NUMA,max.donors,THIN,NUMP,kLL,L,g-1,G[t.ch],G[t.ch],t.transitions, t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,t.initProb[hap[h],],label,t.ndonors,t.donates,t.donatesl,t.flips,c.fors[[h]],c.sumfors[[h]],c.scalefactor[[h]]) # g-1
    cppbackward(hap[h],NUMA,max.donors,THIN,NUMP,L,0,g-1,G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,label,t.ndonors,t.donates,t.donatesr,t.flips,c.backs[[h]],c.scalefactorb[[h]]) # g-1
  }
  return(list(fors=c.fors, sumfors=c.sumfors, backs=c.backs, scalefactor=c.scalefactor, scalefactorb=c.scalefactorb))
}

r_phase_hunt<-function(t.eps.lower, t.ch, t.ind, t.flips, verbose, t.ndonors, t.donates, t.donatesl, t.donatesr, 
		       t.transitions, t.umatch, t.maxmatchsize, t.dw, t.tw, t.gobs, mutmat, maxmiss, t.phase.error.locs, t.initProb, lim=10, PLOT=F, 
		       minbg=0.1, maxbg=1, mult=1.5)
{
  hap<-c(t.ind*2-1,t.ind*2) # t.ind indexes over genotypes, hap the two haplotypes
  # fb calcs moved to here to avoid storing all fors, backs, etc
  t.fors<-list();t.sumfors<-list();t.scalefactor<-list()
  t.backs<-list();t.scalefactorb<-list()
  THIN=ifelse(max.donors==NUMP,F,T)
  for (h in 1:2)
  {
    k=hap[h]
    t.fors[[h]]<-rep(0,G[t.ch]*max.donors*L);t.sumfors[[h]]<-matrix(0,G[t.ch],L);t.scalefactor[[h]]<-rep(0,G[t.ch])
    cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,t.initProb[k,],label,
	       t.ndonors,t.donates,t.donatesl,t.flips,t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
    t.backs[[h]]<-rep(0,G[t.ch]*max.donors*L);t.scalefactorb[[h]]<-rep(0,G[t.ch])
    cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,label,
		t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
  }
  ind.orig.ll<-ind.max.ll<-ind.c.ll<- -sum(log(t.scalefactor[[1]]))-sum(log(t.scalefactor[[2]])) # LL for a single individual
  ind.c.v<-create.proposal(t.ch,max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll)
  c.probs<-1/(1+exp(-ind.c.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
  site.c.probs<-c.probs/sum(c.probs)
  nflips<-1
  t.nflips<-0
  iters=1
  end.eps.lower=t.eps.lower
  #if (PLOT) cat("\n")
  if (PLOT)
  {
    plot(g.loc[[t.ch]][c(1,G[[t.ch]])]*1e-6,c(1,lim-1),xlab="phase flip sites (Mb)",ylab="pass",t='n',axes=F)#,main="flipped positions")
    mp<-axTicks(1,round(axp=c(min(g.loc[[t.ch]])*1e-6,max(g.loc[[t.ch]])*1e-6,5)))
    axis(1,at=mp,labels=signif(mp,3))
    text(x=g.loc[[t.ch]][1]*1e-6,y=1:(lim-1),1:(lim-1))
  }
  if (max(ind.c.v)>t.eps.lower)
    while (nflips>0 & iters<lim) {
      #t.eps.lower=end.eps.lower+0.5*log(lim/iters) # fewer flips at the start in phase hunter, down to 2 as lower threshold for a phase flip
      bg=((2*plogis(mult*iters-mult))-1)*(maxbg-minbg)+minbg # smaller band at the start, increasing quickly
      BG=as.integer(bg*GpcM) # bg is number of centiMorgans to each side of a proposed phase flip around which don't try more flips
      if (max(ind.c.v)<t.eps.lower) 
      {
	if (PLOT) cat("Chromosome ", chrnos[t.ch], "individual", t.ind, "none left to flip ")
	break;
      }
      old.ll<-ind.c.ll
      cand.l<-(ind.c.v>t.eps.lower) # T for all those that would have log-like change greater than t.eps.lower if flipped
      cand.u<-which(cand.l) # which are candidates
      cand<-NULL # take the max spike, block off around it, repeat until none left
      while (max(ind.c.v[cand.u])>eps.lower & sum(cand.l)>0) { # while some are still worth flipping
	tmp<-which.max(ind.c.v[cand.u]) # identify maximum
	cand<-c(cand,cand.u[tmp]) # add to the vector of candidates
	cand.l[max(1,cand.u[tmp]-BG):min(cand.u[tmp]+BG,G[t.ch])]<-F # remove block of BG gridpoints to each side
	cand.u<-which(cand.l)
	if (length(cand.u)==0)
	  break
      }
      cand<-sort(cand,method="quick") # required?
      if (PLOT) 
      {
	points(g.loc[[t.ch]][cand]*1e-6, rep(iters,length(cand)),pch=20) 
      }
      #if (PLOT) for (i in 1:length(cand)) cat(cand[i], ":", ind.c.v[cand[i]], "\n") 
      nflips<-length(cand)
      t.nflips<-t.nflips+nflips
      if (nflips==0)
      {
	#if (PLOT) cat("\n")
	break;
      }
      for (g in cand) 
	t.flips[g:G[t.ch]]<-!t.flips[g:G[t.ch]]
      # full refit required here; potentially lots of re-phasing
      for (h in 1:2) {
	k<-hap[h]
	cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,t.initProb[k,],label,t.ndonors,t.donates,t.donatesl,t.flips,
		   t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
	cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,label,t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
      }
      ind.c.ll<- -sum(sum(log(t.scalefactor[[1]])))-sum(sum(log(t.scalefactor[[2]])))
      #if (PLOT) cat("Chromosome ", chrnos[t.ch], "individual", t.ind, " old:", old.ll, "with ", nflips, "flips:", ind.c.ll)  
      if (ind.c.ll==old.ll) # just break
      {
	#if (PLOT) cat("\n")
	break;
      }
      if (ind.c.ll<old.ll) { # only flip max one
	t.nflips<-t.nflips-nflips+1
	nflips=1
	for (g in cand) 
	  t.flips[g:G[t.ch]]<-!t.flips[g:G[t.ch]] # flip all candidates back
	g=which.max(ind.c.v)
	t.flips[g:G[t.ch]]<-!t.flips[g:G[t.ch]] # flip max candidate only
	# full refit required again here; can't simply re-use half of forwards and backwards as these have changed in attempting multiple flips
	for (h in 1:2) {
	  k<-hap[h]
	  t.ind<-as.integer((k+1)/2)
	  cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,t.initProb[k,],label,t.ndonors,t.donates,t.donatesl,t.flips,
		     t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
	  cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,label,t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
	}
	ind.c.ll<- -sum(sum(log(t.scalefactor[[1]])))-sum(sum(log(t.scalefactor[[2]])))
	if (is.nan(ind.c.ll)|is.infinite(ind.c.ll)) ind.c.ll=old.ll+ind.c.v[g];
	#if (PLOT) cat (" don't flip all these; only max one: ", ind.c.ll)
	#if (PLOT) cat(" largest flip only:", old.ll, old.ll+ind.c.v[g], ind.c.ll, "\n")
	#break; # actually we often find more again after a single flip
      } 
      #if (PLOT) cat("\n")
      # if pass above checks then loglikelihood has improved and we try again based on a new proposal
      ind.c.v<-create.proposal(t.ch,max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll)
      if (max(ind.c.v)<=0)
	break;
      if (verbose) cat("iter:", iters, "BG:", BG, "flips:", nflips, ind.c.ll, "\n")
      iters=iters+1
    }
  #if (PLOT) readline() 
  ind.max.ll<-ind.c.ll
  c.probs<-1/(1+exp(-ind.c.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
  site.c.probs<-c.probs/sum(c.probs) 
  if (verbose) cat(t.nflips, "phase flips made: log-likelihood ", ind.orig.ll, "->", ind.c.ll, " on Chr ", chrnos[t.ch], " for ind ", t.ind, 
		   "after ", iters, "hunting passes", "\n")
  # best so far so save to a file
  return(list(ind.max.ll=ind.max.ll,ind.orig.ll=ind.orig.ll,ind.max.flips=t.flips,nflips=t.nflips,niters=iters))
}

#phase_hunt<-cmpfun(r_phase_hunt,list(optimize=optlevel)) # 
phase_hunt<-r_phase_hunt
require(parallel)
r_phase_mcmc<-function(t.ch, t.ind, M, max.donors, t.initProb, t.flips, verbose, t.ndonors, t.donates, t.donatesl, t.donatesr, 
		       t.transitions, t.umatch, t.maxmatchsize, t.dw, t.tw, t.gobs, mutmat, maxmiss, PLOT=F, mcmcprog, mcmchill=T) 
{
  M=as.integer(M*G[t.ch])
  NNL2=max.donors*L
  if (mcmcprog) pb<-txtProgressBar(min=1,max=M,style=3)
  # fb calcs moved to here to avoid storing all fors, backs, etc
  t.fors<-list();t.sumfors<-list();t.scalefactor<-list()
  t.backs<-list();t.scalefactorb<-list()
  THIN=ifelse(max.donors==NUMP,F,T)
  for (h in 1:2)
  {
    k=hap[h]
    t.fors[[h]]<-rep(0,G[t.ch]*max.donors*L);t.sumfors[[h]]<-matrix(0,G[t.ch],L);t.scalefactor[[h]]<-rep(0,G[t.ch])
    cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,t.initProb[k,],label,
	       t.ndonors,t.donates,t.donatesl,t.flips,t.fors[[h]],t.sumfors[[h]],t.scalefactor[[h]])
    t.backs[[h]]<-rep(0,G[t.ch]*max.donors*L);t.scalefactorb[[h]]<-rep(0,G[t.ch])
    cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[t.ch],G[t.ch],t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,mutmat,maxmiss,label,
		t.ndonors,t.donates,t.donatesr,t.flips,t.backs[[h]],t.scalefactorb[[h]])
  }
  ind.orig.ll<-ind.max.ll<-ind.c.ll<- -sum(log(t.scalefactor[[1]]))-sum(log(t.scalefactor[[2]])) # LL for a single individual
  ind.max.flips<-t.flips
  ind.mcmc.ll<-ind.c.ll # running record of log-likelihood for this individual
  hap<-c(t.ind*2-1,t.ind*2) # t.ind indexes over genotypes, hap the two haplotypes
  ind.mcmc.acc<-0
  ind.c.v<-create.proposal(t.ch,max.donors,t.fors,t.sumfors,t.backs,t.scalefactor,t.scalefactorb,t.ndonors,ind.c.ll)
  c.probs<-1/(1+exp(-ind.c.v)) # ==exp(v)/(exp(v)+exp(p)) # marginal probabilities for each flip if we Gibbs sampled
  if (all(c.probs==0)) {cat("none left to propose at this hill climbing rate\n");return();}
  site.c.probs<-c.probs/sum(c.probs)
  if (verbose) cat("running MCMC chain...\n")
  p.flips<-t.flips
  for (m in 1:M) {
    #if (any(t.flips)) cat(m,1,"\n")
    ind.mcmc.acc<-c(ind.mcmc.acc,0) # default to zero, changes below if accepted
    #if (max(c.probs)==0) cat(range(ind.c.v), range(p.v), "\n") 
    g<-sample(2:G[t.ch],1,prob=c.probs[-1]) 
    if (PLOT) 
    {
      plot(g.loc[[t.ch]]*1e-6, ind.c.v,t='l')
      points(g.loc[[t.ch]][g]*1e-6, ind.c.v[g], col=2, lwd=2)
    }
    # given j, choose to propose to flip or not flip with probability c.probs[j]
    p.flip<-rbinom(1,p=c.probs[g],size=1) # this will be low. flip==1 => flip it
    #cat(ind.c.ll, mean(c.probs), c.probs[g], p.flip, "\n")
    if (p.flip==1) {
      # now update proposal distribution
      ##########################################################################################
      p.flips<-t.flips
      p.flips[g:G[t.ch]]<-!t.flips[g:G[t.ch]] # could have been flipped before in which case will get flipped back
      rephased<-rephaser(t.ch,hap,g,t.transitions,t.umatch,t.dw,t.tw,t.gobs,(t.fors),(t.sumfors),(t.backs),(t.scalefactor),(t.scalefactorb),(p.flips),t.initProb,t.ndonors,t.donates,t.donatesl,t.donatesr)
      #gc()
      n.ll<-flip.ll(g, t.ch, t.ndonors[g], NNL2, t.scalefactor, t.scalefactorb, t.fors, t.backs) # LL for a single individual
      p.v<-create.proposal(t.ch,max.donors,rephased$fors,rephased$sumfors,rephased$backs,rephased$scalefactor,rephased$scalefactorb,t.ndonors,n.ll)
      # below is super fast but incorrect; leads to bad proposals in the long run
      #p.v<-ind.c.v;p.v[g]=flip.ll(g, t.ch, t.ndonors[g], NNL2, rephased$scalefactor, rephased$scalefactorb, rephased$fors, rephased$backs)-n.ll
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
	t.flips[g:G[t.ch]]<-!t.flips[g:G[t.ch]] 
	for (h in 1:2) {
	  t.fors[[h]]<-bit::clone(rephased$fors[[h]])
	  t.sumfors[[h]]<-bit::clone(rephased$sumfors[[h]])
	  t.backs[[h]]<-bit::clone(rephased$backs[[h]])
	  t.scalefactor[[h]]<-bit::clone(rephased$scalefactor[[h]])
	  t.scalefactorb[[h]]<-bit::clone(rephased$scalefactorb[[h]])
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
  if (PLOT)
  {
    plot(ind.mcmc.ll, t='l')
    abline(h=ind.orig.ll,col=2)
    abline(h=ind.max.ll,col=3)
  }
  return(list(ind.max.flips=ind.max.flips,ind.mcmc.ll=ind.mcmc.ll,ind.max.ll=ind.max.ll,ind.mcmc.acc=ind.mcmc.acc))
}
#phase_mcmc<-cmpfun(r_phase_mcmc,list(optimize=optlevel)) # 
phase_mcmc<-r_phase_mcmc
