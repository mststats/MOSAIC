# function that performs the EM updates for MOSAIC's parameters. First some statistics that are used in multiple EM updates are calculated
# using functions in intermediate_calcs.R using Cpp code. These amount to counts of various switch types along the target genomes
update_params=function(t.HPC, t.nchrno, t.donates, t.donatesl, t.donatesr, t.NUMA, t.A, t.max.donors, t.NN, t.NUMP, t.NUMI, t.G, t.dr, t.transitions, 
		       t.flips,t.umatch, t.maxmatchsize, t.d.w,t.t.w,t.gobs,t.mutmat,t.maxmatch,t.maxmiss,t.kLL, t.PI, t.alpha, t.lambda, t.Mu, t.rho, t.theta, 
		       t.ndonors, t.doPI, t.dorho, t.dotheta, t.doMu, t.label, t.NL, t.initProbs, t.commonrho, t.commontheta, t.absorbrho, t.singlePI=FALSE)
{
  E.n<-list()
  if (t.HPC!=2)
  {
    for (ch in 1:t.nchrno)
    {
      if (t.HPC==1)
      {
	donates_chr=getdonates(t.donates[[ch]],t.NUMI)
	donatesl_chr=getdonates(t.donatesl[[ch]],t.NUMI)
	donatesr_chr=getdonates(t.donatesr[[ch]],t.NUMI)
	E.n[[ch]]<-foreach(k=1:t.NUMA) %dopar% 
	{
	  ind<-as.integer((k+1)*0.5)
	  calc_E.n(ch,k,t.max.donors,t.NN,t.NUMP,t.NL,t.NUMA,t.G[ch],t.transitions[[ind]],t.flips[[ind]][[ch]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.kLL,t.A,t.PI,t.rho,t.Mu,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],donatesr_chr[[ind]],t.initProbs,t.label,t.doMu)
	}
      }
      if (!t.HPC)
      {
	E.n[[ch]]<-foreach(k=1:t.NUMA) %dopar% 
	{
	  ind<-as.integer((k+1)*0.5)
	  calc_E.n(ch,k,t.max.donors,t.NN,t.NUMP,t.NL,t.NUMA,t.G[ch],t.transitions[[ind]],t.flips[[ind]][[ch]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.kLL,t.A,t.PI,t.rho,t.Mu,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesl[[ch]][[ind]],t.donatesr[[ch]][[ind]],t.initProbs,t.label,t.doMu)
	}
      }
    }
  }
  if (t.HPC==2)
  {
    tmp<-foreach(ch_k=(1:(t.NUMA*t.nchrno))) %dopar% 
    {
      ch=as.integer((ch_k-0.5)/t.NUMA)+1
      k=(ch_k-1)%%t.NUMA+1
      ind<-as.integer((k+1)*0.5)
      donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
      donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
      donatesr_chr_ind=getdonates_ind(t.donatesr[[ch]][[ind]])
      ans=calc_E.n(ch,k,t.max.donors,t.NN,t.NUMP,t.NL,t.NUMA,t.G[ch],t.transitions[[ind]],t.flips[[ind]][[ch]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],
		   t.gobs[[ch]][[ind]],t.mutmat,t.maxmiss,t.kLL,t.A,t.PI,t.rho,t.Mu,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind,t.initProbs,t.label,t.doMu)
      rm(donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
      ans
    }
      for (ch in 1:t.nchrno)
      {
	E.n[[ch]]=list()
	for (k in 1:t.NUMA)
	  E.n[[ch]][[k]]=tmp[[(ch-1)*t.NUMA+k]]
      }
  }
  cloglike=0;for(ch in 1:t.nchrno) for (k in 1:t.NUMA) cloglike=cloglike+E.n[[ch]][[k]]$loglike
  # parameter updates are very fast given E.n[[]]
  if (t.doMu)
  {
    initi<-array(NaN,c(t.nchrno,t.NUMA,t.A,t.NUMP)) # a-posteriori first probs
    for (k in 1:t.NUMA)
    {
      for (ch in 1:t.nchrno) 
	for (k in 1:t.NUMA) {
	  initi[ch,k,,]<-E.n[[ch]][[k]]$initi
	}
    }
    initg=array(NaN, c(t.nchrno,t.NUMA,t.A,t.kLL));
    for (k in 1:t.kLL) 
      #initg[,,,k]=apply(initi[,,,which(t.label==k)],-4,sum)
      initg[,,,k]=apply(array(initi[,,,which(t.label==k)],c(t.nchrno,t.NUMA,t.A,sum(t.label==k))),-4,sum)
    t.Mu[]<-0
    for (ja in 1:t.A)
    {
      for (jl in 1:t.kLL)
	for (ch in 1:t.nchrno)
	  for (k in 1:t.NUMA)
	    t.Mu[jl,ja]<-t.Mu[jl,ja]+initg[ch,k,ja,jl]+sum(E.n[[ch]][[k]]$a[,jl,ja]) + E.n[[ch]][[k]]$r[jl,ja]
      if (all(is.nan(t.Mu[,ja]))) t.Mu[,ja]=1/t.kLL
      if (all(t.Mu[,ja]==0)) t.Mu[,ja]=1/t.kLL
    }
    t.Mu<-t(t(t.Mu)/colSums(t.Mu))
  }
  if (t.doPI) # note that unlike the other parameters, each individual gets their own t.PI matrix 
  {
    # note that E.n[[k]]$l[i] = sum(E.n[[k]]$na[,i])+sum(E.n[[k]]$a[,,i])
    # note that sum(E.n[[k]]$na[,i]) = sum(E.n[[k]]$n[,i])+sum(E.n[[k]]$r[,i])
    for (ind in 1:t.NUMI)
    {
      if (t.NUMA>1) {hap<-c(ind*2-1,ind*2)} else hap=1
      if (t.singlePI) hap=1:t.NUMA # use all to compute each; i.e. single t.PI, etc
      for (i in 1:t.A){
	denom=0
	for (k in hap) 
	  for (ch in 1:t.nchrno)
	    denom=denom+sum(E.n[[ch]][[k]]$na[,i])+sum(E.n[[ch]][[k]]$a[i,,])
	for (j in 1:t.A) {
	  numer=0
	  for (k in hap) 
	    for (ch in 1:t.nchrno)
	      numer=numer+sum(E.n[[ch]][[k]]$a[i,,j])
	  t.PI[[ind]][i,j]=sum(numer)/sum(denom) 
	}
	if (t.absorbrho) t.PI[[ind]][i,i]=0 # absorb into t.rho
      }
      tmp=which(is.na(t.PI[[ind]]),arr.ind=T) 
      for (i in tmp[,1]) for (j in tmp[,2]) t.PI[[ind]][i,j]=t.alpha[[ind]][j] # if never in an anc, just randomly choose another one w.p. t.alpha
      trans=t.PI[[ind]]-diag(rowSums(t.PI[[ind]])) # strictly speaking this should include an initial condition
      w=prcomp(t(trans),center=F)$rotation[,t.A]
      t.alpha[[ind]]=w/sum(w);t.alpha[[ind]][t.alpha[[ind]]<0]=0
      #t.lambda[[ind]]=-log(1-sum(t.PI[[ind]])+sum(diag(t.PI[[ind]])))/dr 
      tmp=1-t(t(t.PI[[ind]]/t.alpha[[ind]]));tmp[tmp==0]=NaN;tmp[tmp==1]=NaN;tmp[is.infinite(tmp)]=NaN;tmp[tmp<0]=NaN;diag(tmp)=NaN; # gives same off diagonals for t.A=2
      tmp=-log(tmp)/t.dr; # gives same off diagonals for t.A=2
      t.lambda[[ind]]=mean(tmp,na.rm=T)
      #t.lambda[[ind]]=mean(-log(1-t.PI[[ind]])/dr)
    }
  }
  if (t.dorho)
  {
    numer=denom=rep(0,t.A)
    for (ch in 1:t.nchrno)
      for (k in 1:t.NUMA)
      {
	numer=numer+colSums(E.n[[ch]][[k]]$r) # E.n[[ch]][[k]][i,,i] all zero now
	denom=denom+colSums(E.n[[ch]][[k]]$r)+colSums(E.n[[ch]][[k]]$n)+apply(E.n[[ch]][[k]]$a,1,sum)
      }
    if (!t.commonrho) {t.rho<-numer/denom;t.rho[denom==0]=mean(t.rho[denom>0])}
    if (t.commonrho) t.rho[]<-sum(numer)/sum(denom) 
    #cat("!t.commonrho = ", numer/denom, "t.commonrho = ", sum(numer)/sum(denom), "\n")
  }
  if (t.dotheta)
  {
    misses=0;matches=0
    for (k in 1:t.NUMA)
    {
      for (ch in 1:t.nchrno) 
      {
	misses=misses+E.n[[ch]][[k]]$e
	matches=matches+E.n[[ch]][[k]]$h
      }
    }
    if (!t.commontheta) {t.theta<-misses/(misses+matches);t.theta[(misses+matches)==0]=mean(t.theta[(misses+matches)>0])}
    if (t.commontheta) t.theta[]<-sum(misses)/sum(misses+matches) 
    mutmat<-fmutmat(t.theta, t.A, t.maxmiss, t.maxmatch)
  } else mutmat=t.mutmat
  rm(E.n) 
  if (t.dorho || t.doPI || t.doMu) 
  {
    for (ind in 1:t.NUMI)
      t.transitions[[ind]]<-s_trans(t.A,t.kLL,t.PI[[ind]],t.Mu,t.rho,t.NL)
    initProb=initprobs(T,t.NUMA,t.A,t.NUMP,t.kLL,t.PI,t.Mu,t.rho,t.alpha,t.label,t.NL)
  }
  return(list(PI=t.PI, alpha=t.alpha, lambda=t.lambda, Mu=t.Mu, rho=t.rho, theta=t.theta, transitions=t.transitions, mutmat=mutmat, initProb=initProb))
}


# run EM algorithm for total iterations or until convergence
run_EM=function(t.HPC, t.nchrno, t.PI, t.Mu, t.rho, t.theta, t.alpha, t.lambda, t.initProb, t.label, t.mutmat, t.transitions, t.ndonors, t.donates, t.donatesl, 
		t.donatesr, t.NUMA, t.NN, t.NL, t.NUMP, t.kLL, t.A, t.NUMI, t.max.donors, t.G, t.dr, t.gobs, t.maxmatchsize, t.umatch, t.flips, t.maxmatch, t.maxmiss, 
		t.d.w, t.t.w,  t.total, verbose=F, t.len, t.cloglike, t.LOG, t.logfile, t.doPI, t.doMu, t.dotheta, t.dorho, t.commonrho, t.commontheta, t.absorbrho,
		t.old.runtime, t.eps) {
  if (verbose) pb<-txtProgressBar(min=1,max=ITER,style=3)
  runtime<-as.numeric(Sys.time());diff.time<-runtime-t.old.runtime # required in case of break below
  for (ITER in 1:t.total)
  {
    old.Mu<-t.Mu; old.PI<-t.PI; old.lambda<-t.lambda; old.alpha<-t.alpha; old.rho<-t.rho; old.theta<-t.theta
    old.mutmat=t.mutmat;old.transitions=t.transitions;old.initProb=t.initProb;old.cloglike<-t.cloglike
    tmp=update_params(t.HPC, t.nchrno, t.donates, t.donatesl, t.donatesr, t.NUMA, t.A, t.max.donors, t.NN, t.NUMP, t.NUMI, t.G, t.dr, t.transitions, t.flips,t.umatch,
		      t.maxmatchsize,t.d.w,t.t.w,t.gobs,t.mutmat,t.maxmatch,t.maxmiss,t.kLL,t.PI, t.alpha, t.lambda, t.Mu, t.rho, t.theta, t.ndonors, t.doPI, t.dorho, 
		      t.dotheta, t.doMu, t.label, t.NL, t.initProb, t.commonrho, t.commontheta, t.absorbrho)
    t.PI=tmp$PI;t.alpha=tmp$alpha;t.lambda=tmp$lambda;t.Mu=tmp$Mu;t.rho=tmp$rho;t.theta=tmp$theta
    t.transitions=tmp$transitions;t.mutmat=tmp$mutmat;t.initProb=tmp$initProb
    # E-step: extra work here as fors will be calculated next iteration of E.n above
    t.cloglike=get_loglike(t.NUMA, t.nchrno, t.G, t.A, t.kLL, t.max.donors, t.NUMP, t.ndonors, t.donates, t.donatesl, t.transitions, t.maxmatchsize, t.umatch, t.flips, t.mutmat, t.maxmiss, t.initProb,t.d.w,t.t.w,t.gobs,t.label, t.HPC)
    cat(round(100*ITER/t.total), "%: ", t.cloglike, "(", t.cloglike-old.cloglike, ")", "\n")
    if (!is.na(old.cloglike)) 
    {
      if ((old.cloglike - t.cloglike)>1e-3)
      {
	t.Mu<-old.Mu; t.PI<-old.PI; t.lambda<-old.lambda; t.alpha<-old.alpha; t.rho<-old.rho; t.theta<-old.theta
	t.transitions=old.transitions;t.mutmat=old.mutmat;t.initProb=old.initProb;t.cloglike<-old.cloglike
	warning("loglikelihood has decreased; abandoning EM", immediate.=T)
	break
      }
      if ((t.cloglike - old.cloglike)< t.eps) 
	{cat("EM iterations have converged\n");break;}
    }
    runtime<-as.numeric(Sys.time());diff.time<-runtime-t.old.runtime
    if (t.LOG) 
      writelog(t.logfile,"EM",diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,t.cloglike) 
    if (verbose) setTxtProgressBar(pb, m)
  }
  if (verbose) close(pb)
  return(list(PI=t.PI,alpha=t.alpha,lambda=t.lambda,Mu=t.Mu,rho=t.rho,theta=t.theta,runtime=runtime,initProb=t.initProb,cloglike=t.cloglike,transitions=t.transitions,mutmat=t.mutmat))
}
