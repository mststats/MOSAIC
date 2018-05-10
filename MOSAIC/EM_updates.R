# function that performs the EM updates for MOSAIC's parameters. First some statistics that are used in multiple EM updates are calculated
# using functions in intermediate_calcs.R using Cpp code. These amount to counts of various switch types along the target genomes
update_params=function(t.HPC, t.nchrno, t.donates, t.donatesl, t.donatesr, t.NUMA, t.L, t.max.donors, t.NUMP, t.NUMI, t.G, t.transitions, t.flips,t.umatch, t.maxmatchsize,
		       t.d.w,t.t.w,t.gobs,t.mutmat,t.maxmiss,t.kLL, t.PI, t.alpha, t.lambda, t.Mu, t.rho, t.theta, t.ndonors, t.doPI, t.dorho, t.dotheta, t.doMu)
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
	  calc_E.n(ch,k,t.max.donors,t.NUMP,t.NUMA,t.G[ch],t.transitions[[ind]],t.flips[[ind]][[ch]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.kLL,t.L,t.PI,t.rho,t.Mu,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],donatesr_chr[[ind]])
	}
      }
      if (!t.HPC)
      {
	E.n[[ch]]<-foreach(k=1:t.NUMA) %dopar% 
	{
	  ind<-as.integer((k+1)*0.5)
	  calc_E.n(ch,k,t.max.donors,t.NUMP,t.NUMA,t.G[ch],t.transitions[[ind]],t.flips[[ind]][[ch]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.kLL,t.L,t.PI,t.rho,t.Mu,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesl[[ch]][[ind]],t.donatesr[[ch]][[ind]])
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
      ans=calc_E.n(ch,k,t.max.donors,t.NUMP,t.NUMA,t.G[ch],t.transitions[[ind]],t.flips[[ind]][[ch]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],
		   t.gobs[[ch]][[ind]],t.mutmat,t.maxmiss,t.kLL,t.L,t.PI,t.rho,t.Mu,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
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
    initi<-array(NaN,c(t.nchrno,t.NUMA,t.L,t.NUMP)) # a-posteriori first probs
    for (k in 1:t.NUMA)
    {
      for (ch in 1:t.nchrno) 
	for (k in 1:t.NUMA) {
	  initi[ch,k,,]<-E.n[[ch]][[k]]$initi
	}
    }
    initg=array(NaN, c(t.nchrno,t.NUMA,t.L,t.kLL));
    for (k in 1:t.kLL) 
      #initg[,,,k]=apply(initi[,,,which(label==k)],-4,sum)
      initg[,,,k]=apply(array(initi[,,,which(label==k)],c(t.nchrno,t.NUMA,t.L,sum(label==k))),-4,sum)
    t.Mu[]<-0
    for (ja in 1:t.L)
    {
      for (jl in 1:t.kLL)
	for (ch in 1:t.nchrno)
	  for (k in 1:t.NUMA)
	    t.Mu[jl,ja]<-t.Mu[jl,ja]+initg[ch,k,ja,jl]+sum(E.n[[ch]][[k]]$a[,jl,ja]) + E.n[[ch]][[k]]$r[jl,ja]
      if (all(t.Mu[,ja]==0)) t.Mu[,ja]=1/t.kLL
    }
    t.Mu<-t(t(t.Mu)/colSums(t.Mu))
  }
  if (t.doPI) # note that unlike the other parameters, each individual gets their own t.PI matrix 
  {
    # note that E.n[[k]]$l[i] = sum(E.n[[k]]$na[,i])+sum(E.n[[k]]$a[,,i])
    # note that sum(E.n[[k]]$na[,i]) = sum(E.n[[k]]$n[,i])+sum(E.n[[k]]$r[,i])
    if (!exists("singlePI")) singlePI=F
    for (ind in 1:t.NUMI)
    {
      if (t.NUMA>1) {hap<-c(ind*2-1,ind*2)} else hap=1
      if (singlePI) hap=1:t.NUMA # use all to compute each; i.e. single t.PI, etc
      for (i in 1:t.L){
	denom=0
	for (k in hap) 
	  for (ch in 1:t.nchrno)
	    denom=denom+sum(E.n[[ch]][[k]]$na[,i])+sum(E.n[[ch]][[k]]$a[i,,])
	for (j in 1:t.L) {
	  numer=0
	  for (k in hap) 
	    for (ch in 1:t.nchrno)
	      numer=numer+sum(E.n[[ch]][[k]]$a[i,,j])
	  t.PI[[ind]][i,j]=sum(numer)/sum(denom) 
	}
	if (absorbrho) t.PI[[ind]][i,i]=0 # absorb into t.rho
      }
      tmp=which(is.na(t.PI[[ind]]),arr.ind=T) 
      for (i in tmp[,1]) for (j in tmp[,2]) t.PI[[ind]][i,j]=t.alpha[[ind]][j] # if never in an anc, just randomly choose another one w.p. t.alpha
      trans=t.PI[[ind]]-diag(rowSums(t.PI[[ind]])) # strictly speaking this should include an initial condition
      w=prcomp(t(trans),center=F)$rotation[,t.L]
      t.alpha[[ind]]=w/sum(w);t.alpha[[ind]][t.alpha[[ind]]<0]=0
      #t.lambda[[ind]]=-log(1-sum(t.PI[[ind]])+sum(diag(t.PI[[ind]])))/dr 
      tmp=1-t(t(t.PI[[ind]]/t.alpha[[ind]]));tmp[tmp==0]=NaN;tmp[tmp==1]=NaN;tmp[is.infinite(tmp)]=NaN;tmp[tmp<0]=NaN;diag(tmp)=NaN; # gives same off diagonals for t.L=2
      tmp=-log(tmp)/dr; # gives same off diagonals for t.L=2
      t.lambda[[ind]]=mean(tmp,na.rm=T)
      #t.lambda[[ind]]=mean(-log(1-t.PI[[ind]])/dr)
    }
  }
  if (t.dorho)
  {
    numer=denom=rep(0,t.L)
    for (ch in 1:t.nchrno)
      for (k in 1:t.NUMA)
      {
	numer=numer+colSums(E.n[[ch]][[k]]$r) # E.n[[ch]][[k]][i,,i] all zero now
	denom=denom+colSums(E.n[[ch]][[k]]$r)+colSums(E.n[[ch]][[k]]$n)+apply(E.n[[ch]][[k]]$a,1,sum)
      }
    if (!commonrho) {t.rho<-numer/denom;t.rho[denom==0]=mean(t.rho[denom>0])}
    if (commonrho) t.rho[]<-sum(numer)/sum(denom) 
    #cat("!commonrho = ", numer/denom, "commonrho = ", sum(numer)/sum(denom), "\n")
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
    if (!commontheta) {t.theta<-misses/(misses+matches);t.theta[(misses+matches)==0]=mean(t.theta[(misses+matches)>0])}
    if (commontheta) t.theta[]<-sum(misses)/sum(misses+matches) 
    mutmat<-fmutmat(t.theta, t.L, t.maxmiss, maxmatch)
  }
  rm(E.n) 
  if (t.dorho || t.doPI || t.doMu) 
  {
    for (ind in 1:t.NUMI)
      t.transitions[[ind]]<-s_trans(t.L,t.kLL,t.PI[[ind]],t.Mu,t.rho,NL)
    initProb=initprobs(T,t.NUMA,t.L,t.NUMP,t.kLL,t.PI,t.Mu,t.rho,t.alpha,label,NL)
  }
  return(list(PI=t.PI, alpha=t.alpha, lambda=t.lambda, Mu=t.Mu, rho=t.rho, theta=t.theta, transitions=t.transitions, mutmat=mutmat, initProb=initProb))
}


# run EM algorithm for total iterations
run_EM=function(verbose=F) {
  if (verbose) pb<-txtProgressBar(min=1,max=ITER,style=3)
  for (ITER in 1:total)
  {
    old.Mu<-Mu; old.PI<-PI; old.lambda<-lambda; old.alpha<-alpha; old.rho<-rho; old.theta<-theta
    old.mutmat=mutmat;old.transitions=transitions;old.initProb=initProb
    old.cloglike<-cloglike
    tmp=update_params(HPC, nchrno, donates, donatesl, donatesr, NUMA, L, max.donors, NUMP, NUMI, G, transitions, flips,umatch,maxmatchsize,d.w,t.w,gobs,mutmat,maxmiss,kLL,
		      PI, alpha, lambda, Mu, rho, theta, ndonors, doPI, dorho, dotheta, doMu)
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta
    transitions=tmp$transitions;mutmat=tmp$mutmat;oldinitProb=tmp$initProb
    #cat("###### check:", range(mutmat-old.mutmat), " : " , cloglike-old.cloglike,"########\n")
    # E-step: extra work here as fors will be calculated next iteration of E.n above
    cloglike=get_loglike(NUMA, nchrno, G, L, kLL, max.donors, NUMP, donates, donatesl, transitions, maxmatchsize, umatch, flips, mutmat, maxmiss, initProb)
    cat(round(100*ITER/total), "%: ", cloglike, "(", cloglike-old.cloglike, ")", "\n")
    if (!is.na(old.cloglike)) 
    {
      if ((old.cloglike - cloglike)>1e-3)
      {
	Mu<-old.Mu; PI<-old.PI; lambda<-old.lambda; alpha<-old.alpha; rho<-old.rho; theta<-old.theta
	transitions=old.transitions;mutmat=old.mutmat;initProb=old.initProb;
	cloglike<-old.cloglike
	warning("loglikelihood has decreased; abandoning EM", immediate.=T)
	break
      }
      if ((cloglike - old.cloglike)< eps) 
	{cat("EM iterations have converged\n");break;}
    }
    if (LOG) 
    {
      runtime<-as.numeric(Sys.time());diff.time<-runtime-old.runtime;old.runtime<-runtime;
      writelog(EMlogfile,"EM",diff.time,len)
    }
    if (verbose) setTxtProgressBar(pb, m)
  }
  if (verbose) close(pb)
  return(list(PI=PI,alpha=alpha,lambda=lambda,Mu=Mu,rho=rho,theta=theta,runtime=runtime,initProb=initProb,cloglike=cloglike,transitions=transitions,mutmat=mutmat))
}
