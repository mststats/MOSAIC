# functions used to perform EM inference to initialise copying matrix Mu. First chromopainter style copying is inferred then windowed 
# (0.5cM by default). Then calculate the expected number of switches into each donor group in each window. 
# Finally, fit a mixture model (using EMmult.R) where the number of mixtures is the number of hidden ancestries we wish to model.
window_chunks<-function(nswitches, t.dr, t.G, t.kLL, t.NUMA, ww=0.5, min.swiches=0e-4, verbose=FALSE) # windows of 0.5M sum(G)*dr/w*100~=7000 for w=0.5
{
  nchrno=length(nswitches)
  #for (ch in 1:nchrno) nswitches[[ch]]<-nswitches[[ch]][-ignorepanels,,]
  # convert into window lengths in terms of #gridpoints
  #nw=sum(G)*dr/w*100; w=sum(G)/nw 
  w=ww/t.dr/100
  if (verbose) cat("Initialising copying matrix Mu;", w, "gridpoints per", ww, "cM width window\n") 
  nw=as.integer(t.G/w);offsetw=cumsum(nw)-nw
  wmat<-list()
  for (h in 1:t.NUMA)
  {
    wmat[[h]]=matrix(0,t.kLL,sum(nw)) # collate switching across chromosomes and haplotypes
    for (ch in 1:nchrno)
    {
      nswitches[[ch]][nswitches[[ch]]<min.swiches]=NaN
      for (i in 1:nw[ch])
      {
	wind=((i-1)*w+1):(i*w)
	wmat[[h]][,i+offsetw[ch]]=rowSums(nswitches[[ch]][,wind,h],na.rm=T) # sum over window 
      }
      if ((nw[ch]*w)<t.G[ch]) # if any left over just add in to last one
      {
	wind<-(nw[ch]*w+1):t.G[ch]
	wmat[[h]][,i+offsetw[ch]]=wmat[[h]][,i+offsetw[[ch]]]+rowSums(as.matrix(nswitches[[ch]][,wind,h]),na.rm=T) # sum over window 
      }
    }
  }
  return(list(wmat=wmat,w=w))
}

# Fit a mixture model to windowed switching counts of a chromopainter style model fit 
# where the number of mixtures is the number of hidden ancestries we wish to model.
tmpexp=exp(.Machine$double.min.exp/1.4)
zhclust=function(veccounts,k)
{
  hdists<-dist(veccounts)
  clustering=cutree(hclust(hdists,method="ward.D2"), k) # small spherical clusters
  norm=rowMeans(veccounts)
  n=length(clustering)
  groups=sort(unique(clustering))
  groups<-as.numeric(factor(groups))
  clustering<-as.numeric(factor(as.character(clustering)))
  Z=matrix(0, n, k)
  for (j in 1:k) 
    Z[clustering == groups[j], j]=1
  # Z is now the hard clustering matrix; soften to get membership
  membership=Z*0
  for(i in 1:k)
    for (l in 1:n)
    {
      membership[l,i]=1/mean((veccounts[l,]-colMeans(as.matrix(veccounts[Z[,i]==1,],sum(Z[,i]))))^2)
    }
  return(list(membership=membership))
}
ziffy<-function(veccounts, k) 
{
  if (k<as.integer(nrow(veccounts)/2)) 
  {
    tmp=fanny(veccounts, k, memb.exp=1.05);Mu=tmp$membership
  }
  if (k>=as.integer(nrow(veccounts)/2)) 
  {
    tmp=zhclust(veccounts, k);Mu=tmp$membership
  }
  Mu[which(is.nan(Mu))]=1/nrow(Mu) # remove NaNs
  Mu[which(is.infinite(Mu))]=100 # remove infinities
  Mu=t(t(Mu)/colSums(Mu))
  alpha=colSums(Mu)
  alpha=alpha/sum(alpha)
  return(list(Mu=Mu,alpha=alpha))
}

f=function(countscol, t.alpha, logMu, t.kLL, t.A) 
{
  log(sum(t.alpha*exp(.colSums(countscol*logMu,t.kLL,t.A))))
}
loglikeMult<-function(counts, t.alpha, logMu, t.kLL, t.A, t.NUMI)
{
  ans<-foreach(ind=1:t.NUMI) %dopar%  
  {
    hap=c(ind*2-1,ind*2)
    sum(apply(counts[[hap[1]]], 2, f, t.alpha[[ind]], logMu, t.kLL, t.A))+sum(apply(counts[[hap[2]]], 2, f, t.alpha[[ind]], logMu, t.kLL, t.A))
  }
    return(sum(unlist(ans)))
}

r_EMmult<-function(counts, t.A, t.NUMI, t.NUMA, itmax=200, eps=log(1.01), verbose=FALSE, t.singlePI=FALSE) # # i.e. a 1% increase in relative likelihood
{
  t.kLL=nrow(counts[[1]])
  nw=ncol(counts[[1]])
  if (verbose) cat("Fitting mixture model of switch counts in windows\n")
  # initial M-step
  tmp=matrix(0,t.kLL,nw*t.NUMA)
  for (h in 1:t.NUMA)
    tmp[,((h-1)*nw+1):(h*nw)]=counts[[h]] # stack across target haps
  tmp=tmp+1e-3
  tmp=qlogis(tmp/rowSums(tmp))
  tmp<-ziffy(tmp,t.A)
  alpha=list()
  for (h in 1:t.NUMI) 
    alpha[[h]]=tmp$alpha
  Mu=tmp$Mu
  p<-array(NaN,c(nw,t.NUMA,t.A))
  ll=NaN
  #return(list(Mu=Mu,alpha=alpha,ll=ll,p=p)) # for debugging
  for (iter in 1:itmax)
  {
    old.ll=ll
    # E-step
    ans<-foreach(h=1:t.NUMA) %dopar%  
    {
      hp<-matrix(NaN,nw,t.A)
      for (i in 1:t.A) 
	hp[,i]<-log(alpha[[as.integer((h+1)/2)]][i])+.colSums(counts[[h]]*log(Mu[,i]),t.kLL,nw) # prob target is of anc i in window iw
      hp<-hp-apply(hp,1,max);hp<-exp(hp); # yes, it's 1
      hp<-hp/rowSums(hp)
      hp
    }
      for (h in 1:t.NUMA)
	p[,h,]=ans[[h]]
    # M-step
    o.Mu<-Mu
    Mu[]=0
    for (i in 1:t.A) {
      for (k in 1:t.kLL)
	for (h in 1:t.NUMA)
	  Mu[k,i]=Mu[k,i]+sum(p[,h,i]*counts[[h]][k,]) # mean of Multinomial = np
        Mu[which(is.nan(Mu))]=1/t.kLL # remove NaNs
      if (all(Mu[,i]==0)) Mu[,i]=1/t.kLL # if all zero replace with 1/t.kLL
    }
    Mu<-t(t(Mu)/colSums(Mu))
    Mu[Mu<tmpexp]=tmpexp;Mu<-t(t(Mu)/colSums(Mu))
    if (!t.singlePI)
      for (ind in 1:t.NUMI)
      {
	hap=c(ind*2-1,ind*2)
	alpha[[ind]]<-.colSums(p[,hap[1],]+p[,hap[2],], nw, t.A)
      }
    if (t.singlePI)
    {
      tmpalpha=rep(0,t.A)
      for (ind in 1:t.NUMI)
      {
	hap=c(ind*2-1,ind*2)
	tmpalpha=tmpalpha+.colSums(p[,hap[1],]+p[,hap[2],], nw, t.A)
      }
      for (ind in 1:t.NUMI)
	alpha[[ind]]<-tmpalpha
    }
    for (ind in 1:t.NUMI)
      alpha[[ind]]<-alpha[[ind]]/sum(alpha[[ind]])
    logMu=log(Mu)
    print(alpha)
    tmp.ll=loglikeMult(counts, alpha, logMu, t.kLL, t.A, t.NUMI)
    ll=tmp.ll
    if (verbose)
      cat(iter,"/",itmax, " ", ll, " ", ll-old.ll, "\n", sep="")
    if (!is.nan(old.ll))
      if ((ll-old.ll)<eps)
      {
	if (ll<old.ll-eps) warning("########## Decreasing log-likelihood ##########",immediate.=T)
	if (verbose) cat("EM converged in mixture model for initialising copying matrix Mu\n")
	ll=old.ll # shouldn't see LL decreases but...
	break
      }
  }
  return(list(Mu=Mu,alpha=alpha,ll=ll,p=p))
}
EMmult<-r_EMmult
#EMmult<-cmpfun(r_EMmult,list(optimize=3))

cluster_windows<-function(windows,t.dr,t.kLL,t.A,t.NUMI,t.NUMA,t.NL,t.absorbrho,verbose=F)
{
  nw=ncol(windows$wmat[[1]])
  res<-EMmult(windows$wmat, t.A, t.NUMI, t.NUMA, verbose=verbose) # fit the mixture model
  Mu<-res$Mu*t.NL[1:t.kLL];Mu<-t(t(Mu)/colSums(Mu)) # note the use of NL to scale by panel size
  alpha<-res$alpha
  gpw=windows$w # gridpoints per window
  lambda=list()
  PI<-list()
  for (ind in 1:t.NUMI) 
  {
    hap=c(ind,ind);if (t.NUMA>1) hap=c(ind*2-1,ind*2) # haploid or diploid
    PI[[ind]]=matrix(0,t.A,t.A)
    for (i in 1:t.A)
    {
      for (j in 1:t.A)
      {
	# average p(i->j)/p(i) across windows and then average this over both haplotypes
	PI[[ind]][i,j]=0.5*(sum(res$p[-nw,hap[1],i]*res$p[-1,hap[1],j])/sum(res$p[-nw,hap[1],i])+
			    sum(res$p[-nw,hap[2],i]*res$p[-1,hap[2],j])/sum(res$p[-nw,hap[2],i]))
      }
      if (t.absorbrho) PI[[ind]][i,i]=0
    }
    PI[[ind]]=1-(1-PI[[ind]])^(1/gpw) # rescale PI to grid widths from window widths
    PI[[ind]]=PI[[ind]]/5 # ad-hoc reduction as this is typically overestimated
    tmp=which(is.na(PI[[ind]]),arr.ind=T) 
    for (i in tmp[,1]) for (j in tmp[,2]) PI[[ind]][i,j]=alpha[[ind]][j] # if never in an anc, just randomly choose another one w.p. alpha
    trans=PI[[ind]]-diag(rowSums(PI[[ind]])) # strictly speaking this should include an initial condition
    w=prcomp(t(trans),center=F)$rotation[,t.A]
    alpha[[ind]]=w/sum(w);alpha[[ind]][alpha[[ind]]<0]=0 # very close to the above alpha but now consistent with PI estimates
    lambda[[ind]]=-log(1-PI[[ind]])/t.dr
  }
  return(list(Mu=Mu,alpha=alpha,lambda=lambda,PI=PI,ll=res$ll)) 
}
