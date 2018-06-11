# Fit a mixture model to windowed switching counts of a chromopainter style model fit 
# where the number of mixtures is the number of hidden ancestries we wish to model.
require(gtools)
require(doParallel)
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
  if (k<(nrow(veccounts)/2)) 
  {
    require(cluster);tmp=fanny(veccounts, k, memb.exp=1.05);Mu=tmp$membership;alpha=colSums(Mu)
  }
  if (k>=(nrow(veccounts)/2)) 
  {
    tmp=zhclust(veccounts, k);Mu=tmp$membership;alpha=colSums(Mu)
  }
  Mu=t(t(Mu)/colSums(Mu))
  alpha=alpha/sum(alpha)
  return(list(Mu=Mu,alpha=alpha))
}

f=function(countscol, t.alpha, logMu, t.kLL, t.L) 
{
  log(sum(t.alpha*exp(.colSums(countscol*logMu,t.kLL,t.L))))
}
loglikeMult<-function(counts, t.alpha, logMu, t.kLL, t.L)
{
  ans<-foreach(ind=1:NUMI) %dopar%  
  {
    hap=c(ind*2-1,ind*2)
    sum(apply(counts[[hap[1]]], 2, f, t.alpha[[ind]], logMu, t.kLL, t.L))+sum(apply(counts[[hap[2]]], 2, f, t.alpha[[ind]], logMu, t.kLL, t.L))
  }
    return(sum(unlist(ans)))
}

r_EMmult<-function(counts, t.L, itmax=200, eps=log(1.01), verbose=F) # # i.e. a 1% increase in relative likelihood
{
  kLL=nrow(counts[[1]])
  nw=ncol(counts[[1]])
  if (verbose) cat("Fitting mixture model of switch counts in windows\n")
  # initial M-step
  tmp=matrix(0,kLL,nw*NUMA)
  for (h in 1:NUMA)
    tmp[,((h-1)*nw+1):(h*nw)]=counts[[h]] # stack across target haps
  tmp=tmp+1e-3
  tmp=qlogis(tmp/rowSums(tmp))
  tmp<-ziffy(tmp,t.L)
  alpha=list()
  for (h in 1:NUMI) 
    alpha[[h]]=tmp$alpha
  Mu=tmp$Mu
  p<-array(NaN,c(nw,NUMA,t.L))
  ll=NaN
  #return(list(Mu=Mu,alpha=alpha,ll=ll,p=p)) # for debugging
  for (iter in 1:itmax)
  {
    old.ll=ll
    # E-step
    ans<-foreach(h=1:NUMA) %dopar%  
    {
      hp<-matrix(NaN,nw,t.L)
      for (i in 1:t.L) 
	hp[,i]<-log(alpha[[as.integer((h+1)/2)]][i])+.colSums(counts[[h]]*log(Mu[,i]),kLL,nw) # prob target is of anc i in window iw
      hp<-hp-apply(hp,1,max);hp<-exp(hp); # yes, it's 1
      hp<-hp/rowSums(hp)
      hp
    }
      for (h in 1:NUMA)
	p[,h,]=ans[[h]]
    # M-step
    o.Mu<-Mu
    Mu[]=0
    for (k in 1:kLL)
      for (i in 1:t.L)
	for (h in 1:NUMA)
	  Mu[k,i]=Mu[k,i]+sum(p[,h,i]*counts[[h]][k,]) # mean of Multinomial = np
    Mu<-t(t(Mu)/colSums(Mu))
    Mu[Mu<tmpexp]=tmpexp;Mu<-t(t(Mu)/colSums(Mu))
    if (!singlePI)
      for (ind in 1:NUMI)
      {
	hap=c(ind*2-1,ind*2)
	alpha[[ind]]<-.colSums(p[,hap[1],]+p[,hap[2],], nw, t.L)
      }
    if (singlePI)
    {
      tmpalpha=rep(0,t.L)
      for (ind in 1:NUMI)
      {
	hap=c(ind*2-1,ind*2)
	tmpalpha=tmpalpha+.colSums(p[,hap[1],]+p[,hap[2],], nw, t.L)
      }
      for (ind in 1:NUMI)
	alpha[[ind]]<-tmpalpha
    }
    for (ind in 1:NUMI)
      alpha[[ind]]<-alpha[[ind]]/sum(alpha[[ind]])
    logMu=log(Mu)
    tmp.ll=loglikeMult(counts, alpha, logMu, kLL, t.L)
    ll=tmp.ll
    if (verbose)
      cat(iter,"/",itmax, " ", ll, " ", ll-old.ll, "\n", sep="")
    if (!is.nan(old.ll))
      if ((ll-old.ll)<eps)
      {
	if (ll<old.ll-eps) warning("Decreasing log-likelihood",immediate.=T)
	if (verbose) cat("EM converged in mixture model for initialising copying matrix Mu\n")
	ll=old.ll # shouldn't see LL decreases but...
	break
      }
  }
  return(list(Mu=Mu,alpha=alpha,ll=ll,p=p))
}
EMmult<-r_EMmult
#EMmult<-cmpfun(r_EMmult,list(optimize=optlevel))

