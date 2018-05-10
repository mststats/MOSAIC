# script that performs the EM updates for MOSAIC's parameters. First some statistics that are used in multiple EM updates are calculated
# using functions in intermediate_calcs.R using Cpp code. These amount to counts of various switch types along the target genomes
E.n<-list()
if (HPC!=2)
{
  for (ch in 1:nchrno)
  {
    if (HPC==1)
    {
      donates_chr=getdonates(donates[[ch]],NUMI)
      donatesl_chr=getdonates(donatesl[[ch]],NUMI)
      donatesr_chr=getdonates(donatesr[[ch]],NUMI)
      E.n[[ch]]<-foreach(k=1:NUMA) %dopar% 
      {
	ind<-as.integer((k+1)*0.5)
	calc_E.n(ch,k,max.donors,NUMP,NUMA,G[ch],transitions[[ind]],flips[[ind]][[ch]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],
		 mutmat,maxmiss,kLL,L,PI,rho,Mu,ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],donatesr_chr[[ind]])
      }
    }
    if (!HPC)
    {
      E.n[[ch]]<-foreach(k=1:NUMA) %dopar% 
      {
	ind<-as.integer((k+1)*0.5)
	calc_E.n(ch,k,max.donors,NUMP,NUMA,G[ch],transitions[[ind]],flips[[ind]][[ch]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],
		 mutmat,maxmiss,kLL,L,PI,rho,Mu,ndonors[[ch]][[ind]],donates[[ch]][[ind]],donatesl[[ch]][[ind]],donatesr[[ch]][[ind]])
      }
    }
  }
}
if (HPC==2)
{
  tmp<-foreach(ch_k=(1:(NUMA*nchrno))) %dopar% 
  {
    ch=as.integer((ch_k-0.5)/NUMA)+1
    k=(ch_k-1)%%NUMA+1
    ind<-as.integer((k+1)*0.5)
    donates_chr_ind=getdonates_ind(donates[[ch]][[ind]])
    donatesl_chr_ind=getdonates_ind(donatesl[[ch]][[ind]])
    donatesr_chr_ind=getdonates_ind(donatesr[[ch]][[ind]])
    ans=calc_E.n(ch,k,max.donors,NUMP,NUMA,G[ch],transitions[[ind]],flips[[ind]][[ch]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],
		 gobs[[ch]][[ind]],mutmat,maxmiss,kLL,L,PI,rho,Mu,ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
    rm(donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
    ans
  }
    for (ch in 1:nchrno)
    {
      E.n[[ch]]=list()
      for (k in 1:NUMA)
	E.n[[ch]][[k]]=tmp[[(ch-1)*NUMA+k]]
    }
}
cloglike=0;for(ch in 1:nchrno) for (k in 1:NUMA) cloglike=cloglike+E.n[[ch]][[k]]$loglike
# parameter updates are very fast given E.n[[]]
if (doMu)
{
  initi<-array(NaN,c(nchrno,NUMA,L,NUMP)) # a-posteriori first probs
  for (k in 1:NUMA)
  {
    for (ch in 1:nchrno) 
      for (k in 1:NUMA) {
	initi[ch,k,,]<-E.n[[ch]][[k]]$initi
      }
  }
  initg=array(NaN, c(nchrno,NUMA,L,kLL));
  for (k in 1:kLL) 
    #initg[,,,k]=apply(initi[,,,which(label==k)],-4,sum)
    initg[,,,k]=apply(array(initi[,,,which(label==k)],c(nchrno,NUMA,L,sum(label==k))),-4,sum)
  Mu[]<-0
  for (ja in 1:L)
  {
    for (jl in 1:kLL)
      for (ch in 1:nchrno)
	for (k in 1:NUMA)
	  Mu[jl,ja]<-Mu[jl,ja]+initg[ch,k,ja,jl]+sum(E.n[[ch]][[k]]$a[,jl,ja]) + E.n[[ch]][[k]]$r[jl,ja]
    if (all(Mu[,ja]==0)) Mu[,ja]=1/kLL
  }
  Mu<-t(t(Mu)/colSums(Mu))
}
if (doPI) # note that unlike the other parameters, each individual gets their own PI matrix 
{
  # note that E.n[[k]]$l[i] = sum(E.n[[k]]$na[,i])+sum(E.n[[k]]$a[,,i])
  # note that sum(E.n[[k]]$na[,i]) = sum(E.n[[k]]$n[,i])+sum(E.n[[k]]$r[,i])
  if (!exists("singlePI")) singlePI=F
  for (ind in 1:NUMI)
  {
    if (NUMA>1) {hap<-c(ind*2-1,ind*2)} else hap=1
    if (singlePI) hap=1:NUMA # use all to compute each; i.e. single PI, etc
    for (i in 1:L){
      denom=0
      for (k in hap) 
	for (ch in 1:nchrno)
	  denom=denom+sum(E.n[[ch]][[k]]$na[,i])+sum(E.n[[ch]][[k]]$a[i,,])
      for (j in 1:L) {
	numer=0
	for (k in hap) 
	  for (ch in 1:nchrno)
	    numer=numer+sum(E.n[[ch]][[k]]$a[i,,j])
	PI[[ind]][i,j]=sum(numer)/sum(denom) 
      }
      if (absorbrho) PI[[ind]][i,i]=0 # absorb into rho
    }
    tmp=which(is.na(PI[[ind]]),arr.ind=T) 
    for (i in tmp[,1]) for (j in tmp[,2]) PI[[ind]][i,j]=alpha[[ind]][j] # if never in an anc, just randomly choose another one w.p. alpha
    trans=PI[[ind]]-diag(rowSums(PI[[ind]])) # strictly speaking this should include an initial condition
    w=prcomp(t(trans),center=F)$rotation[,L]
    alpha[[ind]]=w/sum(w);alpha[[ind]][alpha[[ind]]<0]=0
    #lambda[[ind]]=-log(1-sum(PI[[ind]])+sum(diag(PI[[ind]])))/dr 
    tmp=1-t(t(PI[[ind]]/alpha[[ind]]));tmp[tmp==0]=NaN;tmp[tmp==1]=NaN;tmp[is.infinite(tmp)]=NaN;tmp[tmp<0]=NaN;diag(tmp)=NaN; # gives same off diagonals for L=2
    tmp=-log(tmp)/dr; # gives same off diagonals for L=2
    lambda[[ind]]=mean(tmp,na.rm=T)
    #lambda[[ind]]=mean(-log(1-PI[[ind]])/dr)
  }
}
if (dorho)
{
  numer=denom=rep(0,L)
  for (ch in 1:nchrno)
    for (k in 1:NUMA)
    {
      numer=numer+colSums(E.n[[ch]][[k]]$r) # E.n[[ch]][[k]][i,,i] all zero now
      denom=denom+colSums(E.n[[ch]][[k]]$r)+colSums(E.n[[ch]][[k]]$n)+apply(E.n[[ch]][[k]]$a,1,sum)
    }
  if (!commonrho) {rho<-numer/denom;rho[denom==0]=mean(rho[denom>0])}
  if (commonrho) rho[]<-sum(numer)/sum(denom) 
  #cat("!commonrho = ", numer/denom, "commonrho = ", sum(numer)/sum(denom), "\n")
}
if (dotheta)
{
  misses=0;matches=0
  for (k in 1:NUMA)
  {
    for (ch in 1:nchrno) 
    {
      misses=misses+E.n[[ch]][[k]]$e
      matches=matches+E.n[[ch]][[k]]$h
    }
  }
  if (!commontheta) {theta<-misses/(misses+matches);theta[(misses+matches)==0]=mean(theta[(misses+matches)>0])}
  if (commontheta) theta[]<-sum(misses)/sum(misses+matches) 
  mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
}
rm(E.n) 
if (dorho || doPI || doMu) 
{
  for (ind in 1:NUMI)
    transitions[[ind]]<-s_trans(L,kLL,PI[[ind]],Mu,rho,NL)
  initProb=initprobs()
}


