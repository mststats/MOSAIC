window_chunks<-function(nswitches, ww=0.5, min.swiches=0e-4, verbose=F) # windows of 0.5M sum(G)*dr/w*100~=7000 for w=0.5
{
  nchrno=length(nswitches)
  #for (ch in 1:nchrno) nswitches[[ch]]<-nswitches[[ch]][-ignorepanels,,]
  # convert into window lengths in terms of #gridpoints
  #nw=sum(G)*dr/w*100; w=sum(G)/nw 
  w=ww/dr/100
  if (verbose) cat("Initialising copying matrix Mu;", w, "gridpoints per", ww, "cM width window\n") 
  nw=as.integer(G/w);offsetw=cumsum(nw)-nw
  wmat<-list()
  for (h in 1:NUMA)
  {
    wmat[[h]]=matrix(0,kLL,sum(nw)) # collate switching across chromosomes and haplotypes
    for (ch in 1:nchrno)
    {
      nswitches[[ch]][nswitches[[ch]]<min.swiches]=NaN
      for (i in 1:nw[ch])
      {
        wind=((i-1)*w+1):(i*w)
	wmat[[h]][,i+offsetw[ch]]=rowSums(nswitches[[ch]][,wind,h],na.rm=T) # sum over window 
      }
      if ((nw[ch]*w)<G[ch]) # if any left over just add in to last one
      {
        wind<-(nw[ch]*w+1):G[ch]
       	wmat[[h]][,i+offsetw[ch]]=wmat[[h]][,i+offsetw[[ch]]]+rowSums(as.matrix(nswitches[[ch]][,wind,h]),na.rm=T) # sum over window 
      }
    }
  }
  return(list(wmat=wmat,w=w))
}
cluster_windows<-function(windows,PLOT=F,t.L=L,verbose=F)
{
  require(gtools)
  source("EMmult.R")
  nw=ncol(windows$wmat[[1]])
  res<-EMmult(windows$wmat, t.L, verbose=verbose) # fit the mixture model
  Mu<-res$Mu*NL[1:kLL];Mu<-t(t(Mu)/colSums(Mu)) # note the use of NL to scale by panel size
  alpha<-res$alpha
  gpw=windows$w # gridpoints per window
  lambda=list()
  #probs=list();probs[[1]]=list() # like one long chromosome here
  #for (h in 1:NUMA)
  #{
  #  probs[[1]][[h]]=matrix(NaN,nw,kLL*t.L)
  #  # ignoring the groups here, using P(a|s); some redundancy required for compatibility with coancs usage
  #  for (k in 1:kLL) for(l in 1:t.L) probs[[1]][[h]][,(l-1)*kLL+k]=res$p[,h,l]/kLL 
  #}
  #source("coancestry.R");acoancs<-create_coancs(probs,dr*gpw,"DIP")
  #if (!exists("singlelambda")) singlelambda=T
  #if (!singlelambda)
  #  for (ind in 1:NUMI)
  #    lambda[[ind]]=plot_coanccurves(acoancs,dr*gpw,PLOT=F,verbose=F,k=ind)$params[,,3]
  #if (singlelambda)
  #{
  #  tmplambda=plot_coanccurves(acoancs,dr*gpw,PLOT=F,verbose=F,k=NULL)$params[,,3]
  #  for (ind in 1:NUMI)
  #    lambda[[ind]]=tmplambda
  #}
  Q<-list()
  for (ind in 1:NUMI) 
  {
    hap=c(ind,ind);if (NUMA>1) hap=c(ind*2-1,ind*2) # haploid or diploid
    Q[[ind]]=matrix(0,t.L,t.L)
    for (i in 1:t.L)
    {
      for (j in 1:t.L)
      {
	# average p(i->j)/p(i) across windows and then average this over both haplotypes
        Q[[ind]][i,j]=0.5*(sum(res$p[-nw,hap[1],i]*res$p[-1,hap[1],j])/sum(res$p[-nw,hap[1],i])+
                           sum(res$p[-nw,hap[2],i]*res$p[-1,hap[2],j])/sum(res$p[-nw,hap[2],i]))
      }
      if (absorbrho) Q[[ind]][i,i]=0
    }
    # diagonal includes non-anc-switches. Need to remove these next; expect switches to same anc roughly the same frequency as switches from another anc
    #for (i in 1:t.L)
    #  Q[[ind]][i,i]=Q[[ind]][i,i]=mean(Q[[ind]][-i,i]) 
    #tmp=-log(1-Q[[ind]])/(dr*gpw);Q[[ind]]=1-exp(-dr*tmp); # rescale Q to grid widths from window widths
    Q[[ind]]=1-(1-Q[[ind]])^(1/gpw) # rescale Q to grid widths from window widths
    Q[[ind]]=Q[[ind]]/5 # ad-hoc reduction as this is typically overestimated
    tmp=which(is.na(Q[[ind]]),arr.ind=T) 
    for (i in tmp[,1]) for (j in tmp[,2]) Q[[ind]][i,j]=alpha[[ind]][j] # if never in an anc, just randomly choose another one w.p. alpha
    trans=Q[[ind]]-diag(rowSums(Q[[ind]])) # strictly speaking this should include an initial condition
    w=prcomp(t(trans),center=F)$rotation[,L]
    alpha[[ind]]=w/sum(w);alpha[[ind]][alpha[[ind]]<0]=0 # very close to the above alpha but now consistent with Q estimates
    lambda[[ind]]=-log(1-Q[[ind]])/dr
  }
  #if (!commonrho) rho=rho-rowMeans(sapply(Q,function(x) (rowSums(x)+colSums(x)))) # all anc switches in and out
  #if (commonrho) rho=rho-mean(sapply(Q,sum))
  rownames(Mu)<-panels[1:kLL]#[-ignorepanels]
  if (PLOT)
  {
    if (!exists("PNG")) PNG=T
    cexa=ifelse(PNG, 3, 1.5)
    if (PNG) png(file=paste0("PLOTS/",target,"_",t.L,"-way","_initial_Mu.png"),width=920,height=1490)
    par(mfrow=c(1,1))
    source("plot_funcs.R")
    plot_Mu(Mu,alpha,NL,MODE="copy",cexa=cexa)
    if (PNG) dev.off()
    if (!PNG) dev.new()
    if (PNG) png(file=paste0("PLOTS/",target,"_",t.L,"-way","_initial_curves.png"),width=1490,height=1490)
    ans=plot_coanccurves(acoancs,dr*windows$w,PLOT=T,verbose=verbose,cexa=cexa)
    if (PNG) dev.off()
  }
  return(list(Mu=Mu,alpha=alpha,lambda=lambda,Q=Q,ll=res$ll)) 
}
