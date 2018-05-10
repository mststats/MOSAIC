# function to calculate the stationary state for the MOSAIC HMM transitions to use as initial probabilities at the first gridpoint in the forward algorithm
initprobs=function(fastcalc_init=T, t.NUMA=NUMA, t.L=L, t.NUMP=NUMP, t.kLL=kLL, t.PI=PI, t.Mu=Mu, t.rho=rho, t.alpha=alpha, t.label=label, t.NL=NL) 
{
  if (!fastcalc_init)
    allhaps_initProb=matrix(0,t.NUMA,t.L*t.NUMP)
  initProb=matrix(0,t.NUMA,t.kLL*t.L)
  for (k in 1:t.NUMA)
  {
    ind=as.integer((k+1)/2)
    if (!fastcalc_init)
    {
      trans_matrix<-matrix(0,t.L*t.NUMP,t.L*t.NUMP) # this must hold all haps initial probabilities
      for (ia in 1:t.L) 
	for (l in 1:t.L) 
	  for (ik in 1:t.NUMP) 
	  {
	    trans_matrix[(ia-1)*t.NUMP+ik,(l-1)*t.NUMP+ik]=explicit.trans(t.PI[[ind]],t.Mu,t.rho,t.NL,ia,T,l,t.label[ik]) #same hap
	    trans_matrix[(ia-1)*t.NUMP+(1:t.NUMP)[-ik],(l-1)*t.NUMP+ik]=explicit.trans(t.PI[[ind]],t.Mu,t.rho,t.NL,ia,F,l,t.label[ik]) #switch hap
	  }
      # from fiveMinuteStats
      lvec=Re(eigen(t(trans_matrix))$vectors) # Get the left eigenvectors of P
      allhaps_initProb[k,]=lvec[,1]/sum(lvec[,1]) # note columns used
      allhaps_initProb[k,]<-allhaps_initProb[k,]-min(allhaps_initProb[k,]) # translate to all positive
      if (sum(allhaps_initProb[k,])<0) allhaps_initProb[k,]=-allhaps_initProb[k,] # sign is arbitrary
      allhaps_initProb[k,]<-allhaps_initProb[k,]/sum(allhaps_initProb[k,])
      for (ia in 1:t.L)
	for (ik in 1:t.kLL)
	  initProb[k,(ia-1)*t.kLL+ik]<-sum(allhaps_initProb[k,t.label[1:t.NUMP]==ik]/t.NL[ik])
    }
    if (fastcalc_init) #  faster version almost the same; checked using fastcalc_init on and off
    {
      for (ia in 1:t.L)
	for (ik in 1:t.kLL)
	  initProb[k,(ia-1)*t.kLL+ik]<-t.Mu[ik,ia]/t.NL[ik]*t.alpha[[ind]][ia] 
    }
    #k_allhaps_initProb=rep(0,t.L*t.NUMP); for (l in 1:t.L) for (ik in 1:t.NUMP) k_allhaps_initProb[(l-1)*t.NUMP+ik]=initProb[k,(l-1)*t.kLL+t.label[ik]]/t.NL[t.label[ik]]
    # fastcalc_init usage checked using range(k_allhaps_initProb%*%trans_matrix-k_allhaps_initProb)
    initProb[k,]=initProb[k,]/sum(initProb[k,])
  }
  return(initProb)
}
