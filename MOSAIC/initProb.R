#source("transitions.R")
fastcalc_init=T
if (!fastcalc_init)
  allhaps_initProb=matrix(0,NUMA,L*NUMP)
initProb=matrix(0,NUMA,kLL*L)
for (k in 1:NUMA)
{
  ind=as.integer((k+1)/2)
  if (!fastcalc_init)
  {
    trans_matrix<-matrix(0,L*NUMP,L*NUMP) # this must hold all haps initial probabilities
    for (ia in 1:L) 
      for (l in 1:L) 
	for (ik in 1:NUMP) 
	{
	  trans_matrix[(ia-1)*NUMP+ik,(l-1)*NUMP+ik]=explicit.trans(PI[[ind]],Mu,rho,NL,ia,T,l,label[ik]) #same hap
	  trans_matrix[(ia-1)*NUMP+(1:NUMP)[-ik],(l-1)*NUMP+ik]=explicit.trans(PI[[ind]],Mu,rho,NL,ia,F,l,label[ik]) #switch hap
	}
    # from fiveMinuteStats
    lvec=Re(eigen(t(trans_matrix))$vectors) # Get the left eigenvectors of P
    allhaps_initProb[k,]=lvec[,1]/sum(lvec[,1]) # note columns used
    allhaps_initProb[k,]<-allhaps_initProb[k,]-min(allhaps_initProb[k,]) # translate to all positive
    if (sum(allhaps_initProb[k,])<0) allhaps_initProb[k,]=-allhaps_initProb[k,] # sign is arbitrary
    allhaps_initProb[k,]<-allhaps_initProb[k,]/sum(allhaps_initProb[k,])
    for (ia in 1:L)
      for (ik in 1:kLL)
	initProb[k,(ia-1)*kLL+ik]<-sum(allhaps_initProb[k,label[1:NUMP]==ik]/NL[ik])
  }
  if (fastcalc_init) #  faster version almost the same; checked using fastcalc_init on and off
  {
    for (ia in 1:L)
      for (ik in 1:kLL)
	initProb[k,(ia-1)*kLL+ik]<-Mu[ik,ia]/NL[ik]*alpha[[ind]][ia] 
  }
  #k_allhaps_initProb=rep(0,L*NUMP); for (l in 1:L) for (ik in 1:NUMP) k_allhaps_initProb[(l-1)*NUMP+ik]=initProb[k,(l-1)*kLL+label[ik]]/NL[label[ik]]
  # checked using range(k_allhaps_initProb%*%trans_matrix-k_allhaps_initProb)
  initProb[k,]=initProb[k,]/sum(initProb[k,])
}
