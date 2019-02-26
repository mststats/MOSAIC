# function to calculate the local ancestry estimates based on the gridded forward-backward probabilities gfbs
get_localanc=function(t.gfbs,t.G,t.L,t.kLL,t.NUMA,t.NUMI,tol=1e-8,t.g.true_anc=NULL) {
  ans=list()
  t.nchrno=length(t.gfbs)
  localanc<-list()
  for (ch in 1:t.nchrno)
  {
    localanc[[ch]]<-array(1/t.L,c(t.L,t.NUMA,t.G[ch]))
    for (k in 1:t.NUMA)
    {
      localanc[[ch]][,k,]<-0
      for (j in 1:t.L) 
	for (jk in 1:t.kLL) 
	  localanc[[ch]][j,k,]<-localanc[[ch]][j,k,] + t.gfbs[[ch]][[k]][,(j-1)*t.kLL+jk]
      #localanc[[ch]][,k,]<-round(localanc[[ch]][,k,],3) # Hapmix does this!
      localanc[[ch]][,k,][localanc[[ch]][,k,]<tol]<-tol;localanc[[ch]][,k,][localanc[[ch]][,k,]>(1-tol)]<-1-tol
      localanc[[ch]][,k,]<-t(t(localanc[[ch]][,k,])/apply(localanc[[ch]][,k,],2,sum))
    }
  }
  ans=list(localanc=localanc)
  # re-order truth to be same labels
  if (!is.null(t.g.true_anc))
  {
    r.ord<-function(ord) 
    {
      tmp.anc<-t.g.true_anc
      for (ch in 1:t.nchrno) tmp.anc[[ch]][,,]=t.g.true_anc[[ch]][ord,,]
      return(dip_fr2(tmp.anc,localanc))
    }
    all.ord<-permn(t.L)
    best.r<-0
    for (i in 1:length(all.ord))
      if (r.ord(all.ord[[i]])>best.r) 
      {
	best.r=r.ord(all.ord[[i]])
	best.ord<-all.ord[[i]]
      }
    for (ch in 1:t.nchrno) 
      t.g.true_anc[[ch]][,,]<-t.g.true_anc[[ch]][best.ord,,]
    ans$g.true_anc=t.g.true_anc
  }
  return(ans)
}

# function to find snp positions values of x based on gridded values of x
grid_to_pos=function(x,pos,glocs) { # arguments are thing-to-map, SNP positions, grid locations
  S=length(pos)
  g.map<-vapply(1:S, function(s) which.min((pos[s]-glocs)^2),0L) # create map from rates to grid
  if (length(dim(x))==2)
    ans=x[,g.map]
  if (length(dim(x))==3)
    ans=x[,,g.map]
}

# calculate proportions from local ancestry information
alpha_from_local=function(y) {ans=apply(sapply(y, function(x) apply(x,1,sum)),1,sum);return(ans/sum(ans))}
