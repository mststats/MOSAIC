# functions for calculating local ancestry without MOSAIC inferred re-phasing and without MOSAIC ancestry switches. 
# i.e. using only knowledge of which donor haplotypes are copied from 
get_ancunaware_localanc=function(t.NUMA,t.A,t.G,t.nchrno,t.noanc_gfbs,t.Mu,t.alpha){ 
  ancunaware_localanc<-list() # note that ancunaware refers to undoing the re-phasing steps made by MOSAIC
  for (ch in 1:t.nchrno) {
    ancunaware_localanc[[ch]]<-array(NaN,c(t.A,t.NUMA,t.G[ch]))
    for (k in 1:t.NUMA) {
      ind=as.integer((k+1)/2)
      for (a in 1:t.A) {
	ancunaware_localanc[[ch]][a,k,]=c(t.noanc_gfbs[[ch]][[k]]%*%t.Mu[,a])*t.alpha[[ind]][a]
      }
      ancunaware_localanc[[ch]][,k,]=t(t(ancunaware_localanc[[ch]][,k,])/colSums(ancunaware_localanc[[ch]][,k,]))
    }
  }
  return(ancunaware_localanc)
}
get_ancunaware_gfbs=function(t.NUMA,t.A,t.G,t.nchrno,t.noanc_gfbs,t.Mu,t.alpha){ 
  ancunaware_gfbs=list()
  for (ch in 1:t.nchrno){
    ancunaware_gfbs[[ch]]=list()
    for (k in 1:t.NUMA) {
      ind=as.integer((k+1)/2)
      ancunaware_gfbs[[ch]][[k]]<-matrix(NaN,t.G[ch],t.A*kLL)
      for (ll in 1:kLL) 
	for (a in 1:t.A) {
	  ancunaware_gfbs[[ch]][[k]][,(a-1)*kLL+ll]=t.Mu[ll,a]*t.noanc_gfbs[[ch]][[k]][,ll]*t.alpha[[ind]][a]
	}
      ancunaware_gfbs[[ch]][[k]]=ancunaware_gfbs[[ch]][[k]]/rowSums(ancunaware_gfbs[[ch]][[k]])
    }
  }
  return(ancunaware_gfbs)
}

