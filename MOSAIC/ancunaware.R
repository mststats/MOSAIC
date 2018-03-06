# script for calculating local ancestry without MOSAIC inferred re-phasing and without MOSAIC ancestry switches. 
# i.e. using only knowledge of which donor haplotypes are copied from 
#for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # this is only required if not reloading but working with current session
get_ancunaware_localanc=function(){ 
  ancunaware_localanc<-list() # note that ancunaware refers to undoing the re-phasing steps made by MOSAIC
  for (ch in 1:nchrno) {
    ancunaware_localanc[[ch]]<-array(NaN,c(L,NUMA,G[ch]))
    for (k in 1:NUMA) {
      ind=as.integer((k+1)/2)
      for (a in 1:L) {
	ancunaware_localanc[[ch]][a,k,]=c(noanc_gfbs[[ch]][[k]]%*%Mu[,a])*alpha[[ind]][a]
      }
      ancunaware_localanc[[ch]][,k,]=t(t(ancunaware_localanc[[ch]][,k,])/colSums(ancunaware_localanc[[ch]][,k,]))
    }
  }
  return(ancunaware_localanc)
}
get_ancunaware_gfbs=function(){
  ancunaware_gfbs=list()
  for (ch in 1:nchrno){
    ancunaware_gfbs[[ch]]=list()
    for (k in 1:NUMA) {
      ind=as.integer((k+1)/2)
      ancunaware_gfbs[[ch]][[k]]<-matrix(NaN,G[ch],L*kLL)
      for (ll in 1:kLL) 
	for (a in 1:L) {
	  ancunaware_gfbs[[ch]][[k]][,(a-1)*kLL+ll]=Mu[ll,a]*noanc_gfbs[[ch]][[k]][,ll]*alpha[[ind]][a]
	}
      ancunaware_gfbs[[ch]][[k]]=ancunaware_gfbs[[ch]][[k]]/rowSums(ancunaware_gfbs[[ch]][[k]])
    }
  }
  return(ancunaware_gfbs)
}

