# calculate the local ancestry estimates based on the gridded forward-backward probabilities gfbs
if (!exists("G")) G=sapply(gfbs, function(x) nrow(x[[1]]))
localanc<-list()
for (ch in 1:nchrno)
{
  localanc[[ch]]<-array(1/L,c(L,NUMA,G[ch]))
  for (k in 1:NUMA)
  {
    localanc[[ch]][,k,]<-0
    for (j in 1:L) 
      for (jk in 1:kLL) 
	localanc[[ch]][j,k,]<-localanc[[ch]][j,k,] + gfbs[[ch]][[k]][,(j-1)*kLL+jk]
    #localanc[[ch]][,k,]<-round(localanc[[ch]][,k,],3) # Hapmix does this!
    localanc[[ch]][,k,][localanc[[ch]][,k,]<tol]<-tol;localanc[[ch]][,k,][localanc[[ch]][,k,]>(1-tol)]<-1-tol
    localanc[[ch]][,k,]<-t(t(localanc[[ch]][,k,])/apply(localanc[[ch]][,k,],2,sum))
  }
}
