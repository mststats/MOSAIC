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
# re-order truth to be same labels
if (target=="simulated")
{
  require(combinat)
  source("calc_r2.R")
  r.ord<-function(ord) 
  {
    tmp.anc<-g.true_anc
    for (ch in 1:nchrno) tmp.anc[[ch]][,,]=g.true_anc[[ch]][ord,,]
    return(dip_fr2(tmp.anc,localanc))
  }
  all.ord<-permn(L)
  best.r<-0
  for (i in 1:length(all.ord))
    if (r.ord(all.ord[[i]])>best.r) 
    {
      best.r=r.ord(all.ord[[i]])
      best.ord<-all.ord[[i]]
    }
  for (ch in 1:nchrno) 
    g.true_anc[[ch]][,,]<-g.true_anc[[ch]][best.ord,,]
  for (k in 1:NUMI) 
    sim.alpha[[k]]<-sim.alpha[[k]][best.ord]
}
