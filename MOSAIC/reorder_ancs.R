#alpha    chrnos   dr       g.loc    kLL      L        lambda   localanc
# Mu       nchrno   NL       NUMA     Q        rho      theta    tol 
if (length(neworder)!=L)
  stop(paste("fewer than "), L, " ancestries specified")
o.Mu<-o.Mu[,neworder]
for (ind in 1:NUMI)
{
  o.alpha[[ind]]<-o.alpha[[ind]][neworder]
  o.Q[[ind]]<-o.Q[[ind]][neworder,neworder]
}
o.rho<-o.rho[neworder]
o.theta<-o.theta[neworder]
Mu<-Mu[,neworder]
for (ind in 1:NUMI)
{
  alpha[[ind]]<-alpha[[ind]][neworder]
  Q[[ind]]<-Q[[ind]][neworder,neworder]
}
rho<-rho[neworder]
theta<-theta[neworder]
for (ch in 1:nchrno) localanc[[ch]]<-localanc[[ch]][neworder,,]
if (!exists("reordertruth")) reordertruth=T
if (target=="simulated" & reordertruth) 
{
  for (ch in 1:nchrno) g.true_anc[[ch]]<-g.true_anc[[ch]][neworder,,]
}
if (exists("gfbs"))
{
  tmp=expand.grid(1:kLL,neworder);tmp=c((tmp[,2]-1)*kLL+tmp[,1])
  for (ch in 1:nchrno) for (k in 1:NUMA) gfbs[[ch]][[k]]<-gfbs[[ch]][[k]][,tmp]
}
if (exists("colvec")) colvec=colvec[neworder]
