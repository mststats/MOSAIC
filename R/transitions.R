# function to calculate elements of the MOSAIC HMM transition matrix
explicit.trans<-function(t.PI,t.Mu,t.rho,t.NL,t.i,t.kequalsn,t.l,t.ll)
{
  ans<-NaN
  if (t.i!=t.l) # p(anc)*p(hap)
    return(t.PI[t.i,t.l]*t.Mu[t.ll,t.l]/t.NL[t.ll])
  if (t.i==t.l)
  {
    if (!t.kequalsn) # (p(no anc)*p(hap)+p(anc to same))*p(hap)
      ans=((1-sum(t.PI[t.i,]))*t.rho[t.i]+t.PI[t.i,t.i])*t.Mu[t.ll,t.i]/t.NL[t.ll]
    if (t.kequalsn) # (p(no anc)*p(hap)+p(anc to same))*p(hap) + p(no anc)*p(no hap)
      ans=((1-sum(t.PI[t.i,]))*t.rho[t.i]+t.PI[t.i,t.i])*t.Mu[t.ll,t.i]/t.NL[t.ll] + (1-sum(t.PI[t.i,]))*(1-t.rho[t.i])
  }
  #if(ans<1e-3) ans<-1e-3
  return(ans)
}
