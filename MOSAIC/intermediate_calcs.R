# calculations related to EM updates used in EM_updates.R
require(compiler)
# observed number of switches from (ia,il) to (ja,jl) 
r.calc_E.n<-function(ch,k,maxdonors,NUMP,NUMA,t.G,t.transitions,t.flips,t.umatch,t.maxmatchsize,t.dw,t.tw,r_gobs,mutmat,maxmiss,kLL,L,PI,rho,Mu,t.ndonors,t.donates,t.donatesl,t.donatesr) 
{
  THIN=ifelse(maxdonors==NUMP,F,T)
  a<-array(0,c(L,kLL,L)) # need to track which groups are switched to for Mu updates
  r<-matrix(0,kLL,L) # need to track which groups are switched to for Mu updates
  na<-matrix(0,kLL,L) # no anc switch, may or may not switch hap
  n<-matrix(0,kLL,L) # nothing happens; no anc switch and no hap switch
  # fb calcs moved to here to avoid storing all fors, backs, etc
  t.fors<-rep(0,t.G*maxdonors*L);t.sumfors<-matrix(0,t.G,L);t.scalefactor<-rep(0,t.G);
  cppforward(k,NUMA,maxdonors,THIN,NUMP,kLL,L,0,t.G,t.G,t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,r_gobs,mutmat,maxmiss,initProb[k,],label,
	     t.ndonors,t.donates,t.donatesl,t.flips,t.fors,t.sumfors,t.scalefactor)
  t.backs<-rep(0,t.G*maxdonors*L);t.scalefactorb<-rep(0,t.G);
  cppbackward(k,NUMA,maxdonors,THIN,NUMP,L,0,t.G,t.G,t.transitions,t.umatch,t.maxmatchsize,t.dw,t.tw,r_gobs,mutmat,maxmiss,label,
	      t.ndonors,t.donates,t.donatesr,t.flips,t.backs,t.scalefactorb)
  probs<-cppprobs(k,NUMA,maxdonors,THIN,L,kLL,NN,NUMP,t.G,label,t.fors,t.sumfors,t.backs,t.transitions,t.flips,mutmat,maxmiss,t.umatch,t.maxmatchsize,t.dw,t.tw,r_gobs,t.ndonors,t.donates,t.donatesl)
  for (type in 1:length(probs)) probs[[type]][probs[[type]]<0]=0
  #probs$switches are switches that go from (il,ia) to (jl,ja) but doesn't include same haps i.e. must be a switch to a different hap
  #probs$self are both switches and non-switches that go to same ancestry / hap pair
  #sum(probs$switches)+sum(probs$self)=G-1
  ind=as.integer((k+1)/2)
  for (ia in 1:L) {
    for (ja in {1:L}[-ia]) # 1. must have changed ancestry => psa=1
    {
      for (jl in 1:kLL)
	a[ia,jl,ja]<-a[ia,jl,ja]+probs$switches[ia,jl,ja]
    }
    # now do ja=ia and jk!=ik; use current paras to weight which switch type may have occurred:
    # (1) switched to same anc (2) no anc switch but switched hap
    for (jl in (1:kLL)) 
    {
      psa<-PI[[ind]][ia,ia]
      npsa<-1-sum(PI[[ind]][ia,]) 
      psr<-rho[ia]
      npsr<-1-rho[ia]
      spsa<-psa/(psa+psr*npsa) # relative contribution of ancestry switches
      sweights<-c(spsa, 1-spsa)*probs$switches[ia,jl,ia] #, 0) # psr=1 only when different hap but probs$switches already doesn't contain that move
      a[ia,jl,ia]<-a[ia,jl,ia]+sweights[1] # switches ancestry 
      na[jl,ia]<-na[jl,ia]+sweights[2] # doesn't switch ancestry but switches hap 
      r[jl,ia]<-r[jl,ia]+sweights[2] # doesn't switch ancestry but switches hap
      # now do jk=ik 
      spsa<-(     psa*Mu[jl,ia]/NL[jl])/((psa+psr*npsa)*Mu[jl,ia]/NL[jl]+npsr*npsa) # relative contribution of ancestry switches
      spsr<-(npsa*psr*Mu[jl,ia]/NL[jl])/((psa+psr*npsa)*Mu[jl,ia]/NL[jl]+npsr*npsa) # relative contribution of recombination switches
      for (ik in which(label==jl)) # 1 may have switched to same anc 2 may have switched to same hap w/ no anc switch 3 may not have switched
      {
	sweights<-c(spsa, spsr, 1-spsa-spsr)*probs$self[ik,ia] #, no switch prob is {1-spsa}*{1-spsr}) = 1-spsa-spsr
	a[ia,jl,ia]<-a[ia,jl,ia]+sweights[1] # switches ancestry
	r[jl,ia]<-r[jl,ia]+sweights[2] # doesn't switch ancestry and switches hap  
	na[jl,ia]<-na[jl,ia]+sweights[2]+sweights[3] # doesn't switch ancestry and switches hap  
	n[jl,ia]<-n[jl,ia]+sweights[3] # doesn't switch ancestry and doesn't switch hap
      }
    }
  }
  initi=NaN
  if (doMu) 
    initi<-t(matrix(cppforback(maxdonors,THIN,NUMP,L,1,t.ndonors,t.donates,t.fors,t.backs),NUMP)) # a-posteriori first gridpoint probs
  l<-apply(probs$switches,3,sum)+apply(probs$self,2,sum) # all switches into anc=1:L times the distance from the last locus
  loglike=-sum(log(t.scalefactor))
  a[a<0]=0;r[r<0]=0;n[n<0]=0;na[na<0]=0
  a[is.nan(a)]=0;r[is.nan(r)]=0;n[is.nan(n)]=0;na[is.nan(na)]=0
  list(a=a,r=r,na=na,n=n,l=l,e=probs$errors,h=probs$hits,initi=initi,loglike=loglike) 
}
#  calc_E.n<-cmpfun(r.calc_E.n,list(optimize=optlevel))
calc_E.n<-r.calc_E.n
#sum(E.n[[k]]$a)+sum(E.n[[k]]$r)+sum(E.n[[k]]$n)=G-1

create_PI<-function(alpha,lambda,L,dr,NUMI)
{
  PI<-list()
  for (ind in 1:NUMI)
  {
    PI[[ind]]<-matrix(0,L,L)
    for (i in 1:L)
      for (j in 1:L)
	PI[[ind]][i,j]=(1-exp(-dr*lambda[[ind]]))*alpha[[ind]][j] # rate of actual ancestry switching
  }
  PI
}

