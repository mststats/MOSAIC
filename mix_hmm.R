# transition probabilities function and function to return forward-backward probabilities using Cpp code in forbacks.cpp
#### function to return transition probabilities ####
s_trans<-function(t.L,t.kLL,t.PI,t.Mu,t.rho,t.NL) 
{
  knvec<-c(F,T)
  ans<-array(NaN, c(t.L,2,t.L,t.kLL)) 
  for (t.i in 1:t.L) 
    for (t.kn in 1:2)
      for (t.l in 1:t.L) 
	for (t.ll in 1:t.kLL) 
	  ans[t.i,t.kn,t.l,t.ll]<-explicit.trans(t.PI,t.Mu,t.rho,t.NL,t.i,knvec[t.kn],t.l,t.ll) 
  ans
}

get_gfbs<-function(t.NUMP, t.nchrno, t.max.donors, t.donates, t.donatesl, t.donatesr, t.NUMA, t.L, t.G, t.kLL, t.transitions, t.umatch, t.maxmatchsize, t.d.w, t.t.w, t.gobs, 
		   t.mutmat, t.maxmiss, t.initProb, t.label, t.ndonors, t.flips, t.HPC)
{
  ans<-list()
  THIN=ifelse(t.max.donors==t.NUMP, F, T)
  for (ch in 1:t.nchrno)
  {
    if (t.HPC==1)
    {
      donates_chr=getdonates(t.donates[[ch]],NUMI)
      donatesl_chr=getdonates(t.donatesl[[ch]],NUMI)
      donatesr_chr=getdonates(t.donatesr[[ch]],NUMI)
      ans[[ch]]<-list()
      tmp<-foreach(k=1:t.NUMA) %dopar% 
      {
	ind=as.integer((k+1)/2)
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.L);t.sumfors<-matrix(0,t.G[ch],t.L);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	t.backs<-rep(0,t.G[ch]*t.max.donors*t.L);t.scalefactorb<-rep(0,t.G[ch]);
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		    t.mutmat,t.maxmiss,t.label,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesr_chr[[ind]],t.flips[[ind]][[ch]],t.backs,t.scalefactorb)
	cppgforback(t.max.donors,THIN,t.kLL,t.NUMP,t.label,t.L,t.G[ch],t.ndonors[[ch]][[ind]],donates_chr[[ind]],t.fors,t.backs)
      }
    }
    if (t.HPC==2)
    {
      ans[[ch]]<-list()
      tmp<-foreach(k=1:t.NUMA) %dopar% 
      {
	ind=as.integer((k+1)/2)
	donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
	donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
	donatesr_chr_ind=getdonates_ind(t.donatesr[[ch]][[ind]])
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.L);t.sumfors<-matrix(0,t.G[ch],t.L);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	t.backs<-rep(0,t.G[ch]*t.max.donors*t.L);t.scalefactorb<-rep(0,t.G[ch]);
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		    t.mutmat,t.maxmiss,t.label,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesr_chr_ind,t.flips[[ind]][[ch]],t.backs,t.scalefactorb)
	ans=cppgforback(t.max.donors,THIN,t.kLL,t.NUMP,t.label,t.L,t.G[ch],t.ndonors[[ch]][[ind]],donates_chr_ind,t.fors,t.backs)
	rm(donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
	ans
      }
    }
    if (!t.HPC)
    {
      ans[[ch]]<-list()
      tmp<-foreach(k=1:t.NUMA) %dopar% 
      {
	ind=as.integer((k+1)/2)
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.L);t.sumfors<-matrix(0,t.G[ch],t.L);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesl[[ch]][[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	t.backs<-rep(0,t.G[ch]*t.max.donors*t.L);t.scalefactorb<-rep(0,t.G[ch]);
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.w[[ch]],t.gobs[[ch]][[ind]],
		    t.mutmat,t.maxmiss,t.label,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesr[[ch]][[ind]],t.flips[[ind]][[ch]],t.backs,t.scalefactorb)
	cppgforback(t.max.donors,THIN,t.kLL,t.NUMP,t.label,t.L,t.G[ch],t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.fors,t.backs)
      }
    }
    ans[[ch]]<-tmp
  }
  return(ans)
}

