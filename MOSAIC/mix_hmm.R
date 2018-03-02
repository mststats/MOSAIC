source("transitions.R")
#### functions to return transition probabilities ####
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
require(Rcpp)
# note that these functions cannot be re-loaded in a new R session but must be re-compiled again
#sourceCpp("fors.cpp")
#sourceCpp("backs.cpp")
#sourceCpp("forbacks.cpp")

get_gfbs<-function()
{
  ans<-list()
  THIN=ifelse(max.donors==NUMP, F, T)
  for (ch in 1:nchrno)
  {
    if (HPC==1)
    {
      donates_chr=getdonates(donates[[ch]],NUMI)
      donatesl_chr=getdonates(donatesl[[ch]],NUMI)
      donatesr_chr=getdonates(donatesr[[ch]],NUMI)
      ans[[ch]]<-list()
      tmp<-foreach(k=1:NUMA) %dopar% 
      {
        ind=as.integer((k+1)/2)
        t.fors<-rep(0,G[ch]*max.donors*L);t.sumfors<-matrix(0,G[ch],L);t.scalefactor<-rep(0,G[ch]);
        cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,initProb[k,],label,
		   ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
        t.backs<-rep(0,G[ch]*max.donors*L);t.scalefactorb<-rep(0,G[ch]);
        cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,label,
	  	    ndonors[[ch]][[ind]],donates_chr[[ind]],donatesr_chr[[ind]],flips[[ind]][[ch]],t.backs,t.scalefactorb)
        cppgforback(max.donors,THIN,kLL,NUMP,label,L,G[ch],ndonors[[ch]][[ind]],donates_chr[[ind]],t.fors,t.backs)
      }
    }
    if (HPC==2)
    {
      ans[[ch]]<-list()
      tmp<-foreach(k=1:NUMA) %dopar% 
      {
        ind=as.integer((k+1)/2)
        donates_chr_ind=getdonates_ind(donates[[ch]][[ind]])
        donatesl_chr_ind=getdonates_ind(donatesl[[ch]][[ind]])
        donatesr_chr_ind=getdonates_ind(donatesr[[ch]][[ind]])
        t.fors<-rep(0,G[ch]*max.donors*L);t.sumfors<-matrix(0,G[ch],L);t.scalefactor<-rep(0,G[ch]);
        cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,initProb[k,],label,
		   ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
        t.backs<-rep(0,G[ch]*max.donors*L);t.scalefactorb<-rep(0,G[ch]);
        cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,label,
	  	    ndonors[[ch]][[ind]],donates_chr_ind,donatesr_chr_ind,flips[[ind]][[ch]],t.backs,t.scalefactorb)
        ans=cppgforback(max.donors,THIN,kLL,NUMP,label,L,G[ch],ndonors[[ch]][[ind]],donates_chr_ind,t.fors,t.backs)
	rm(donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
	ans
      }
    }
    if (!HPC)
    {
      ans[[ch]]<-list()
      tmp<-foreach(k=1:NUMA) %dopar% 
      {
        ind=as.integer((k+1)/2)
        t.fors<-rep(0,G[ch]*max.donors*L);t.sumfors<-matrix(0,G[ch],L);t.scalefactor<-rep(0,G[ch]);
        cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,initProb[k,],label,
		   ndonors[[ch]][[ind]],donates[[ch]][[ind]],donatesl[[ch]][[ind]],flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
        t.backs<-rep(0,G[ch]*max.donors*L);t.scalefactorb<-rep(0,G[ch]);
        cppbackward(k,NUMA,max.donors,THIN,NUMP,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,label,
	  	    ndonors[[ch]][[ind]],donates[[ch]][[ind]],donatesr[[ch]][[ind]],flips[[ind]][[ch]],t.backs,t.scalefactorb)
        cppgforback(max.donors,THIN,kLL,NUMP,label,L,G[ch],ndonors[[ch]][[ind]],donates[[ch]][[ind]],t.fors,t.backs)
      }
    }
    ans[[ch]]<-tmp
  }
  return(ans)
}

