# function to quickly calculate the log-likelihood for the current MOSAIC HMM fit
get_loglike=function(t.NUMA, t.nchrno, t.G, t.L, t.kLL, t.max.donors, t.NUMP, t.donates, t.donatesl, t.transitions, t.maxmatchsize, t.umatch, t.flips,
		     t.mutmat, t.maxmiss, t.initProb) {
  kcloglike<-matrix(0,t.nchrno,t.NUMA)
  THIN=ifelse(t.max.donors==t.NUMP, F, T)
  for (ch in 1:t.nchrno) 
  {
    if (HPC==1)
    {
      donates_chr=getdonates(t.donates[[ch]],NUMI)
      donatesl_chr=getdonates(t.donatesl[[ch]],NUMI)
      tmp<-foreach(k=1:t.NUMA) %dopar%
      {
	ind=as.integer((k+1)/2)
	# fb calcs moved to here to avoid storing all fors, backs, etc
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.L);t.sumfors<-matrix(0,t.G[ch],t.L);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],label,ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	-sum(log(t.scalefactor))
      }
    }
    if (HPC==2)
    {
      tmp<-foreach(k=1:t.NUMA) %dopar%
      {
	ind=as.integer((k+1)/2)
	donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
	donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
	# fb calcs moved to here to avoid storing all fors, backs, etc
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.L);t.sumfors<-matrix(0,t.G[ch],t.L);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],label,ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	-sum(log(t.scalefactor))
      }
    }
    if (!HPC)
    {
      tmp<-foreach(k=1:t.NUMA) %dopar%
      {
	ind=as.integer((k+1)/2)
	# fb calcs moved to here to avoid storing all fors, backs, etc
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.L);t.sumfors<-matrix(0,t.G[ch],t.L);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.L,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],label,ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesl[[ch]][[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	-sum(log(t.scalefactor))
      }
    }
    kcloglike[ch,]=unlist(tmp)
  }
  cloglike<-sum(kcloglike)
  return(cloglike)
}
