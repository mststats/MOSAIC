# function to quickly calculate the log-likelihood for the current MOSAIC HMM fit
get_loglike=function(t.NUMA, t.nchrno, t.G, t.A, t.kLL, t.max.donors, t.NUMP, t.ndonors, t.donates, t.donatesl, t.transitions, t.maxmatchsize, t.umatch, t.flips,
		     t.mutmat, t.maxmiss, t.initProb, t.d.w, t.t.w, t.gobs, t.label, t.HPC) {
  kcloglike<-matrix(0,t.nchrno,t.NUMA)
  k=NULL # placeholder
  THIN=ifelse(t.max.donors==t.NUMP, FALSE, TRUE)
  t.NUMI=t.NUMA/2 # assuming diploid
  for (ch in 1:t.nchrno) 
  {
    if (t.HPC==1)
    {
      donates_chr=getdonates(t.donates[[ch]],t.NUMI)
      donatesl_chr=getdonates(t.donatesl[[ch]],t.NUMI)
      tmp<-foreach(k=1:t.NUMA) %dopar%
      {
	ind=as.integer((k+1)/2)
	# fb calcs moved to here to avoid storing all fors, backs, etc
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.A);t.sumfors<-matrix(0,t.G[ch],t.A);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	-sum(log(t.scalefactor))
      }
    }
    if (t.HPC==2)
    {
      tmp<-foreach(k=1:t.NUMA) %dopar%
      {
	ind=as.integer((k+1)/2)
	donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
	donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
	# fb calcs moved to here to avoid storing all fors, backs, etc
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.A);t.sumfors<-matrix(0,t.G[ch],t.A);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	-sum(log(t.scalefactor))
      }
    }
    if (!t.HPC)
    {
      tmp<-foreach(k=1:t.NUMA) %dopar%
      {
	ind=as.integer((k+1)/2)
	# fb calcs moved to here to avoid storing all fors, backs, etc
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.A);t.sumfors<-matrix(0,t.G[ch],t.A);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesl[[ch]][[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	-sum(log(t.scalefactor))
      }
    }
    kcloglike[ch,]=unlist(tmp)
  }
  cloglike<-sum(kcloglike)
  return(cloglike)
}
