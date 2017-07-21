kcloglike<-matrix(0,nchrno,NUMA)
THIN=ifelse(max.donors==NUMP, F, T)
for (ch in 1:nchrno) 
{
  if (HPC==1)
  {
    donates_chr=getdonates(donates[[ch]],NUMI)
    donatesl_chr=getdonates(donatesl[[ch]],NUMI)
    tmp<-foreach(k=1:NUMA) %dopar%
    {
      ind=as.integer((k+1)/2)
      # fb calcs moved to here to avoid storing all fors, backs, etc
      t.fors<-rep(0,G[ch]*max.donors*L);t.sumfors<-matrix(0,G[ch],L);t.scalefactor<-rep(0,G[ch]);
      cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,initProb[k,],label,
	         ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
      -sum(log(t.scalefactor))
    }
  }
  if (HPC==2)
  {
    tmp<-foreach(k=1:NUMA) %dopar%
    {
      ind=as.integer((k+1)/2)
      donates_chr_ind=getdonates_ind(donates[[ch]][[ind]])
      donatesl_chr_ind=getdonates_ind(donatesl[[ch]][[ind]])
      # fb calcs moved to here to avoid storing all fors, backs, etc
      t.fors<-rep(0,G[ch]*max.donors*L);t.sumfors<-matrix(0,G[ch],L);t.scalefactor<-rep(0,G[ch]);
      cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,initProb[k,],label,
	         ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
      -sum(log(t.scalefactor))
    }
  }
  if (!HPC)
  {
    tmp<-foreach(k=1:NUMA) %dopar%
    {
      ind=as.integer((k+1)/2)
      # fb calcs moved to here to avoid storing all fors, backs, etc
      t.fors<-rep(0,G[ch]*max.donors*L);t.sumfors<-matrix(0,G[ch],L);t.scalefactor<-rep(0,G[ch]);
      cppforward(k,NUMA,max.donors,THIN,NUMP,kLL,L,0,G[ch],G[ch],transitions[[ind]],umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],mutmat,maxmiss,initProb[k,],label,
	         ndonors[[ch]][[ind]],donates[[ch]][[ind]],donatesl[[ch]][[ind]],flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
      -sum(log(t.scalefactor))
    }
  }
  kcloglike[ch,]=unlist(tmp)
}
cloglike<-sum(kcloglike)
