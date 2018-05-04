# script to run iterations of the MCMC algorithm for re-phasing target genomes based on MOSAIC fit. 
# phase-hunter is far more efficient and finds approximate best phasing in practice
max.ll<-mcmc.ll<-mcmc.acc<-list()
for (ch in 1:nchrno) 
{
  max.ll[[ch]]<-mcmc.ll[[ch]]<-mcmc.acc[[ch]]<-list()
}
if (HPC!=2)
{
  for (ch in 1:nchrno)
  {
    if (HPC==1)
    {  
      donates_chr=getdonates(donates[[ch]],NUMI)
      donatesl_chr=getdonates(donatesl[[ch]],NUMI)
      donatesr_chr=getdonates(donatesr[[ch]],NUMI)
      tmp<-foreach(ind=1:NUMI) %dopar%
        phase_mcmc(ch,ind,M,max.donors,initProb,flips[[ind]][[ch]], verbose,ndonors[[ch]][[ind]], 
	  	   donates_chr[[ind]], donatesl_chr[[ind]], donatesr_chr[[ind]], transitions[[ind]], umatch[[ch]], maxmatchsize[ch], d.w[[ch]], t.w[[ch]], gobs[[ch]][[ind]], mutmat, maxmiss, PLOT=PLOT, mcmcprog)
    }
    if (!HPC)
    {
      tmp<-foreach(ind=1:NUMI) %dopar%
      phase_mcmc(ch,ind,M,max.donors,initProb,flips[[ind]][[ch]], verbose,ndonors[[ch]][[ind]], 
		 donates[[ch]][[ind]], donatesl[[ch]][[ind]], donatesr[[ch]][[ind]], transitions[[ind]], umatch[[ch]], maxmatchsize[ch], d.w[[ch]], t.w[[ch]], gobs[[ch]][[ind]], mutmat, maxmiss, PLOT=PLOT, mcmcprog)
    }
    for (ind in 1:NUMI)
    {
      if (length(tmp[[ind]]$ind.mcmc.acc)>0) 
      {
	flips[[ind]][[ch]]<-tmp[[ind]]$ind.max.flips
	mcmc.ll[[ch]][[ind]]<-tmp[[ind]]$ind.mcmc.ll
	mcmc.acc[[ch]][[ind]]<-tmp[[ind]]$ind.mcmc.acc
	max.ll[[ch]][[ind]]<-tmp[[ind]]$ind.max.ll
      }
    }
  }
}
if (HPC==2)
{  
  tmp<-foreach(ch_ind=1:(nchrno*NUMI)) %dopar%
  {
    ch=as.integer((ch_ind-0.5)/NUMI)+1
    ind=(ch_ind-1)%%NUMI+1
    donates_chr_ind=getdonates_ind(donates[[ch]][[ind]])
    donatesl_chr_ind=getdonates_ind(donatesl[[ch]][[ind]])
    donatesr_chr_ind=getdonates_ind(donatesr[[ch]][[ind]])
    ans=phase_mcmc(ch,ind,M,max.donors,initProb,flips[[ind]][[ch]], verbose,ndonors[[ch]][[ind]], 
	       donates_chr_ind, donatesl_chr_ind, donatesr_chr_ind, transitions[[ind]], umatch[[ch]], maxmatchsize[ch], d.w[[ch]], t.w[[ch]], gobs[[ch]][[ind]], mutmat, maxmiss, PLOT=PLOT, mcmcprog)
    ans
  }
  for (ch in 1:nchrno)
  {
    for (ind in 1:NUMI)
    {
      ch_ind=(ch-1)*NUMI+ind
      if (length(tmp[[ch_ind]]$ind.mcmc.acc)>0) 
      {
	flips[[ind]][[ch]]<-tmp[[ch_ind]]$ind.max.flips
	mcmc.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.mcmc.ll
	mcmc.acc[[ch]][[ind]]<-tmp[[ch_ind]]$ind.mcmc.acc
	max.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.max.ll
      }
    }
  }
}
rm(tmp)
cloglike<-sum(unlist((max.ll)))
if (LOG) 
{
  runtime<-as.numeric(Sys.time());diff.time<-runtime-old.runtime;old.runtime<-runtime;
  writelog(EMlogfile,"phasemcmc",diff.time,len)
}
