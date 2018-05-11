# script to run the fast phase-hunting algorithm to re-phase target genomes using MOSAIC fit
# not written as a function to avoid copying of large flips object
orig.ll<-max.ll<-c.ll<-list()
nflips=0
niters=0
if (verbose) cat("re-phasing... ")
if (HPC!=2)
{
  for (ch in 1:nchrno) 
  {
    if (HPC==1)
    {
      donates_chr=getdonates(donates[[ch]],NUMI)
      donatesl_chr=getdonates(donatesl[[ch]],NUMI)
      donatesr_chr=getdonates(donatesr[[ch]],NUMI)
      orig.ll[[ch]]<-max.ll[[ch]]<-c.ll[[ch]]<-list()
      tmp<-foreach(ind=1:NUMI) %dopar%
	phase_hunt(eps.lower,ch,ind,flips[[ind]][[ch]], FALSE, ndonors[[ch]][[ind]], donates_chr[[ind]], donatesl_chr[[ind]], donatesr_chr[[ind]], 
		   transitions[[ind]], umatch[[ch]], maxmatchsize[ch], d.w[[ch]], t.w[[ch]], gobs[[ch]][[ind]], mutmat, maxmiss, phase.error.locs, 
		   initProb, PLOT=PLOT, minbg=min.bg, maxbg=max.bg)
    }
    if (!HPC)
    {
      orig.ll[[ch]]<-max.ll[[ch]]<-c.ll[[ch]]<-list()
      tmp<-foreach(ind=1:NUMI) %dopar%
	phase_hunt(eps.lower,ch,ind,flips[[ind]][[ch]], FALSE, ndonors[[ch]][[ind]], donates[[ch]][[ind]], donatesl[[ch]][[ind]], donatesr[[ch]][[ind]], 
		   transitions[[ind]], umatch[[ch]], maxmatchsize[ch], d.w[[ch]], t.w[[ch]], gobs[[ch]][[ind]], mutmat, maxmiss, phase.error.locs, 
		   initProb, PLOT=PLOT, minbg=min.bg, maxbg=max.bg)
    }
    for (ind in 1:(NUMI)) 
    {
      flips[[ind]][[ch]]<-tmp[[ind]]$ind.max.flips
      max.ll[[ch]][[ind]]<-c.ll[[ch]][[ind]]<-tmp[[ind]]$ind.max.ll
      orig.ll[[ch]][[ind]]<-tmp[[ind]]$ind.orig.ll
      nflips=nflips+tmp[[ind]]$nflips
      niters=niters+tmp[[ind]]$niters
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
    ans=phase_hunt(eps.lower,ch,ind,flips[[ind]][[ch]], F, ndonors[[ch]][[ind]], donates_chr_ind, donatesl_chr_ind, donatesr_chr_ind, 
		   transitions[[ind]], umatch[[ch]], maxmatchsize[ch], d.w[[ch]], t.w[[ch]], gobs[[ch]][[ind]], mutmat, maxmiss, phase.error.locs, 
		   initProb, PLOT=PLOT, minbg=min.bg, maxbg=max.bg)
    ans
  }
    for (ch in 1:nchrno)
    {
      orig.ll[[ch]]<-max.ll[[ch]]<-c.ll[[ch]]<-list()
      for (ind in 1:NUMI)
      {
	ch_ind=(ch-1)*NUMI+ind
	flips[[ind]][[ch]]<-tmp[[ch_ind]]$ind.max.flips
	max.ll[[ch]][[ind]]<-c.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.max.ll
	orig.ll[[ch]][[ind]]<-tmp[[ch_ind]]$ind.orig.ll
	nflips=nflips+tmp[[ch_ind]]$nflips
	niters=niters+tmp[[ch_ind]]$niters
      }
    }
}
rm(tmp)
cloglike<-sum(unlist(max.ll))
if (verbose)
  cat(nflips, " phase flips made after an average of ", niters/NUMI/nchrno, "hunts/ind/chromosome: log-likelihood", sum(unlist(orig.ll)), "-> ", cloglike, "\n")
if (LOG) 
{
  runtime<-as.numeric(Sys.time());diff.time<-runtime-old.runtime;old.runtime<-runtime;
  writelog(EMlogfile,"phasehunt",diff.time,len,Mu,rho,PI,alpha,lambda,theta,cloglike) 
}
