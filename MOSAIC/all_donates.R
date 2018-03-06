# script to estimate the most useful donors for all target admixed genomes; other donors will not be used in the HMM
# note that this is done locally to each gridpoint allowing for changing top donors along each target genome
# first compute the no-ancestry equivalent parameters Mu, rho, and theta. One for each ind.
ind.Mu=ind.rho=ind.theta=list()
for (ind in 1:NUMI) 
{
  # use current ind specific parameters from ancestry aware model
  ind.Mu[[ind]]=matrix(rowSums(t(t(Mu)*alpha[[ind]])),kLL) # p(g) = sum_a(p(g|a)p(a))
  # note that if commonrho=T then this will just be rho 
  ind.rho[[ind]]=rho%*%alpha[[ind]]+sum(diag(PI[[ind]])) # include all ancestry self-switches as these impose a hap switch (zero if absorbrho)
  #tmp.d=-log(1-ind.rho[[ind]])*max.donors;ind.rho[[ind]]=1-exp(-tmp.d/NUMP) # ~= max.donors/NUMP for small rho 
  # note that if commontheta=T then this will just be theta
  ind.theta[[ind]]=theta%*%alpha[[ind]]
}
require(parallel)
if (verbose & !get_switches & max.donors<NUMP)
  cat("thinning to at most", max.donors, "donors at each gridpoint")
if (NUMA==1) {H=1} else H=2
ndonors<-donates<-donatesl<-donatesr<-list()
if (get_switches) noanc_gswitches<-list()
for (ch in 1:nchrno)
{
  ndonors[[ch]]<-list()
  donates[[ch]]<-donatesl[[ch]]<-donatesr[[ch]]<-list()
  if (get_switches) noanc_gswitches[[ch]]<-array(0L,c(kLL,G[ch],NUMA))
}
if (HPC!=2)
{
  for (ch in 1:nchrno)
  {
    # if using all then a vector will suffice for each donates, donatesl, donatesr rather than a matrix with G[ch] columns
    NvecsG=ifelse(max.donors==NUMP, 1, G[ch]) 
    if (HPC==1)
    {
      tmp<-foreach(ind=1:NUMI) %dopar%
      {
	tmp2=create_donates(get_switches,ch,ind,umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],flips[[ind]][[ch]],kLL,ind.Mu[[ind]],
			    ind.rho[[ind]],ind.theta[[ind]],HPC,prethin=prethin) 
	ans_ndonors=tmp2$ndonors
	ans_donates=ff(tmp2$donates,vmode="integer",dim=c(max.donors,NvecsG),filename=paste0(ffpath,target,"_donates_",ch,"_",ind,".ff"),overwrite=T);close(ans_donates)
	ans_donatesl=ff(tmp2$donatesl,vmode="integer",dim=c(max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesl_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesl)
	ans_donatesr=ff(tmp2$donatesr,vmode="integer",dim=c(max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesr_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesr)
	ans_switches=list()
	if (get_switches)
	  for (h in 1:H)
	  {
	    ans_switches[[h]]=list()
	    for (i in 1:kLL)
	      ans_switches[[h]][[i]]=rowSums(tmp2$switches[[h]][,label[KNOWN]==i]) # sum switches over panels; sum to 1 at each gridpoint for each target hap
	  }
	rm(tmp2)
	gc()
	list(ndonors=ans_ndonors,donates=ans_donates,donatesl=ans_donatesl,donatesr=ans_donatesr,switches=ans_switches)
      }
    }
    if (!HPC)
    {
      tmp<-foreach(ind=1:NUMI) %dopar%
      {
	ans=create_donates(get_switches,ch,ind,umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],flips[[ind]][[ch]],kLL,ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],HPC,prethin=prethin)
	if (get_switches)
	{
	  tmpswitches=ans$switches
	  for (h in 1:H)
	  {
	    ans$switches[[h]]=list()
	    for (i in 1:kLL)
	      ans$switches[[h]][[i]]=rowSums(tmpswitches[[h]][,label[KNOWN]==i]) # sum switches over panels; sum to 1 at each gridpoint for each target hap
	  }
	}
	ans
      }
    }
    for (ind in 1:NUMI)
    {
      ndonors[[ch]][[ind]]<-tmp[[ind]]$ndonors
      donates[[ch]][[ind]]<-tmp[[ind]]$donates
      donatesl[[ch]][[ind]]<-tmp[[ind]]$donatesl
      donatesr[[ch]][[ind]]<-tmp[[ind]]$donatesr
      if (get_switches) 
      {
	if (NUMA>1) hap<-c(ind*2-1,ind*2)
	if (NUMA==1) hap<-1
	for (h in 1:H)
	{
	  if (!HPC)
	    for (i in 1:kLL)
	      noanc_gswitches[[ch]][i,,hap[h]]<-tmp[[ind]]$switches[[h]][[i]]
	  if (HPC==1)
	  {
	    for (i in 1:kLL)
	      noanc_gswitches[[ch]][i,,hap[h]]<-tmp[[ind]]$switches[[h]][[i]]
	  }
	}
      }
    }
    rm(tmp)
  } 
}
if (HPC==2)
{
  tmp<-foreach(ch_ind=1:(nchrno*NUMI)) %dopar% 
  {
    ch=as.integer((ch_ind-0.5)/NUMI)+1
    ind=(ch_ind-1)%%NUMI+1
    NvecsG=ifelse(max.donors==NUMP, 1, G[ch]) 
    tmp2=create_donates(get_switches,ch,ind,umatch[[ch]],maxmatchsize[ch],d.w[[ch]],t.w[[ch]],gobs[[ch]][[ind]],flips[[ind]][[ch]],kLL,ind.Mu[[ind]],
			ind.rho[[ind]],ind.theta[[ind]],HPC,prethin=prethin)
    ans_ndonors=tmp2$ndonors
    ans_donates=ff(tmp2$donates,vmode="integer",dim=c(max.donors,NvecsG),filename=paste0(ffpath,target,"_donates_",ch,"_",ind,".ff"),overwrite=T);close(ans_donates)
    ans_donatesl=ff(tmp2$donatesl,vmode="integer",dim=c(max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesl_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesl)
    ans_donatesr=ff(tmp2$donatesr,vmode="integer",dim=c(max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesr_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesr)
    ans_switches=list()
    if (get_switches)
      for (h in 1:H)
      {
	ans_switches[[h]]=list()
	for (i in 1:kLL)
	  ans_switches[[h]][[i]]=rowSums(tmp2$switches[[h]][,label[KNOWN]==i]) # sum switches over panels; sum to 1 at each gridpoint for each target hap
      }
    rm(tmp2)
    gc()
    list(ndonors=ans_ndonors,donates=ans_donates,donatesl=ans_donatesl,donatesr=ans_donatesr,switches=ans_switches)
  }
    for (ch in 1:nchrno)
    {
      for (ind in 1:NUMI)
      {
	ch_ind=(ch-1)*NUMI+ind
	ndonors[[ch]][[ind]]<-tmp[[ch_ind]]$ndonors
	donates[[ch]][[ind]]<-tmp[[ch_ind]]$donates
	donatesl[[ch]][[ind]]<-tmp[[ch_ind]]$donatesl
	donatesr[[ch]][[ind]]<-tmp[[ch_ind]]$donatesr
	close(donates[[ch]][[ind]]);close(donatesl[[ch]][[ind]]);close(donatesr[[ch]][[ind]])
	if (get_switches) 
	{
	  if (NUMA>1) hap<-c(ind*2-1,ind*2)
	  if (NUMA==1) hap<-1
	  for (h in 1:H)
	  {
	    for (i in 1:kLL)
	      noanc_gswitches[[ch]][i,,hap[h]]<-tmp[[ch_ind]]$switches[[h]][[i]] # sum switches over panels; sum to 1 at each gridpoint for each target hap
	  }
	}
      }
    }
  rm(tmp)
}
if (LOG)
{
  runtime<-as.numeric(Sys.time());diff.time<-runtime-old.runtime;old.runtime<-runtime;
  if (verbose & !get_switches & max.donors<NUMP) 
    cat(": log-likelihood", cloglike, "-> ")
  source("klikelihood.R") # some overhead in this so only run if asked for i.e. LOG=T
  writelog("thinning")
  if (verbose & !get_switches & max.donors<NUMP) 
    cat(cloglike)
} else cloglike=NaN
if (verbose) cat("\n")

