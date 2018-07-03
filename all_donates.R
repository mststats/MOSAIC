# function to estimate the most useful donors for all target admixed genomes; other donors will not be used in the HMM
# note that this is done locally to each gridpoint allowing for changing top donors along each target genome
# first compute the no-ancestry equivalent parameters Mu, rho, and theta. One for each ind.
all_donates=function(t.NUMI, t.Mu, t.alpha, t.kLL, t.PI, t.rho, t.lambda, t.theta, verbose=T, t.get_switches, t.max.donors, t.NUMP, t.G, t.umatch, t.maxmatchsize, t.d.w, 
		     t.t.w, t.gobs, t.flips, t.label, t.KNOWN, t.HPC, prethin=F, t.NUMA, t.nchrno, t.initProb, t.old.runtime, t.len, t.LOG, t.transitions, t.mutmat) {
  ind.Mu=ind.rho=ind.theta=list()
  for (ind in 1:t.NUMI) 
  {
    # use current ind specific parameters from ancestry aware model
    ind.Mu[[ind]]=matrix(rowSums(t(t(t.Mu)*t.alpha[[ind]])),t.kLL) # p(g) = sum_a(p(g|a)p(a))
    # note that if commonrho=T then this will just be t.rho 
    ind.rho[[ind]]=t.rho%*%t.alpha[[ind]]+sum(diag(t.PI[[ind]])) # include all ancestry self-switches as these impose a hap switch (zero if absorbrho)
    #tmp.d=-log(1-ind.rho[[ind]])*t.max.donors;ind.rho[[ind]]=1-exp(-tmp.d/t.NUMP) # ~= t.max.donors/t.NUMP for small t.rho 
    # note that if commontheta=T then this will just be t.theta
    ind.theta[[ind]]=t.theta%*%t.alpha[[ind]]
  }
  if (verbose & !t.get_switches & t.max.donors<t.NUMP)
    cat("thinning to at most", t.max.donors, "donors at each gridpoint")
  if (t.NUMA==1) {H=1} else H=2
  ndonors<-donates<-donatesl<-donatesr<-list()
  if (t.get_switches) noanc_gswitches<-list()
  for (ch in 1:t.nchrno)
  {
    ndonors[[ch]]<-list()
    donates[[ch]]<-donatesl[[ch]]<-donatesr[[ch]]<-list()
    if (t.get_switches) noanc_gswitches[[ch]]<-array(0L,c(t.kLL,t.G[ch],t.NUMA))
  }
  if (t.HPC!=2)
  {
    for (ch in 1:t.nchrno)
    {
      # if using all then a vector will suffice for each donates, donatesl, donatesr rather than a matrix with t.G[ch] columns
      NvecsG=ifelse(t.max.donors==t.NUMP, 1, t.G[ch]) 
      if (t.HPC==1)
      {
	tmp<-foreach(ind=1:t.NUMI) %dopar%
	{
	  tmp2=create_donates(t.get_switches,ch,ind,t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],t.flips[[ind]][[ch]],
			      t.kLL,ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],t.HPC,prethin=prethin,t.max.donors,t.NUMP) 
	  ans_ndonors=tmp2$ndonors
	  ans_donates=ff(tmp2$donates,vmode="integer",dim=c(t.max.donors,NvecsG),filename=paste0(ffpath,target,"_donates_",ch,"_",ind,".ff"),overwrite=T)
	  close(ans_donates)
	  ans_donatesl=ff(tmp2$donatesl,vmode="integer",dim=c(t.max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesl_",ch,"_",ind,".ff"),overwrite=T)
	  close(ans_donatesl)
	  ans_donatesr=ff(tmp2$donatesr,vmode="integer",dim=c(t.max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesr_",ch,"_",ind,".ff"),overwrite=T)
	  close(ans_donatesr)
	  ans_switches=list()
	  if (t.get_switches)
	    for (h in 1:H)
	    {
	      ans_switches[[h]]=list()
	      for (i in 1:t.kLL)
		# sum switches over panels; sum to 1 at each gridpoint for each target hap
		ans_switches[[h]][[i]]=rowSums(tmp2$switches[[h]][,t.label[t.KNOWN]==i]) 
	    }
	  rm(tmp2)
	  gc()
	  list(ndonors=ans_ndonors,donates=ans_donates,donatesl=ans_donatesl,donatesr=ans_donatesr,switches=ans_switches)
	}
      }
      if (!t.HPC)
      {
	tmp<-foreach(ind=1:t.NUMI) %dopar%
	{
	  ans=create_donates(t.get_switches,ch,ind,t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],t.flips[[ind]][[ch]],
			     t.kLL,ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],t.HPC,prethin=prethin,t.max.donors,t.NUMP)
	  if (t.get_switches)
	  {
	    tmpswitches=ans$switches
	    for (h in 1:H)
	    {
	      ans$switches[[h]]=list()
	      for (i in 1:t.kLL)
		# sum switches over panels; sum to 1 at each gridpoint for each target hap
		ans$switches[[h]][[i]]=rowSums(tmpswitches[[h]][,t.label[t.KNOWN]==i]) 
	    }
	  }
	  ans
	}
      }
      for (ind in 1:t.NUMI)
      {
	ndonors[[ch]][[ind]]<-tmp[[ind]]$ndonors
	donates[[ch]][[ind]]<-tmp[[ind]]$donates
	donatesl[[ch]][[ind]]<-tmp[[ind]]$donatesl
	donatesr[[ch]][[ind]]<-tmp[[ind]]$donatesr
	if (t.get_switches) 
	{
	  if (t.NUMA>1) hap<-c(ind*2-1,ind*2)
	  if (t.NUMA==1) hap<-1
	  for (h in 1:H)
	  {
	    if (!t.HPC)
	      for (i in 1:t.kLL)
		noanc_gswitches[[ch]][i,,hap[h]]<-tmp[[ind]]$switches[[h]][[i]]
	    if (t.HPC==1)
	    {
	      for (i in 1:t.kLL)
		noanc_gswitches[[ch]][i,,hap[h]]<-tmp[[ind]]$switches[[h]][[i]]
	    }
	  }
	}
      }
      rm(tmp)
    } 
  }
  if (t.HPC==2)
  {
    tmp<-foreach(ch_ind=1:(t.nchrno*t.NUMI)) %dopar% 
    {
      ch=as.integer((ch_ind-0.5)/t.NUMI)+1
      ind=(ch_ind-1)%%t.NUMI+1
      NvecsG=ifelse(t.max.donors==t.NUMP, 1, t.G[ch]) 
      tmp2=create_donates(t.get_switches,ch,ind,t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],t.flips[[ind]][[ch]],t.kLL,
			  ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],t.HPC,prethin=prethin,t.max.donors,t.NUMP)
      ans_ndonors=tmp2$ndonors
      ans_donates=ff(tmp2$donates,vmode="integer",dim=c(t.max.donors,NvecsG),filename=paste0(ffpath,target,"_donates_",ch,"_",ind,".ff"),overwrite=T)
      close(ans_donates)
      ans_donatesl=ff(tmp2$donatesl,vmode="integer",dim=c(t.max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesl_",ch,"_",ind,".ff"),overwrite=T)
      close(ans_donatesl)
      ans_donatesr=ff(tmp2$donatesr,vmode="integer",dim=c(t.max.donors,NvecsG),filename=paste0(ffpath,target,"_donatesr_",ch,"_",ind,".ff"),overwrite=T)
      close(ans_donatesr)
      ans_switches=list()
      if (t.get_switches)
	for (h in 1:H)
	{
	  ans_switches[[h]]=list()
	  for (i in 1:t.kLL)
	    # sum switches over panels; sum to 1 at each gridpoint for each target hap
	    ans_switches[[h]][[i]]=rowSums(tmp2$switches[[h]][,t.label[t.KNOWN]==i]) 
	}
      rm(tmp2)
      gc()
      list(ndonors=ans_ndonors,donates=ans_donates,donatesl=ans_donatesl,donatesr=ans_donatesr,switches=ans_switches)
    }
      for (ch in 1:t.nchrno)
      {
	for (ind in 1:t.NUMI)
	{
	  ch_ind=(ch-1)*t.NUMI+ind
	  ndonors[[ch]][[ind]]<-tmp[[ch_ind]]$ndonors
	  donates[[ch]][[ind]]<-tmp[[ch_ind]]$donates
	  donatesl[[ch]][[ind]]<-tmp[[ch_ind]]$donatesl
	  donatesr[[ch]][[ind]]<-tmp[[ch_ind]]$donatesr
	  close(donates[[ch]][[ind]]);close(donatesl[[ch]][[ind]]);close(donatesr[[ch]][[ind]])
	  if (t.get_switches) 
	  {
	    if (t.NUMA>1) hap<-c(ind*2-1,ind*2)
	    if (t.NUMA==1) hap<-1
	    for (h in 1:H)
	    {
	      for (i in 1:t.kLL)
		# sum switches over panels; sum to 1 at each gridpoint for each target hap
		noanc_gswitches[[ch]][i,,hap[h]]<-tmp[[ch_ind]]$switches[[h]][[i]] 
	    }
	  }
	}
      }
    rm(tmp)
  }
  runtime<-as.numeric(Sys.time());diff.time<-runtime-t.old.runtime; # even if not required i.e. NaN
  if (t.LOG)
  {
    if (verbose & !t.get_switches & t.max.donors<t.NUMP) 
      cat(": log-likelihood", cloglike, "-> ")
    # some overhead in this so only run if asked for i.e. t.LOG==TRUE
    cloglike=get_loglike(t.NUMA, t.nchrno, t.G, L, t.kLL, t.max.donors, t.NUMP, ndonors, donates, donatesl, t.transitions, t.maxmatchsize, t.umatch, t.flips, t.mutmat, maxmiss, t.initProb)
    writelog(EMlogfile,"thinning",diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,cloglike) 
    if (verbose & !t.get_switches & t.max.donors<t.NUMP) 
      cat(cloglike)
  } else cloglike=NaN
  if (verbose) cat("\n")
  if (!t.get_switches) return(list(ndonors=ndonors, donates=donates, donatesl=donatesl, donatesr=donatesr, runtime=runtime, cloglike=cloglike)) 
  if (t.get_switches) return(list(ndonors=ndonors, donates=donates, donatesl=donatesl, donatesr=donatesr, noanc_gswitches=noanc_gswitches, runtime=runtime, cloglike=cloglike)) 
}

