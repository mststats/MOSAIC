# functions to calculate the locally specified list of top (most useful) donors to pass to the HMM thus reducing computation time and memory usage
# if we have a version impervious to phasing it will last longer without changing during thin / phase / EM
# code to select enough donors to capture prop.don of the copying 
create_donates<-function(getswitches,ch,ind,t.NUMA,t.umatch,t.maxmatchsize,t.maxmatch,t.maxmiss,t.d.w,t.t.w,t.gobs,t.flips,t.kLL,t.Mu,t.rho,t.theta,t.HPC,prethin=FALSE,
			 t.min.donors,t.max.donors,t.prop.don,t.NUMP,t.NL,t.label,t.G_ch,
			 prethin_ndonors=NULL, prethin_donates=NULL, prethin_donatesl=NULL, prethin_donatesr=NULL)
{
  THIN=F
  if (t.NUMA==1) H=1 else H=2
  hap<-c(ind*2-1,ind*2)
  # we need these (temporarily) to calculate E[switches] and (bizarrely) to calculate the thinned version of themselves
  t.donates<-0:(t.NUMP-1) # Holds vectors to indicate which donors are copied from at each gridpoint. 
  t.donatesl<-t.donatesr<-t.donates # left and right one locations of aligned indices
  t.ndonors<-rep(t.NUMP,t.G_ch) # too large; should limit this really
  if (t.max.donors==t.NUMP & !getswitches) 
    return(list(ndonors=t.ndonors,donates=t.donates,donatesl=t.donatesl,donatesr=t.donatesr)) # these are vectors; dealt with automatically in cpp calls
  # note that we always need this when thinning; see probmass below
  probmass<-matrix(0,t.G_ch,t.NUMP)
  switches<-list() # returns empty list if not needed; serves as placeholder
  # fit a HMM with no latent ancestry
  A<-1;tmpPI<-matrix(0,1,1);
  noanc_initProb<-matrix(0,H,t.kLL); for (h in 1:H) noanc_initProb[h,]<-t.Mu/t.NL[1:t.kLL]
  noanc_mutmat<-fmutmat(t.theta, A, t.maxmiss, t.maxmatch)
  noanc_transitions<-s_trans(A,t.kLL,tmpPI,t.Mu,t.rho,t.NL)
  noanc_fors<-noanc_sumfors<-noanc_backs<-noanc_scalefactor<-noanc_scalefactorb<-list()
  for (h in 1:H) { 
    noanc_fors[[h]]<-rep(0,t.G_ch*t.NUMP)
    noanc_sumfors[[h]]<-matrix(0,t.G_ch,A)
    noanc_backs[[h]]<-rep(0,t.G_ch*t.NUMP)
    noanc_scalefactor[[h]]<-rep(0,t.G_ch)
    noanc_scalefactorb[[h]]<-rep(0,t.G_ch)
  }
  if (prethin & !getswitches) # only need to use this if pre-thinning and not finding switch counts based on all donors
  {
    THIN=T
    t.ndonors=prethin_ndonors[[ch]][[ind]]
    if (t.HPC)
    {
      open(prethin_donates[[ch]][[ind]]);open(prethin_donatesl[[ch]][[ind]]);open(prethin_donatesr[[ch]][[ind]]);
      t.donates=prethin_donates[[ch]][[ind]][];t.donatesl=prethin_donatesl[[ch]][[ind]][];t.donatesr=prethin_donatesr[[ch]][[ind]][]
      close(prethin_donates[[ch]][[ind]]);close(prethin_donatesl[[ch]][[ind]]);close(prethin_donatesr[[ch]][[ind]]);
    }
    if (!t.HPC)
    {
      t.donates=prethin_donates[[ch]][[ind]];t.donatesl=prethin_donatesl[[ch]][[ind]];t.donatesr=prethin_donatesr[[ch]][[ind]]
    }
  }
  for (h in 1:H)
  {
    k=hap[h]
    cppforward(k,t.NUMA,t.NUMP,THIN,t.NUMP,t.kLL,A,0,t.G_ch,t.G_ch,noanc_transitions,t.umatch,t.maxmatchsize,t.d.w,t.t.w,t.gobs,noanc_mutmat,t.maxmiss,noanc_initProb[h,],
	       t.label,t.ndonors,t.donates,t.donatesl,t.flips,noanc_fors[[h]],noanc_sumfors[[h]],noanc_scalefactor[[h]])
    cppbackward(k,t.NUMA,t.NUMP,THIN,t.NUMP,A,0,t.G_ch,t.G_ch,noanc_transitions,t.umatch,t.maxmatchsize,t.d.w,t.t.w,t.gobs,noanc_mutmat,t.maxmiss,
		t.label,t.ndonors,t.donates,t.donatesr,t.flips,noanc_backs[[h]],noanc_scalefactorb[[h]])
    probmass<-probmass+cppforback(t.NUMP,THIN,t.NUMP,A,t.G_ch,t.ndonors,t.donates,noanc_fors[[h]],noanc_backs[[h]]) # 1 as first argument b/c using all now
    if (getswitches) 
      switches[[h]]<-t(matrix(cppswitches(h,t.NUMA,t.NUMP,THIN,t.NUMP,t.G_ch,t.NL,t.label,noanc_sumfors[[h]],noanc_backs[[h]],
					  noanc_transitions,t.flips,noanc_mutmat,t.maxmiss,t.umatch,t.maxmatchsize,t.d.w,t.t.w,t.gobs,t.ndonors,t.donates)$switches,t.NUMP))
  }
  if (t.max.donors==t.NUMP & getswitches) 
    return(list(ndonors=t.ndonors,donates=t.donates,donatesl=t.donatesl,donatesr=t.donatesr,switches=switches)) 
  #print(range(t.ndonors));print(range(noanc_fors[[1]]));readline()
  # next lines tested and absolutely required
  if (t.NUMA>1)
    for (h in 1:2) # both haps fors and backs must be calculated before the next line is called
      probmass<-probmass+cppforback(t.NUMP,THIN,t.NUMP,A,t.G_ch,t.ndonors,t.donates,noanc_fors[[h]],noanc_backs[[h+ifelse(h%%2,1,-1)]]) # 1 as first argument b/c using all now
  # if t.NUMA==1 probmass is now f[1]*b[1] 
  # if t.NUMA>1 probmass is now f[1]*b[1], f[2]*b[2], f[1]*b[2], f[2]*b[1]
  probmass<-t(t(probmass)/rowSums(probmass))
  f<-function(x)
  {
    tmp<-order(x,decreasing=T)[1:t.max.donors] # descending order of donors; can't use partial sorting as index can't be returned
    quants<-cumsum(x[tmp]);quants[t.max.donors]=1
    donors<-tmp
    cutoff<-max(which(quants>t.prop.don)[1], t.min.donors) # index of first one to go over t.prop.don but take at least t.min.donors
    donors<-donors[1:t.max.donors] # easier to return also the not needed ones
    #donors=1:t.max.donors;cutoff=t.max.donors;# debugging
    return(list(donors=donors, ndonors=cutoff))
  }
  g.donors<-apply(probmass,1,f) 
  t.ndonors<-vapply(g.donors,function(x) x$ndonors,FUN.VALUE=0L)
  t.donates<-matrix(vapply(g.donors,function(x) x$donors,FUN.VALUE=rep(0L,t.max.donors)),t.max.donors)
  # now calculate the right shifted and left shifted locations of the matching indices 
  t.donatesl<-t.donatesr<-matrix(0L,t.max.donors,t.G_ch) # re-sizing and 0s first
  t.donatesl[,2:t.G_ch]<-matrix(unlist(vapply(2:t.G_ch, function(g) match(t.donates[,g], t.donates[,g-1]),
					     FUN.VALUE=rep(0L,t.max.donors)),use.names=F),t.max.donors)
  t.donatesr[,1:(t.G_ch-1)]<-matrix(unlist(vapply(1:(t.G_ch-1), function(g) match(t.donates[,g], t.donates[,g+1]),
						 FUN.VALUE=rep(0L,t.max.donors)),use.names=F),t.max.donors)
  t.donatesl[is.na(t.donatesl)]<-0 # match() returns NA where no match occurs, replace with 0s
  t.donatesr[is.na(t.donatesr)]<-0 # match() returns NA where no match occurs, replace with 0s
  return(list(ndonors=t.ndonors,donates=t.donates-1,donatesl=t.donatesl-1,donatesr=t.donatesr-1,switches=switches)) # -1 to convert to cpp indexing
}
# functions to load up the donors already created by create_donates etc
getdonates<-function(t.donates,t.NUMI) # note that this is used for donates, donatesl, and donatesr
{
  donates_chr=list()
  for(ind in 1:t.NUMI) 
  {
    open(t.donates[[ind]])
    donates_chr[[ind]]=t.donates[[ind]][] # open all donors on this chromosome
    close(t.donates[[ind]])
  }
  donates_chr
}
getdonates_ind<-function(t.donates) # note that this is used for donates, donatesl, and donatesr
{
  open(t.donates)
  donates_chr_ind=t.donates[] # open all donors on this chromosome
  close(t.donates)
  donates_chr_ind
}


# function to estimate the most useful donors for all target admixed genomes; other donors will not be used in the HMM
# note that this is done locally to each gridpoint allowing for changing top donors along each target genome
# first compute the no-ancestry equivalent parameters Mu, rho, and theta. One for each ind.
all_donates=function(target, t.A, t.NUMI, t.Mu, t.alpha, t.kLL, t.PI, t.rho, t.lambda, t.theta, verbose=T, t.get_switches, t.min.donors, t.max.donors, t.prop.don, t.NUMP, 
		     t.NL, t.G, t.umatch, t.maxmatchsize, t.maxmatch, t.maxmiss, t.d.w, t.t.w, t.gobs, t.flips, t.label, t.KNOWN, t.HPC, prethin=F, t.NUMA, 
		     t.nchrno, t.initProb, t.old.runtime, t.len, t.LOG, t.transitions, t.mutmat, t.cloglike, t.logfile, ffpath) {
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
	  tmp2=create_donates(t.get_switches,ch,ind,t.NUMA,t.umatch[[ch]],t.maxmatchsize[ch],t.maxmatch,t.maxmiss,t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],t.flips[[ind]][[ch]],
			      t.kLL,ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],t.HPC,prethin=prethin,t.min.donors,t.max.donors,t.prop.don,t.NUMP,t.NL,t.label,t.G[ch]) 
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
	  ans=create_donates(t.get_switches,ch,ind,t.NUMA,t.umatch[[ch]],t.maxmatchsize[ch],t.maxmatch,t.maxmiss,t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],t.flips[[ind]][[ch]],
			     t.kLL,ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],t.HPC,prethin=prethin,t.min.donors,t.max.donors,t.prop.don,t.NUMP,t.NL,t.label,t.G[ch])
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
      tmp2=create_donates(t.get_switches,ch,ind,t.NUMA,t.umatch[[ch]],t.maxmatchsize[ch],t.maxmatch,t.maxmiss,t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],t.flips[[ind]][[ch]],t.kLL,
			  ind.Mu[[ind]],ind.rho[[ind]],ind.theta[[ind]],t.HPC,prethin=prethin,t.min.donors,t.max.donors,t.prop.don,t.NUMP,t.NL,t.label,t.G[ch])
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
      cat(": log-likelihood", t.cloglike, "-> ")
    # some overhead in this so only run if asked for i.e. t.LOG==TRUE
    cloglike=get_loglike(t.NUMA, t.nchrno, t.G, t.A, t.kLL, t.max.donors, t.NUMP, ndonors, donates, donatesl, t.transitions, t.maxmatchsize, t.umatch, t.flips, t.mutmat, t.maxmiss, t.initProb,t.d.w,t.t.w,t.gobs,t.label, t.HPC)
    writelog(t.logfile,"thinning",diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,cloglike) 
    if (verbose & !t.get_switches & t.max.donors<t.NUMP) 
      cat(cloglike)
  } else cloglike=NaN
  if (verbose) cat("\n")
  ans=list()
  ans$ndonors=ndonors
  ans$donates=donates
  ans$donatesl=donatesl
  ans$donatesr=donatesr
  ans$runtime=runtime
  ans$cloglike=cloglike
  if (t.get_switches) 
    ans$noanc_gswitches=noanc_gswitches
  return(ans)
}

