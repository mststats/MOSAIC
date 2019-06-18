# need to run this otherwise so that have approx. correct parameters to start, else donates usage is very poor
# just fit noanc on a couple of targets and a couple of chromosomes
fit_noanc_model=function(target, t.samp_chrnos, t.chrnos, t.NUMA, t.NUMP, t.kLL, t.A, t.KNOWN, t.label, t.NL, t.NN, t.umatch, t.G, t.dr, t.flips, t.gobs,
			 t.PI, t.Mu, t.rho, t.theta, t.alpha, t.lambda, t.prop.don, t.min.donors, t.max.donors, t.maxmatchsize, t.maxmatch, t.maxmiss, 
			 t.initProb, t.d.w, t.t.w, t.subNUMA, t.subNL, HPC, t.runtime, t.resultsdir, t.GpcM, t.eps, t.cloglike, doMu, doPI, dorho, dotheta, ffpath,
			 firstind, EM=T, getnoancgfbs=FALSE, t.LOG=T, t.HPC=2, verbose=TRUE) {
  # subNUMA=t.NUMA=>use all; number of target haps used in no-ancestry initial fit; don't use less than min(2,t.NUMA)
  ans=list()
  nchrno=length(t.chrnos)
  o.nchrno=nchrno;o.chrnos=t.chrnos;t.chrnos=t.samp_chrnos;nchrno=length(t.samp_chrnos);
  o.NUMA=t.NUMA;t.NUMA=min(o.NUMA,t.subNUMA);t.NUMI=max(t.NUMA/2,1)
  o.NUMP<-t.NUMP;o.label<-t.label;o.KNOWN<-t.KNOWN;o.NL<-t.NL;o.NN<-t.NN
  dons<-NULL;for (k in 1:t.kLL) dons<-c(dons,sort(sample(which(t.label==k),min(t.subNL,sum(t.label==k)))))
  t.NUMP<-length(dons);t.label<-c(t.label[dons],t.label[!t.KNOWN]);t.NL<-c(table(t.label));t.NN<-sum(t.NL);t.KNOWN<-c(t.KNOWN[dons],t.KNOWN[!t.KNOWN])
  # if required, use subset of targets and subset of donors on subset of chromosomes
  if (nchrno!=o.nchrno | t.NUMA!=o.NUMA | t.NUMP!=o.NUMP) 
  {
    o.umatch=t.umatch
    o.maxmatchsize=t.maxmatchsize
    o.G=t.G;o.flips=t.flips;o.gobs=t.gobs;
    t.G=t.G[match(t.samp_chrnos, o.chrnos)]; t.flips=t.gobs=list()
    for (ch in 1:nchrno)
    {
      t.umatch[[ch]]=t.umatch[[which(t.samp_chrnos[ch]==o.chrnos)]]
      t.maxmatchsize[ch]=o.maxmatchsize[which(t.samp_chrnos[ch]==o.chrnos)]
      t.d.w[[ch]]=t.d.w[[which(t.samp_chrnos[ch]==o.chrnos)]]
      t.t.w[[ch]]=t.t.w[[which(t.samp_chrnos[ch]==o.chrnos)]]
      if (t.NUMP!=o.NUMP)
	for (g in 1:t.G[ch])
	  t.d.w[[ch]]$w[[g]]=t.d.w[[ch]]$w[[g]][dons] # subset of the donors to fit the no latent ancestry model parameters 
      if (t.NUMA!=o.NUMA)
	for (g in 1:t.G[ch])
	  t.t.w[[ch]]$w[[g]]=t.t.w[[ch]]$w[[g]][1:t.NUMA] # subset of the targets to fit the no latent ancestry model parameters 
    }
    for (ind in 1:t.NUMI) {t.flips[[ind]]=list(); for (ch in 1:nchrno) t.flips[[ind]][[ch]]=o.flips[[ind]][[which(t.samp_chrnos[ch]==o.chrnos)]]}
    for (ch in 1:nchrno) {t.gobs[[ch]]=list(); for (ind in 1:t.NUMI) t.gobs[[ch]][[ind]]=o.gobs[[which(t.samp_chrnos[ch]==o.chrnos)]][[ind]]}
  }
  ###############################  fit no anc model and run EM #######################
  o.L<-t.A;t.A<-1
  noanc.rho=t.rho<-mean(t.rho)
  noanc.theta=t.theta<-mean(t.theta)
  # first compute the no-ancestry equivalent parameters t.Mu, t.rho, and theta. One for each ind.
  ind.Mu=list()
  for (ind in 1:t.NUMI) 
  {
    # use current ind specific parameters from ancestry aware model
    ind.Mu[[ind]]=matrix(rowSums(t(t(t.Mu)*t.alpha[[ind]])),t.kLL) # p(g) = sum_a(p(g|a)p(a))
  }
  noanc.Mu=t.Mu=matrix(rowSums(t.Mu%*%Reduce("+",t.alpha)/t.NUMI),t.kLL) 
  o.doMu<-doMu;doMu<-T;o.dotheta<-dotheta;dotheta<-T;o.dorho<-dorho;dorho<-T;
  o.doPI<-doPI;doPI<-F;o.PI<-t.PI;t.PI=list();for (ind in 1:t.NUMI) t.PI[[ind]]<-matrix(0,1,1)
  o.alpha=t.alpha;t.alpha=list();for (ind in 1:t.NUMI) t.alpha[[ind]]=1
  o.lambda=t.lambda;t.lambda=list();for (ind in 1:t.NUMI) t.lambda[[ind]]=0
  o.prop.don<-t.prop.don;o.max.donors<-t.max.donors
  t.prop.don=1;t.max.donors=t.NUMP # use all donor haplotypes here 
  transitions=list();for (ind in 1:t.NUMI) transitions[[ind]]<-s_trans(t.A,t.kLL,t.PI[[ind]],ind.Mu[[ind]],t.rho,t.NL)
  mutmat<-fmutmat(t.theta, t.A, t.maxmiss, t.maxmatch) # possibly overkill / some redundancy as t.maxmiss and t.maxmatch may have fallen for this subset
  # dummy run; this will return all donors at all gridpoints and is not affected by parameter values
  tmp=all_donates(target, t.A, t.NUMI, t.Mu, t.alpha, t.kLL, t.PI, t.rho, t.lambda, t.theta, verbose=verbose, t.get_switches=F, t.min.donors, t.max.donors, t.prop.don, t.NUMP, t.NL, t.G, 
		  t.umatch, t.maxmatchsize, t.maxmatch,t.maxmiss,t.d.w, t.t.w, t.gobs, t.flips, t.label, t.KNOWN, t.HPC, prethin=F, t.NUMA, nchrno, t.initProb, 
		  t.runtime, len, F, transitions, mutmat, NaN, NULL, ffpath)
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=t.runtime=tmp$runtime;cloglike=tmp$cloglike
  t.initProb=initprobs(T,t.NUMA,t.A,t.NUMP,t.kLL,t.PI,t.Mu,t.rho,t.alpha,t.label,t.NL)

  if(verbose) 
    cat("Fitting no-ancestry model\n") 
  t.runtime=old.runtime=tmp$rtime;diff.time=0
  if (t.LOG) 
  {
    tmp=create_logfile(t.resultsdir,target,t.kLL,t.A,t.NUMI,firstind,t.chrnos,nchrno,t.NN,t.GpcM)
    len=tmp$len
    noanclogfile=tmp$logfile
  }
  total=50 # only estimating some of the parameters, not required to be super accurate
  #stop("wait")
  if (EM) {
    # no anc fit and all donors included; should remove EM output
    tmp=run_EM(t.HPC, nchrno, t.PI, t.Mu, t.rho, t.theta, t.alpha, t.lambda, t.initProb, t.label, mutmat, transitions, ndonors, donates, donatesl, donatesr,
	       t.NUMA, t.NN, t.NL, t.NUMP, t.kLL, t.A, t.NUMI, t.max.donors, t.G, t.dr, t.gobs, t.maxmatchsize, t.umatch, t.flips, t.maxmatch, t.maxmiss, t.d.w, t.t.w,  
	       total, verbose=F, len, t.cloglike, t.LOG, noanclogfile, doPI, doMu, dotheta, dorho, TRUE, TRUE, TRUE, t.runtime, t.eps) 
    t.PI=tmp$t.PI;t.alpha=tmp$t.alpha;t.lambda=tmp$t.lambda;t.Mu=tmp$Mu;t.rho=tmp$rho;t.theta=tmp$theta;t.runtime=tmp$runtime;t.initProb=tmp$initProb;
    t.cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
  } 
  if (getnoancgfbs)
    noanc_gfbs=get_gfbs(t.NUMP, nchrno, t.max.donors, donates, donatesl, donatesr, t.NUMA, t.A, t.G, t.kLL, transitions, t.umatch, t.maxmatchsize, t.d.w, t.t.w, t.gobs, mutmat, 
			t.maxmiss, t.initProb, t.label, ndonors, t.flips, t.HPC)
  t.A<-o.L
  # return parameters, etc to correct sizes
  doMu<-o.doMu;dorho=o.dorho;dotheta=o.dotheta
  doPI<-o.doPI;t.PI<-o.PI;t.alpha=o.alpha;t.lambda=o.lambda
  noanc.rho=t.rho;t.rho<-rep(noanc.rho,t.A) # note that this includes all the latent ancestry switches
  noanc.Mu<-t.Mu;t.Mu<-NULL; for (l in 1:t.A) t.Mu<-cbind(t.Mu,noanc.Mu)
  noanc.theta=t.theta;t.theta<-rep(noanc.theta,t.A)
  # next line gets called if some groups dropped but it's fast so potential redundancy is ok
  for (ind in 1:t.NUMI) transitions[[ind]]<-s_trans(t.A,t.kLL,t.PI[[ind]],t.Mu,t.rho,t.NL)
  mutmat<-fmutmat(t.theta, t.A, t.maxmiss, t.maxmatch)
  if (!getnoancgfbs)
  ans$transitions=transitions
  ans$mutmat=mutmat
  ans$Mu=t.Mu
  ans$theta=t.theta
  ans$rho=t.rho
  ans$ndonors=ndonors
  ans$donates=donates
  ans$donatesl=donatesl
  ans$donatesr=donatesr
  if (getnoancgfbs)
    ans$noanc_gfbs=noanc_gfbs
  return(ans)
}
