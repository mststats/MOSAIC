run_mosaic=function(ANC,chrnos,datasource,doMu,doPI,dorho,dotheta,EM,ffpath,firstind,L,MC,nchrno,NUMA,PLOT,target,verbose) {
  ans=list()
  tmp=setup_data_etc(NUMA,target,nchrno) # sets default parameters, sets up some objects required later, reads in data, and initialises model.
  # should replace the below with assign() usage
  resultsdir=tmp$resultsdir;PHASE=tmp$PHASE;HPC=tmp$HPC;GpcM=tmp$GpcM;LOG=tmp$LOG
  mcmcprog=tmp$mcmcprog;absorbrho=tmp$absorbrho;commonrho=tmp$commonrho;commontheta=tmp$commontheta;prethin=tmp$prethin
  s.M=tmp$s.M;M=tmp$M;PI.total=tmp$PI.total;s.total=tmp$s.total;REPS=tmp$REPS
  eps.lower=tmp$eps.lower;min.bg=tmp$min.bg;max.bg=tmp$max.bg;samp_chrnos=tmp$samp_chrnos;dr=tmp$dr
  FLAT=tmp$FLAT;maxmatch=tmp$maxmatch;maxmiss=tmp$maxmiss;umatch=tmp$umatch;d.w=tmp$d.w
  t.w=tmp$t.w;ans$g.loc=tmp$g.loc;gobs=tmp$gobs;NUMP=tmp$NUMP;NUMI=tmp$NUMI
  label=tmp$label;KNOWN=tmp$KNOWN;kLL=tmp$kLL;NL=tmp$NL;G=tmp$G
  NN=tmp$NN;maxmatchsize=tmp$maxmatchsize;panels=tmp$panels;min.donors=tmp$min.donors;
  theta=tmp$theta;rho=tmp$rho;lambda=tmp$lambda;alpha=tmp$alpha;PI=tmp$PI;Mu=tmp$MU
  transitions=tmp$transitions;total=tmp$total
  o.PI=tmp$o.PI;o.alpha=tmp$o.alpha;o.rho=tmp$o.rho;o.lambda=tmp$o.lambda;o.theta=tmp$o.theta;o.phi.theta=tmp$o.phi.theta 
  flips=tmp$flips;prop.don=tmp$prop.don;max.donors=tmp$max.donors
  if (target=="simulated")
    ans$g.true_anc=tmp$g.true_anc
  rm(tmp)
  old.runtime<-as.numeric(Sys.time())
  o.total=total
  writelog<-function(t.logfile,t.alg,t.diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,t.cloglike) # single consistent function to write to EMlogfile
    write(file=t.logfile,c(t.alg,signif(t.diff.time,4),signif(t(t.Mu),4),signif(t.rho,4),c(sapply(t.PI, function(x) signif(t(x),4))),
			   sapply(t.alpha, function(x) signif(x,4)),sapply(t.lambda,function(x) round(x,4)),signif(t.theta,4),round(t.cloglike,4)),ncol=t.len,append=T)
  eps=log(1.01) # i.e. a 1% increase in relative likelihood
  Mu<-matrix(rep(1/kLL,L*kLL),kLL);for (ind in 1:NUMI) alpha[[ind]]=rep(1/L,L) # flatten out w.r.t. ancestry
  runtime=NaN
  # always need to run noanc.R b/c need good paras for init_Mu
  tmp=fit_noanc_model(samp_chrnos, chrnos, NUMA, NUMP, kLL, L, KNOWN, label, NL, NN, umatch, G, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
		      prop.don, max.donors, maxmatch, maxmiss, initProb, d.w, t.w, NUMA, 100, HPC) 
  transitions=tmp$t.transitions;mutmat=tmp$mutmat;Mu=tmp$Mu;theta=tmp$theta;rho=tmp$rho
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
  initProb=initprobs(T,NUMA,L,NUMP,kLL,PI,Mu,rho,alpha,label,NL)
  runtime<-as.numeric(Sys.time())
  if (kLL>L) # otherwise can't cluster kLL things into L clusters
  {
    # use this to get #switches in noanc model w/o writelog
    tmp=all_donates(NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=T, max.donors, NUMP, G, umatch, maxmatchsize, d.w, 
		    t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len,F,transitions,mutmat)
    ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
    noanc_gswitches=tmp$noanc_gswitches
    rm(tmp)
    windowed_copying<-window_chunks(nswitches=noanc_gswitches,ww=0.5,verbose=verbose) # similar in the same windows
    rm(noanc_gswitches) 
    tmp<-cluster_windows(windowed_copying,t.L=L,verbose=verbose)
    Mu<-tmp$Mu
    alpha<-tmp$alpha
    PI<-tmp$PI
    all.o.lambda=tmp$lambda
    #lambda=as.list(sapply(tmp$lambda,mean,na.rm=T)) # doesn't work well!
    lambda=o.lambda # re-use o.lambda from above
    rm(windowed_copying,tmp)
  } else {diag(Mu)=10*diag(Mu);Mu=t(t(Mu)/colSums(Mu))}
  rownames(Mu)<-panels[1:kLL]
  o.Mu<-Mu;o.alpha<-alpha;o.lambda=lambda;o.PI=PI # these are the official starting ancestry related parameters now
  mutmat<-fmutmat(theta, L, maxmiss, maxmatch); for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,PI[[ind]],Mu,rho,NL)
  o.M<-M;M<-s.M
  if (EM) {
    tmp=create_logfile(resultsdir,target,kLL,L,NUMI,firstind,chrnos,nchrno,NN,GpcM)
    runtime=old.runtime=tmp$rtime;diff.time=0;len=tmp$len
    EMlogfile=tmp$logfile
  }
  cloglike=NaN 
  # decide on donor set using initial parameters
  tmp=all_donates(NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=F, max.donors, NUMP, G, umatch, maxmatchsize, d.w, 
		  t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len, F, transitions, mutmat)
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike

  ############ a few PI only updates first; very useful to do before first re-phasing
  if (PI.total>0 & EM)
  {
    o.doMu=doMu;o.dotheta=dotheta;o.dorho=dorho;o.doPI=doPI;doPI=T;dorho=dotheta=doMu=F;
    if (verbose) cat("Inferring ancestry switching rates holding other parameters fixed\n");
    total=PI.total
    tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NUMP, kLL, L,
	       NUMI, max.donors, G, gobs, maxmatchsize, umatch, flips, maxmiss, d.w, t.w,  total, verbose=F, len, cloglike, LOG, EMlogfile, doPI, doMu, 
	       dotheta, dorho)
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
    cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
    if (!absorbrho | !commonrho | !commontheta) 
    {
      tmp=all_donates(NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=F, max.donors, NUMP, G, umatch, maxmatchsize, d.w, 
		      t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len, LOG, transitions, mutmat)
      ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
    }
    doMu=o.doMu;dorho=o.dorho;dotheta=o.dotheta;doPI=o.doPI
  }

  total<-s.total # number of EM steps with repeated loop
  if (EM) 
  {
    for (reps in 1:REPS) 
    {
      cat("######################## round ", reps, "of ",  REPS, "#######################\n")
      if (reps==REPS) 
	M=o.M # on last rep do more MCMC phasing
      # location of this an issue. If above thin&phase, low lambda. If below then first rep has lowered log-like
      tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NUMP, kLL, L,
		 NUMI, max.donors, G, gobs, maxmatchsize, umatch, flips, maxmiss, d.w, t.w,  total, verbose=F, len, cloglike, LOG, EMlogfile, doPI, doMu, 
		 dotheta, dorho) 
      PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
      cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
      old.kLL=kLL
      if (max.donors<NUMP & (kLL==old.kLL))
      {
	# decide on donor set using updated parameters
	tmp=all_donates(NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=F, max.donors, NUMP, G, umatch, maxmatchsize, d.w, 
			t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len, LOG, transitions, mutmat)
	ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
      }
      if (PHASE) 
      {
	flips=phase_hunt_all(donates, donatesl, donatesr, ndonors, nchrno, NUMI, flips, eps.lower, 
			     transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, 
			     initProb, PLOT, min.bg, max.bg, HPC, verbose, LOG) 
	if (M>0)  # now do MCMC re-phasing w/o output to console (ugly due to txtProgressBar)
	  flips=phase_mcmc_all(donates, donatesl, donatesr, ndonors, nchrno, NUMI, flips, 
			       transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, 
			       initProb, PLOT, HPC, verbose, LOG) 
      }
    }
    total=o.total # do longer run EM on last rep
    if (verbose)
      cat("run one final round of EM\n")
    tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NUMP, kLL, L,
	       NUMI, max.donors, G, gobs, maxmatchsize, umatch, flips, maxmiss, d.w, t.w,  total, verbose=F, len, cloglike, LOG, EMlogfile, 
	       doPI, doMu, dotheta, dorho) 
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
    cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
  }
  ans$final.flips=flips
  EM=F;getnoancgfbs=T;eps=log(1.01);LOG=F;PLOT=F;
  a.Mu=Mu;a.rho=rho;a.theta=theta;a.PI=PI;a.alpha=alpha;a.lambda=lambda
  a.o.Mu=o.Mu;a.o.rho=o.rho;a.o.theta=o.theta;a.o.PI=o.PI;a.o.alpha=o.alpha;a.o.lambda=o.lambda
  samp_chrnos=chrnos;

  ######### fully Mosaic curves with Mosaic phasing ############
  ans$gfbs=get_gfbs(NUMP, max.donors, donates, donatesl, donatesr, NUMA, L, G, kLL, transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, initProb, 
		label, ndonors, flips)
  if (verbose) cat("saving localanc results to file\n")
  if (target!="simulated") tmp=get_localanc(ans$gfbs,G,L,kLL,NUMA,NUMI)
  if (target=="simulated") tmp=get_localanc(ans$gfbs,G,L,kLL,NUMA,NUMI,t.g.true_anc=g.true_anc)
  ans$localanc=tmp$localanc
  if (target=="simulated") 
    ans$g.true_anc=tmp$g.true_anc
  save(file=paste0(resultsdir,"localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), ans$localanc, ans$final.flips, ans$g.loc)
  if (target=="simulated")
    save(file=paste0(resultsdir,"localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		     "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), ans$localanc, ans$g.true_anc, ans$final.flips, ans$g.loc)
  save(file=paste0(resultsdir,"gfbs_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), ans$gfbs)

  if (verbose) cat("calculating ancestry aware re-phased coancestry curves\n"); acoancs=create_coancs(ans$localanc,dr,"DIP");
  ######## GlobeTrotter style curves original phasing ##########
  for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # undo phase flips
  tmp=fit_noanc_model(samp_chrnos, chrnos, NUMA, NUMP, kLL, L, KNOWN, label, NL, NN, umatch, G, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
		      prop.don, max.donors, maxmatch, maxmiss, initProb, d.w, t.w, NUMA, max(NL), HPC, getnoancgfbs=TRUE) 
  transitions=tmp$t.transitions;mutmat=tmp$mutmat;Mu=tmp$Mu;theta=tmp$theta;rho=tmp$rho
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
  noanc_gfbs=tmp$noanc_gfbs
  if (HPC) cleanup_ff_files(donates, donatesl, donatesr, nchrno, NUMI, ffpath, FALSE)
  Mu=a.Mu;rho=a.rho;theta=a.theta;PI=a.PI;alpha=a.alpha;lambda=a.lambda
  o.Mu=a.o.Mu;o.rho=a.o.rho;o.theta=a.o.theta;o.PI=a.o.PI;o.alpha=a.o.alpha;o.lambda=a.o.lambda
  ans$noanc_unphased_localanc=get_ancunaware_localanc(NUMA,L,G,nchrno,noanc_gfbs,Mu,alpha) # works off noanc_gfbs
  save(file=paste0(resultsdir,"noanc_unphased_localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), ans$noanc_unphased_localanc, flips, ans$g.loc)
  if (verbose) cat("calculating ancestry unaware input phasing coancestry curves\n"); coancs=create_coancs(ans$noanc_unphased_localanc,dr,"DIP")

  if (verbose) cat("saving final results to file\n")
  save(file=paste0(resultsdir,"",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",
		   GpcM,"_",prop.don,"_",max.donors,".RData"), target, o.Mu, o.lambda, o.theta, o.alpha, o.PI, o.rho, 
       Mu, lambda, theta, alpha, PI, rho, L, NUMA, nchrno, chrnos, dr, NL, kLL, acoancs, coancs)
  ans$target=target;ans$o.Mu=o.Mu;ans$o.lambda=o.lambda;ans$o.theta=o.theta;ans$o.alpha=o.alpha;ans$o.PI=o.PI;ans$o.rho=o.rho
  ans$Mu=Mu;ans$lambda=lambda;ans$theta=theta;ans$alpha=alpha;ans$PI=PI;ans$rho=rho
  ans$L=L;ans$NUMA=NUMA;ans$nchrno=nchrno;ans$chrnos=chrnos;ans$dr=dr;ans$NL=NL;ans$kLL=kLL;ans$acoancs=acoancs;ans$coancs=coancs
  return(ans)
}
