run_mosaic=function(target,datasource,chrnos,A,NUMI,pops=NULL,mask=NULL,PLOT=FALSE,doFst=TRUE,PHASE=TRUE,gens=0,ratios=NULL,EM=TRUE,
			 ffpath="/dev/shm/",MC=0,return.res=TRUE,REPS=0,GpcM=60,nl=1000,max.donors=100,prop.don=0.99,
			 doMu=TRUE,doPI=TRUE,dorho=TRUE,dotheta=TRUE,firstind=1,verbose=TRUE,Ne=9e4) {
  nchrno=length(chrnos) # number of chromosomes for these target haplotypes
  # sets default parameters, sets up some objects required later, reads in data, and initialises model.
  if (A<2)
    stop("need to fit at least a 2-way model\n")
  if (target=="simulated" & length(pops)<A)
    stop("Please provide ", A, " groups to simulate from\n")
  tmp=setup_data_etc(NUMI,target,chrnos,pops,A,datasource,EM,gens,ratios,MC,REPS=REPS,GpcM=GpcM,nl=nl,mask=mask,PHASE=PHASE,Ne=Ne) 
  resultsdir=tmp$resultsdir;PHASE=tmp$PHASE;HPC=tmp$HPC;GpcM=tmp$GpcM;LOG=tmp$LOG
  mcmcprog=tmp$mcmcprog;absorbrho=tmp$absorbrho;commonrho=tmp$commonrho;commontheta=tmp$commontheta;prethin=tmp$prethin
  s.M=tmp$s.M;M=tmp$M;PI.total=tmp$PI.total;s.total=tmp$s.total;REPS=tmp$REPS
  eps.lower=tmp$eps.lower;min.bg=tmp$min.bg;max.bg=tmp$max.bg;samp_chrnos=tmp$samp_chrnos;dr=tmp$dr
  FLAT=tmp$FLAT;maxmatch=tmp$maxmatch;maxmiss=tmp$maxmiss;umatch=tmp$umatch;d.w=tmp$d.w
  t.w=tmp$t.w;g.loc=tmp$g.loc;gobs=tmp$gobs;NUMP=tmp$NUMP;NUMI=tmp$NUMI;NUMA=tmp$NUMA
  label=tmp$label;KNOWN=tmp$KNOWN;kLL=tmp$kLL;NL=tmp$NL;G=tmp$G;
  NN=tmp$NN;maxmatchsize=tmp$maxmatchsize;panels=tmp$panels;min.donors=tmp$min.donors;
  theta=tmp$theta;rho=tmp$rho;lambda=tmp$lambda;alpha=tmp$alpha;PI=tmp$PI;Mu=tmp$MU
  transitions=tmp$transitions;total=tmp$total
  o.PI=tmp$o.PI;o.alpha=tmp$o.alpha;o.rho=tmp$o.rho;o.lambda=tmp$o.lambda;o.theta=tmp$o.theta;o.phi.theta=tmp$o.phi.theta 
  flips=tmp$flips;prop.don=tmp$prop.don;max.donors=tmp$max.donors
  if (target=="simulated")
    g.true_anc=tmp$g.true_anc
  rm(tmp)
  old.runtime<-as.numeric(Sys.time())
  o.total=total
  eps=log(1.01) # i.e. a 1% increase in relative likelihood
  Mu<-matrix(rep(1/kLL,A*kLL),kLL);for (ind in 1:NUMI) alpha[[ind]]=rep(1/A,A) # flatten out w.r.t. ancestry
  runtime=NaN
  # always need to run noanc.R b/c need good paras for init_Mu
  tmp=fit_noanc_model(target, samp_chrnos, chrnos, NUMA, NUMP, kLL, A, KNOWN, label, NL, NN, umatch, G, dr, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
		      prop.don, min.donors, max.donors, maxmatchsize, maxmatch, maxmiss, initProb, d.w, t.w, NUMA, 100, HPC, runtime, resultsdir, GpcM, eps, NaN,
		      doMu, doPI, dorho, dotheta, ffpath, firstind, EM) 
  transitions=tmp$t.transitions;mutmat=tmp$mutmat;Mu=tmp$Mu;theta=tmp$theta;rho=tmp$rho
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
  initProb=initprobs(T,NUMA,A,NUMP,kLL,PI,Mu,rho,alpha,label,NL)
  runtime<-as.numeric(Sys.time())
  if (kLL>A) # otherwise can't cluster kLL things into A clusters
  {
    # use this to get #switches in noanc model w/o writelog
    tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=T, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, maxmatchsize,
		    maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, NULL,F,transitions,mutmat,NaN,NULL,ffpath)
    ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
    noanc_gswitches=tmp$noanc_gswitches
    rm(tmp)
    windowed_copying<-window_chunks(nswitches=noanc_gswitches,dr,G,kLL,NUMA,ww=0.5,verbose=verbose) # similar in the same windows
    rm(noanc_gswitches) 
    tmp<-cluster_windows(windowed_copying,dr,kLL,A,NUMI,NUMA,NL,absorbrho,verbose=verbose)
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
  mutmat<-fmutmat(theta, A, maxmiss, maxmatch); for (ind in 1:NUMI) transitions[[ind]]<-s_trans(A,kLL,PI[[ind]],Mu,rho,NL)
  o.M<-M;M<-s.M
  logfile=NULL
  tmp=create_logfile(resultsdir,target,kLL,A,NUMI,firstind,chrnos,nchrno,NN,GpcM)
  runtime=old.runtime=tmp$rtime;diff.time=0;len=tmp$len
  logfile=tmp$logfile
  # decide on donor set using initial parameters
  tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=F, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, maxmatchsize, 
		  maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len, F, transitions, mutmat,NaN,logfile,ffpath)
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike

  ############ a few PI only updates first; very useful to do before first re-phasing
  if (PI.total>0 & EM & doPI)
  {
    o.doMu=doMu;o.dotheta=dotheta;o.dorho=dorho;o.doPI=doPI;doPI=T;dorho=dotheta=doMu=F;
    if (verbose) cat("Inferring ancestry switching rates holding other parameters fixed\n");
    total=PI.total
    tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, label, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NN, NL, NUMP, kLL, A,
	       NUMI, max.donors, G, dr, gobs, maxmatchsize, umatch, flips, maxmatch, maxmiss, d.w, t.w,  total, verbose=F, len, cloglike, LOG, logfile, doPI, doMu, 
	       dotheta, dorho, commonrho, commontheta, absorbrho, runtime, eps)
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
    cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
    if (!absorbrho | !commonrho | !commontheta) 
    {
      tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=F, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, 
		      maxmatchsize, maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len, LOG, transitions, 
		      mutmat,cloglike,logfile,ffpath)
      ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
    }
    doMu=o.doMu;dorho=o.dorho;dotheta=o.dotheta;doPI=o.doPI
  }

  total<-s.total # number of EM steps with repeated loop
    for (reps in 1:REPS) 
    {
      cat("######################## round ", reps, "of ",  REPS, "#######################\n")
      if (reps==REPS) 
	M=o.M # on last rep do more MCMC phasing
  if (EM) {
      tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, label, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NN, NL, NUMP, kLL, A,
		 NUMI, max.donors, G, dr, gobs, maxmatchsize, umatch, flips, maxmatch, maxmiss, d.w, t.w,  total, verbose=F, len, cloglike, LOG, logfile, doPI, doMu, 
		 dotheta, dorho, commonrho, commontheta, absorbrho, runtime, eps) 
      PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
      cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
  }
      old.kLL=kLL
      if (max.donors<NUMP & (kLL==old.kLL))
      {
	# decide on donor set using updated parameters
	tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=F, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, 
			maxmatchsize, maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len, LOG, transitions, 
			mutmat,cloglike,logfile,ffpath)
	ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
      }
      if (PHASE) 
      {
	tmp=phase_hunt_all(donates, donatesl, donatesr, ndonors, NUMP, G, A, GpcM, max.donors, nchrno, NUMA, NUMI, kLL, flips, eps.lower, 
			     transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, 
			     initProb, label, min.bg, max.bg, len,Mu,rho,PI,alpha,lambda,theta,runtime, logfile, HPC, verbose, LOG) 
	flips=tmp$flips
	runtime=tmp$runtime
	cloglike=tmp$cloglike
	rm(tmp)
	if (M>0)  # now do MCMC re-phasing w/o output to console (ugly due to txtProgressBar)
	{
	  tmp=phase_mcmc_all(donates, donatesl, donatesr, ndonors, NUMP, NUMA, G, A, max.donors, nchrno, NUMI, kLL, flips, 
			       transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, 
			       initProb, label, len,Mu,rho,PI,alpha,lambda,theta,runtime, logfile, HPC, verbose, LOG) 
	flips=tmp$flips
	runtime=tmp$runtime
	cloglike=tmp$cloglike
	rm(tmp)
	}
      }
    }
  if (EM) {
    total=o.total # do longer run EM on last rep
    if (verbose)
      cat("run one final round of EM\n")
    tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, label, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NN, NL, NUMP, kLL, A,
	       NUMI, max.donors, G, dr, gobs, maxmatchsize, umatch, flips, maxmatch, maxmiss, d.w, t.w,  total, verbose=F, len, cloglike, LOG, logfile, 
	       doPI, doMu, dotheta, dorho, commonrho, commontheta, absorbrho, runtime, eps) 
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
    cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
  }
  final.flips=flips
  EM=F;getnoancgfbs=T;eps=log(1.01);LOG=F;
  a.Mu=Mu;a.rho=rho;a.theta=theta;a.PI=PI;a.alpha=alpha;a.lambda=lambda
  a.o.Mu=o.Mu;a.o.rho=o.rho;a.o.theta=o.theta;a.o.PI=o.PI;a.o.alpha=o.alpha;a.o.lambda=o.lambda
  samp_chrnos=chrnos;

  ######### MOSAIC curves with MOSAIC phasing ############
  gfbs=get_gfbs(NUMP, nchrno, max.donors, donates, donatesl, donatesr, NUMA, A, G, kLL, transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, initProb, 
		label, ndonors, flips, HPC)
  if (verbose) cat("saving localanc results to file\n")
  if (target!="simulated") tmp=get_localanc(gfbs,G,A,kLL,NUMA,NUMI)
  if (target=="simulated") tmp=get_localanc(gfbs,G,A,kLL,NUMA,NUMI,t.g.true_anc=g.true_anc)
  localanc=tmp$localanc
  if (target=="simulated") 
    g.true_anc=tmp$g.true_anc

  save(file=paste0(resultsdir,"localanc_",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), localanc, final.flips, g.loc)
  if (target=="simulated")
    save(file=paste0(resultsdir,"localanc_",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		     "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), localanc, g.true_anc, final.flips, g.loc)
  save(file=paste0(resultsdir,"gfbs_",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), gfbs)

  if (verbose) cat("calculating ancestry aware re-phased coancestry curves\n"); acoancs=create_coancs(localanc,dr,"DIP");
  all_Fst=NULL
  if (doFst) {
    if (verbose) cat("calculating Fst values\n")
    flocalanc=phase_localanc(localanc,final.flips) 
    if (target=="simulated")
      write_admixed_summary(target,NL,targetdatasource=resultsdir,datasource=datasource,g.loc=g.loc,t.localanc=flocalanc,chrnos=chrnos)
    if (target!="simulated")
      write_admixed_summary(target,NL,targetdatasource=datasource,datasource=datasource,g.loc=g.loc,t.localanc=flocalanc,chrnos=chrnos)
    write_panel_summaries(panels=rownames(Mu),datasource=datasource,,chrnos=chrnos)
    all_Fst=Fst_combos(target, A, sum(NL), rownames(Mu))
  }

  ######## GlobeTrotter style curves original phasing ##########
  for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # undo phase flips
  tmp=fit_noanc_model(target, samp_chrnos, chrnos, NUMA, NUMP, kLL, A, KNOWN, label, NL, NN, umatch, G, dr, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
		      prop.don, min.donors, max.donors, maxmatchsize, maxmatch, maxmiss, initProb, d.w, t.w, NUMA, max(NL), HPC, runtime, resultsdir, 
		      GpcM, eps, NaN, doMu, doPI, dorho, dotheta, ffpath, firstind, EM, getnoancgfbs=TRUE) 
  transitions=tmp$t.transitions;mutmat=tmp$mutmat;Mu=tmp$Mu;theta=tmp$theta;rho=tmp$rho
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
  noanc_gfbs=tmp$noanc_gfbs
  if (HPC) cleanup_ff_files(donates, donatesl, donatesr, nchrno, NUMI, ffpath, FALSE)
  Mu=a.Mu;rho=a.rho;theta=a.theta;PI=a.PI;alpha=a.alpha;lambda=a.lambda
  o.Mu=a.o.Mu;o.rho=a.o.rho;o.theta=a.o.theta;o.PI=a.o.PI;o.alpha=a.o.alpha;o.lambda=a.o.lambda
  noanc_unphased_localanc=get_ancunaware_localanc(NUMA,A,G,nchrno,noanc_gfbs,Mu,alpha) # works off noanc_gfbs
  save(file=paste0(resultsdir,"noanc_unphased_localanc_",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), noanc_unphased_localanc, flips, g.loc)
  if (verbose) cat("calculating ancestry unaware input phasing coancestry curves\n"); coancs=create_coancs(noanc_unphased_localanc,dr,"DIP")

  if (verbose) cat("saving final results to file\n")
  save(file=paste0(resultsdir,"",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",
		   GpcM,"_",prop.don,"_",max.donors,".RData"), target, logfile, #o.Mu, o.lambda, o.theta, o.alpha, o.PI, o.rho, 
       Mu, lambda, theta, alpha, PI, rho, A, NUMA, nchrno, chrnos, dr, NL, kLL, acoancs, coancs, all_Fst)

  cat("Expected r-squared (genomewide):", dip_expected_fr2(localanc),"\n")
  if (target=="simulated") 
    cat("Actual r-squared (genomewide):", dip_fr2(localanc,g.true_anc),"\n")
  cat("Fst between mixing groups:\n")
  print(all_Fst$ancs)
  cat("Rst between mixing groups:\n")
  print(all_Fst$Rst)


  if (return.res & target!="simulated")
    return(list(g.loc=g.loc,localanc=localanc,logfile=logfile,final.flips=final.flips,dr=dr,A=A,kLL=kLL,NUMP=NUMP,NN=NN,NL=NL,label=label,chrnos=chrnos,
		NUMA=NUMA,NUMI=NUMI,GpcM=GpcM,PI=PI,lambda=lambda,alpha=alpha,Mu=Mu,theta=theta,rho=rho,acoancs=acoancs,coancs=coancs,target=target,
		prop.don=prop.don,max.donors=max.donors,all_Fst=all_Fst))
  if (return.res & target=="simulated")
    return(list(g.loc=g.loc,localanc=localanc,logfile=logfile,g.true_anc=g.true_anc,final.flips=final.flips,dr=dr,A=A,kLL=kLL,NUMP=NUMP,NN=NN,NL=NL,label=label,chrnos=chrnos,
		NUMA=NUMA,NUMI=NUMI,GpcM=GpcM,PI=PI,lambda=lambda,alpha=alpha,Mu=Mu,theta=theta,rho=rho,acoancs=acoancs,coancs=coancs,target=target,
		prop.don=prop.don,max.donors=max.donors,all_Fst=all_Fst))
}
