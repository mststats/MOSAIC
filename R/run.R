run_mosaic=function(target,datasource,chrnos,A,NUMI,pops=NULL,mask=NULL,PLOT=TRUE,doFst=TRUE,PHASE=TRUE,gens=NULL,ratios=NULL,EM=TRUE, 
		    ffpath=tempdir(),MC=0,return.res=TRUE,REPS=0,GpcM=60,nl=1000,max.donors=100,prop.don=0.99,
		    doMu=TRUE,doPI=TRUE,dorho=TRUE,dotheta=TRUE,firstind=1,verbose=TRUE,Ne=9e4,MODE="DIP",singlePI=FALSE, 
		    init.rho=NULL, init.theta=NULL, init.Mu=NULL, init.PI=NULL, commonrho=TRUE, commontheta=TRUE, resultsdir="MOSAIC_RESULTS") {
  if (!dir.exists(ffpath))
    stop(paste("Requested location ", ffpath, " for fast-file storage not found\n"))
  nchrno=length(chrnos) # number of chromosomes for these target haplotypes
  # sets default parameters, sets up some objects required later, reads in data, and initialises model.
  if (A<2)
    stop("need to fit at least a 2-way model\n")
  if (target=="simulated" & length(pops)<A)
    stop("Please provide ", A, " groups to simulate from\n")
  if (MODE=="HAP") {
    warning("########## Haploid mode is under construction; for now ancestry recombination rates are equal for consecutive haplotypes as per diploid runs  ##########", immediate.=TRUE)
  }
  if (MODE=="HAP" & PHASE)
    stop("are you trying to perform phasing on haploid data?")
  if (!is.null(init.Mu)) {
    if (ncol(init.Mu)!=A)
      stop(paste("There should be A columns of Mu. You have provided", ncol(init.Mu), "columns for", A, "ancestries"))
    if (!all.equal(colSums(init.Mu),rep(1,A), tolerance=0.01)) 
      warning("########## column sums of Mu must all be one: rescaling supplied values ##########", immediate.=TRUE)
    init.Mu=t(t(init.Mu)/colSums(init.Mu))
  }
  if (!is.null(init.PI)) {
    if (length(init.PI)!=NUMI)
      stop(paste("There should be NUMI PI matrices. You have provided", length(init.PI), "matrices for", NUMI, "individuals"))
    if (any(sapply(init.PI, dim)!=A))
      stop(paste("All PI matrices should be AxA but you have A =", A, "ancestries and PI dimensions of", toString(sapply(init.PI, function(x) paste0(nrow(x),"x",ncol(x))))))
  }
  setup=setup_data_etc(NUMI,target,chrnos,pops,A,datasource,EM,gens,ratios,MC,prop.don=prop.don,max.donors=max.donors,
		       firstind=firstind,REPS=REPS,GpcM=GpcM,nl=nl,mask=mask,PHASE=PHASE,Ne=Ne,singlePI=singlePI, 
		       init.rho=init.rho, init.theta=init.theta, init.PI=init.PI, commonrho=commonrho, commontheta=commontheta, resultsdir=resultsdir) 
  resultsdir=setup$resultsdir;PHASE=setup$PHASE;HPC=setup$HPC;GpcM=setup$GpcM;LOG=setup$LOG
  mcmcprog=setup$mcmcprog;absorbrho=setup$absorbrho;commonrho=setup$commonrho;commontheta=setup$commontheta;prethin=setup$prethin
  s.M=setup$s.M;M=setup$M;PI.total=setup$PI.total;s.total=setup$s.total;REPS=setup$REPS
  eps.lower=setup$eps.lower;min.bg=setup$min.bg;max.bg=setup$max.bg;samp_chrnos=setup$samp_chrnos;dr=setup$dr
  FLAT=setup$FLAT;maxmatch=setup$maxmatch;maxmiss=setup$maxmiss;umatch=setup$umatch;d.w=setup$d.w
  t.w=setup$t.w;g.loc=setup$g.loc;gobs=setup$gobs;NUMP=setup$NUMP;NUMI=setup$NUMI;NUMA=setup$NUMA
  label=setup$label;KNOWN=setup$KNOWN;kLL=setup$kLL;NL=setup$NL;G=setup$G;
  NN=setup$NN;maxmatchsize=setup$maxmatchsize;panels=setup$panels;min.donors=setup$min.donors;
  theta=setup$theta;rho=setup$rho;lambda=setup$lambda;alpha=setup$alpha;PI=setup$PI;Mu=setup$MU
  transitions=setup$transitions;total=setup$total
  o.PI=setup$o.PI;o.alpha=setup$o.alpha;o.rho=setup$o.rho;o.lambda=setup$o.lambda;o.theta=setup$o.theta;o.phi.theta=setup$o.phi.theta 
  flips=setup$flips;prop.don=setup$prop.don;max.donors=setup$max.donors
  if (!is.null(init.Mu)) {
    if (nrow(init.Mu)!=kLL)
      stop(paste("There should be same number of row of Mu as panels. You have provided", nrow(init.Mu), "rows for", kLL, "panels"))
  }
  if (target=="simulated")
    g.true_anc=setup$g.true_anc
  rm(setup)
  old.runtime<-as.numeric(Sys.time())
  o.total=total
  eps=log(1.01) # i.e. a 1% increase in relative likelihood
  Mu<-matrix(rep(1/kLL,A*kLL),kLL)
  #for (ind in 1:NUMI) alpha[[ind]]=rep(1/A,A) # flatten out w.r.t. ancestry 
  runtime=NaN
  if (any(is.null(init.Mu),is.null(init.theta),is.null(init.rho),is.null(init.PI))) {
    if (verbose) cat("Initialise parameters of MOSAIC based on ancestry unaware copying probabilities\n");
    # run noanc.R b/c need good paras for init_Mu
    tmp=fit_noanc_model(target, samp_chrnos, chrnos, NUMA, NUMP, kLL, A, KNOWN, label, NL, NN, umatch, G, dr, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
			prop.don, min.donors, max.donors, maxmatchsize, maxmatch, maxmiss, initProb, d.w, t.w, NUMA, 100, HPC, runtime, resultsdir, GpcM, eps, NaN,
			doMu, doPI, dorho, dotheta, ffpath, firstind, EM) 
    transitions=tmp$t.transitions;mutmat=tmp$mutmat;Mu=tmp$Mu
    if (is.null(init.theta)) theta=tmp$theta
    if (is.null(init.rho)) rho=tmp$rho
    ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
  }
  initProb=initprobs(TRUE,NUMA,A,NUMP,kLL,PI,Mu,rho,alpha,label,NL)
  runtime<-as.numeric(Sys.time())
  if (is.null(init.Mu)) {
    if (kLL>A) # otherwise can't cluster kLL things into A clusters
    {
      # use this to get #switches in noanc model w/o writelog
      tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=TRUE, t.get_switches=TRUE, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, maxmatchsize,
		      maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=FALSE, NUMA, nchrno, initProb, runtime, NULL,FALSE,transitions,mutmat,NaN,NULL,ffpath)
      ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike
      noanc_gswitches=tmp$noanc_gswitches
      rm(tmp)
      windowed_copying<-window_chunks(nswitches=noanc_gswitches,dr,G,kLL,NUMA,ww=0.5,verbose=verbose) # similar in the same windows
      rm(noanc_gswitches) 
      tmp<-cluster_windows(windowed_copying,dr,kLL,A,NUMI,NUMA,NL,absorbrho,verbose=verbose)
      Mu<-tmp$Mu
      if (is.null(init.PI)) {
	PI<-tmp$PI
	tmp=alphalambda_from_PI(PI,dr)
	if (is.null(ratios)) alpha<-tmp$alpha
	if (is.null(gens)) {
	  lambda=tmp$lambda # or could start with deliberately lower values as above method overestimates lambda
	  #lambda=o.lambda # re-use o.lambda from above
	}
      }
      rm(windowed_copying,tmp)
    } else {diag(Mu)=10*diag(Mu);Mu=t(t(Mu)/colSums(Mu))}
  } else Mu=init.Mu
  rownames(Mu)<-panels[1:kLL]
  o.Mu<-Mu;o.alpha<-alpha;o.lambda=lambda;o.PI=PI # these are the official starting ancestry related parameters now
  mutmat<-fmutmat(theta, A, maxmiss, maxmatch); for (ind in 1:NUMI) transitions[[ind]]<-s_trans(A,kLL,PI[[ind]],Mu,rho,NL)
  o.M<-M;M<-s.M
  logfile=NULL
  tmp=create_logfile(resultsdir,target,kLL,A,NUMI,firstind,chrnos,nchrno,NN,GpcM)
  runtime=old.runtime=tmp$rtime;diff.time=0;len=tmp$len
  logfile=tmp$logfile
  # decide on donor set using initial parameters
  tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=TRUE, t.get_switches=FALSE, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, maxmatchsize, 
		  maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=FALSE, NUMA, nchrno, initProb, runtime, len, FALSE, transitions, mutmat,NaN,logfile,ffpath)
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike

  ############ a few PI only updates first; very useful to do before first re-phasing
  if (PI.total>0 & EM & doPI)
  {
    o.doMu=doMu;o.dotheta=dotheta;o.dorho=dorho;o.doPI=doPI;doPI=TRUE;dorho=dotheta=doMu=FALSE;
    if (verbose) cat("Inferring ancestry switching rates holding other parameters fixed\n");
    total=PI.total
    tmp=run_EM(HPC, nchrno, PI, Mu, rho, theta, alpha, lambda, initProb, label, mutmat, transitions, ndonors, donates, donatesl, donatesr, NUMA, NN, NL, NUMP, kLL, A,
	       NUMI, max.donors, G, dr, gobs, maxmatchsize, umatch, flips, maxmatch, maxmiss, d.w, t.w,  total, verbose=FALSE, len, cloglike, LOG, logfile, doPI, doMu, 
	       dotheta, dorho, commonrho, commontheta, absorbrho, singlePI, runtime, eps)
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
    cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
    if (!absorbrho | !commonrho | !commontheta) 
    {
      tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=TRUE, t.get_switches=FALSE, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, 
		      maxmatchsize, maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=FALSE, NUMA, nchrno, initProb, runtime, len, LOG, transitions, 
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
		 NUMI, max.donors, G, dr, gobs, maxmatchsize, umatch, flips, maxmatch, maxmiss, d.w, t.w,  total, verbose=FALSE, len, cloglike, LOG, logfile, doPI, doMu, 
		 dotheta, dorho, commonrho, commontheta, absorbrho, singlePI, runtime, eps) 
      PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
      cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
    }
    old.kLL=kLL
    if (max.donors<NUMP & (kLL==old.kLL))
    {
      # decide on donor set using updated parameters
      tmp=all_donates(target, A, NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=TRUE, t.get_switches=FALSE, min.donors, max.donors, prop.don, NUMP, NL, G, umatch, 
		      maxmatchsize, maxmatch, maxmiss, d.w, t.w, gobs, flips, label, KNOWN, HPC, prethin=FALSE, NUMA, nchrno, initProb, runtime, len, LOG, transitions, 
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
	       NUMI, max.donors, G, dr, gobs, maxmatchsize, umatch, flips, maxmatch, maxmiss, d.w, t.w,  total, verbose=FALSE, len, cloglike, LOG, logfile, 
	       doPI, doMu, dotheta, dorho, commonrho, commontheta, absorbrho, singlePI, runtime, eps) 
    PI=tmp$PI;alpha=tmp$alpha;lambda=tmp$lambda;Mu=tmp$Mu;rho=tmp$rho;theta=tmp$theta;runtime=tmp$runtime;initProb=tmp$initProb;
    cloglike=tmp$cloglike;transitions=tmp$transitions;mutmat=tmp$mutmat
  }
  final.flips=flips

  ######### MOSAIC curves with MOSAIC phasing ############
  if (verbose) cat("saving localanc results to file\n")
  if (target!="simulated") tmp=get_localanc(NUMP, nchrno, max.donors, donates, donatesl, donatesr, transitions, umatch, maxmatchsize, d.w, t.w, 
					    gobs, mutmat, maxmiss, initProb, label, ndonors, flips, HPC, G, A, kLL, NUMA, NUMI)
  if (target=="simulated") tmp=get_localanc(NUMP, nchrno, max.donors, donates, donatesl, donatesr, transitions, umatch, maxmatchsize, d.w, t.w, 
					    gobs, mutmat, maxmiss, initProb, label, ndonors, flips, HPC, G, A, kLL, NUMA, NUMI, t.g.true_anc=g.true_anc)
  localanc=tmp$localanc
  if (target=="simulated") 
    g.true_anc=tmp$g.true_anc

  if (target!="simulated")
    save(file=file.path(resultsdir,paste0("localanc_",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		     "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData")), localanc, final.flips, g.loc)
  if (target=="simulated")
    save(file=file.path(resultsdir,paste0("localanc_",target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		     "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData")), localanc, g.true_anc, final.flips, g.loc)
  if (verbose) cat("calculating ancestry aware re-phased coancestry curves\n"); acoancs=create_coancs(localanc,dr,MODE);
  all_Fst=NULL
  if (doFst) {
    if (verbose) cat("calculating Fst values\n")
    flocalanc=phase_localanc(localanc,final.flips) 
    if (target=="simulated")
      write_admixed_summary(target,NL,pathout=file.path(resultsdir,"FREQS"),targetdatasource=resultsdir,datasource=datasource,g.loc=g.loc,t.localanc=flocalanc,chrnos=chrnos)
    if (target!="simulated")
      write_admixed_summary(target,NL,pathout=file.path(resultsdir,"FREQS"),targetdatasource=datasource,datasource=datasource,g.loc=g.loc,t.localanc=flocalanc,chrnos=chrnos)
    write_panel_summaries(pathout=file.path(resultsdir,"FREQS"), panels=rownames(Mu), datasource=datasource,chrnos=chrnos)
    all_Fst=Fst_combos(target, A, sum(NL), rownames(Mu), pathin=file.path(resultsdir, "FREQS"))
  }


  if (verbose) cat("saving final results to file\n")
  save(file=file.path(resultsdir,paste0(target,"_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",
		   GpcM,"_",prop.don,"_",max.donors,".RData")), target, logfile, #o.Mu, o.lambda, o.theta, o.alpha, o.PI, o.rho, 
       Mu, lambda, theta, alpha, PI, rho, A, NUMA, nchrno, chrnos, dr, NL, kLL, acoancs, all_Fst, GpcM)

  cat("Expected r-squared (genomewide):", dip_expected_fr2(localanc),"\n")
  if (target=="simulated") 
    cat("Actual r-squared (genomewide):", dip_fr2(localanc,g.true_anc),"\n")
  if (!is.null(all_Fst)) {
    if (all(!is.nan(all_Fst$ancs))) {
      cat("Fst between mixing groups:\n")
      print(all_Fst$ancs)
    } 
    if (all(!is.nan(all_Fst$Rst))) {
      cat("Rst between mixing groups:\n")
      print(all_Fst$Rst)
    } 
    if (any(is.nan(all_Fst$ancs)) | any(is.nan(all_Fst$Rst))) {
      warning("########## cannot estimate Fst: insufficient loci mapped to different ancestries across target individuals ##########", immediate.=TRUE)
      cat("average alpha=", Reduce("+",alpha)/length(alpha),"\n")
    }
  }
  if (PLOT) {
    if (verbose) cat("saving plots to ", resultsdir, "folder\n")
    if (target!="simulated")
      plot_all_mosaic(resultsdir,target,EM,PHASE,t.GpcM=GpcM,t.all_Fst=all_Fst,t.A=A,t.NUMA=NUMA,
		      t.Mu=Mu,t.chrnos=chrnos,t.alpha=alpha,t.NL=NL,t.acoancs=acoancs,t.dr=dr,t.logfile=logfile,localanc, g.loc)
    if (target=="simulated")
      plot_all_mosaic(resultsdir,target,EM,PHASE,t.GpcM=GpcM,t.all_Fst=all_Fst,t.A=A,t.NUMA=NUMA,
		      t.Mu=Mu,t.chrnos=chrnos,t.alpha=alpha,t.NL=NL,t.acoancs=acoancs,t.dr=dr,t.logfile=logfile,localanc,g.loc,g.true_anc)
  }

  if (return.res & target!="simulated")
    return(list(g.loc=g.loc,localanc=localanc,logfile=logfile,final.flips=final.flips,dr=dr,A=A,kLL=kLL,NUMP=NUMP,NN=NN,NL=NL,label=label,chrnos=chrnos,
		NUMA=NUMA,NUMI=NUMI,GpcM=GpcM,PI=PI,lambda=lambda,alpha=alpha,Mu=Mu,theta=theta,rho=rho,acoancs=acoancs,target=target,
		prop.don=prop.don,max.donors=max.donors,all_Fst=all_Fst))
  if (return.res & target=="simulated")
    return(list(g.loc=g.loc,localanc=localanc,logfile=logfile,g.true_anc=g.true_anc,final.flips=final.flips,dr=dr,A=A,kLL=kLL,NUMP=NUMP,NN=NN,NL=NL,label=label,chrnos=chrnos,
		NUMA=NUMA,NUMI=NUMI,GpcM=GpcM,PI=PI,lambda=lambda,alpha=alpha,Mu=Mu,theta=theta,rho=rho,acoancs=acoancs,target=target,
		prop.don=prop.don,max.donors=max.donors,all_Fst=all_Fst))
}
