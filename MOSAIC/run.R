# script to run all code required to fit MOSAIC. Reads in data, initialises, performs thin->phase->EM cycles and outputs results.
# example usage: Rscript run.R SanKhomani HGDP/ 4 1 60 16
require(mosaicpackage)
######################## first set some options ###############################
# important things to set
shargs<-commandArgs(trailingOnly=TRUE) # read in arguments from the command line; 
target=shargs[1]  # e.g. SpainPopn_1
datasource=shargs[2]; # i.e. spanish/ ; note the / at the end is crucial
L=as.integer(shargs[3]) # number of admixing groups
firstind<-as.integer(shargs[4]); # which target individual to start from. If NUMA=2 then only this ind is run
NUMA=as.integer(shargs[5]) # total number of target admixed haplotypes 
MC=as.integer(shargs[6]) # number of cores to use for parallelized code
chrnos=1:22 # which chromosomes to run on
chrnos=21:22;firstind=1;NUMA=4;L=2;datasource="example_data/";target="Moroccan";ANC=NULL
chrnos=22:22;firstind=1;NUMA=2;L=2;datasource="HGDP/";target="simulated";RPE=0.0;ANC=T;
nchrno=length(chrnos) # number of chromosomes for these target haplotypes
HPC=2 # whether to use ff() chromosome-by-chromosome (HPC=1) or chromosomeXind-by-chromsomeXind(HPC=2) or not at all (HPC=F);
ffpath="/dev/shm/" # location of fast-files
if (!exists("PHASE")) PHASE=T
# the rest are mostly used in debugging, etc
if (!exists("ANC")) ANC=NULL; # no a-priori knowledge of which panels to use for which ancestries
verbose=T # print certain statements of progress as algorithm runs?
EM=T # run EM algorithm?
doMu=T # update copying matrix parameters?
doPI=T # update ancestry switching parameters parameters?
dorho=T # update recombination w/in same ancestry parameters? 
dotheta=T # update error / mutation parameters?

source("setup.R") # sets default parameters, sets up some objects required later, reads in data, and initialises model.
old.runtime<-as.numeric(Sys.time())
o.total=total
writelog<-function(t.logfile,t.alg,t.diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,t.cloglike) # single consistent function to write to EMlogfile
  write(file=t.logfile,c(t.alg,signif(t.diff.time,4),signif(t(t.Mu),4),signif(t.rho,4),c(sapply(t.PI, function(x) signif(t(x),4))),
			 sapply(t.alpha, function(x) signif(x,4)),sapply(t.lambda,function(x) round(x,4)),signif(t.theta,4),round(t.cloglike,4)),ncol=t.len,append=T)
eps=log(1.01) # i.e. a 1% increase in relative likelihood
Mu<-matrix(rep(1/kLL,L*kLL),kLL);for (ind in 1:NUMI) alpha[[ind]]=rep(1/L,L) # flatten out w.r.t. ancestry
runtime=NaN
# always need to run noanc.R b/c need good paras for init_Mu
tmp=fit_noanc_model(samp_chrnos, chrnos, NUMA, NUMP, kLL, L, KNOWN, label, umatch, G, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
		    prop.don, max.donors, maxmatch, maxmiss, initProb, d.w, t.w) 
transitions=tmp$t.transitionsmutmat=tmp$mutmat;Mu=tmp$Mu;theta=tmp$theta;rho=tmp$rho
ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
initProb=initprobs(T,NUMA,L,NUMP,kLL,PI,Mu,rho,alpha,label,NL)
runtime<-as.numeric(Sys.time())
if (kLL>L) # otherwise can't cluster kLL things into L clusters
{
  source("init_Mu.R")
  # use this to get #switches in noanc model w/o writelog
  tmp=all_donates(NUMI, Mu, alpha, kLL, PI, rho, lambda, theta, verbose=T, t.get_switches=T, max.donors, NUMP, G, umatch, maxmatchsize, d.w, 
		     t.w, gobs, flips, label, KNOWN, HPC, prethin=F, NUMA, nchrno, initProb, runtime, len,F,transitions,mutmat)
  ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;old.runtime=runtime=tmp$runtime;cloglike=tmp$cloglike;noanc_gswitches=tmp$noanc_gswitches
  windowed_copying<-window_chunks(nswitches=noanc_gswitches,ww=0.5,verbose=verbose) # similar in the same windows
  if (initonly)
  {
    save(file=paste0(resultsdir,"init_",target,"_",L, "way_", firstind, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,".RData"), 
	 target, panels, Mu, rho, theta, alpha, lambda, PI, windowed_copying, L, NUMA, nchrno, chrnos, g.loc, tol, dr, NL, kLL, 0)
    source("cleanup.R")
    stop("saving initialisation and quitting")
  }
  rm(noanc_gswitches) 
  tmp<-cluster_windows(windowed_copying,PLOT=F,t.L=L,verbose=verbose)
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
      source("phase_hunt.R") 
      if (M>0) 
	source("phase_mcmc.R") # now do MCMC re-phasing w/o output to console (ugly due to txtProgressBar)
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
final.flips=flips
source("ancunaware.R") # functions for getting localanc and gfbs that are ancestry unaware
EM=F;getnoancgfbs=T;eps=log(1.01);LOG=F;PLOT=F;
a.Mu=Mu;a.rho=rho;a.theta=theta;a.PI=PI;a.alpha=alpha;a.lambda=lambda
a.o.Mu=o.Mu;a.o.rho=o.rho;a.o.theta=o.theta;a.o.PI=o.PI;a.o.alpha=o.alpha;a.o.lambda=o.lambda
samp_chrnos=chrnos;subNUMA=NUMA;subNL=max(NL) # use them all
source("coancestry.R")

######### fully Mosaic curves with Mosaic phasing ############
gfbs=get_gfbs(NUMP, max.donors, donates, donatesl, donatesr, NUMA, L, G, kLL, transitions, umatch, maxmatchsize, d.w, t.w, gobs, mutmat, maxmiss, initProb, 
	      label, ndonors, flips)
source("localanc.R") 
if (verbose) cat("saving localanc results to file\n")
if (target!="simulated")
  save(file=paste0(resultsdir,"localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), localanc, final.flips, g.loc)
if (target=="simulated")
  save(file=paste0(resultsdir,"localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		   "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), localanc, g.true_anc, final.flips, g.loc)
save(file=paste0(resultsdir,"gfbs_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		 "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), gfbs)

if (verbose) cat("calculating ancestry aware re-phased coancestry curves\n"); acoancs=create_coancs(localanc,dr,"DIP");
######## GlobeTrotter style curves original phasing ##########
for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # undo phase flips
tmp=fit_noanc_model(samp_chrnos, chrnos, NUMA, NUMP, kLL, L, KNOWN, label, umatch, G, flips, gobs, PI, Mu, rho, theta, alpha, lambda, 
		    prop.don, max.donors, maxmatch, maxmiss, initProb, d.w, t.w) 
transitions=tmp$t.transitionsmutmat=tmp$mutmat;Mu=tmp$Mu;theta=tmp$theta;rho=tmp$rho
ndonors=tmp$ndonors;donates=tmp$donates;donatesl=tmp$donatesl;donatesr=tmp$donatesr;
source("cleanup.R")
Mu=a.Mu;rho=a.rho;theta=a.theta;PI=a.PI;alpha=a.alpha;lambda=a.lambda
o.Mu=a.o.Mu;o.rho=a.o.rho;o.theta=a.o.theta;o.PI=a.o.PI;o.alpha=a.o.alpha;o.lambda=a.o.lambda
noanc_unphased_localanc=get_ancunaware_localanc(NUMA,L,G,nchrno,noanc_gfbs,Mu,alpha) # works off noanc_gfbs
save(file=paste0(resultsdir,"noanc_unphased_localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		 "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), noanc_unphased_localanc, flips, g.loc)
if (verbose) cat("calculating ancestry unaware input phasing coancestry curves\n"); coancs=create_coancs(noanc_unphased_localanc,dr,"DIP")

if (verbose) cat("saving final results to file\n")
save(file=paste0(resultsdir,"",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",
		 GpcM,"_",prop.don,"_",max.donors,".RData"), target, phase.error.locs, o.Mu, o.lambda, o.theta, o.alpha, o.PI, o.rho, 
                 Mu, lambda, theta, alpha, PI, rho, L, NUMA, nchrno, chrnos, tol, dr, NL, kLL, acoancs, coancs)

