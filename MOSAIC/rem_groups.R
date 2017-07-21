source("calc_r2.R")
gfbs<-get_gfbs();source("localanc.R") 
buf=0e-6
old.cloglike<-cloglike<-NaN # so that next EM call doesn't immediately fail
init.kLL<-kLL
full.max.donors<-max.donors;
# first drop group with haps being used less than min.Mu.ratio/NUMP
rem.group<-which(apply(Mu/NL[1:kLL],1,max) < min.Mu.ratio/NUMP)  # max copying from each anc
# Given a true ancestry, we would like to see increase in r2 when a panel is removed and a (potential) decrease when one is kept (proposed and not removed)
# So count how many correct (remove+increase or keep+decrease) and how many wrong (remove+decrease or keep+increase) decisions are made
# (fine to remove and no change to r2)
count.correct=count.wrong=0
if (length(rem.group)>0)
{
  if (verbose) 
    cat("removing as seldom used:", panels[rem.group], "... ")
  if (target=="simulated")
  {
    r2.before<-round(dip_fr2(g.true_anc,localanc),4)
  }
  full.gfbs<-gfbs
  current.anc.r2<-round(dip_expected_fr2(localanc),4)
  source("drop_group_report.R") # calculates drop.anc.r2
  if (target=="simulated")
  {
    r2.after<-round(dip_fr2(g.true_anc, localanc),4)
    r2.diff<-round(r2.after-r2.before,4)
    if (r2.diff>=0) count.correct=count.correct+1 else count.wrong=count.wrong+1
  } 
  if (target=="simulated")
    cat("expected(r^2) = (", current.anc.r2, "->", drop.anc.r2, ") ", 
	"actual(r^2)=(", r2.before, "->", r2.after, ") ", ifelse(r2.diff<0,"(-)","(+)"),sep="", "\n")
  if (target!="simulated")
    cat("removing: expected(r^2) = (", current.anc.r2, "->", drop.anc.r2, ")\n", sep="") 
} else if (verbose) cat("none removed due to seldom used\n")
orig.anc.r2<-current.anc.r2<-round(dip_expected_fr2(localanc),4)
if (target=="simulated") 
  orig.r2<-round(dip_fr2(g.true_anc,localanc),4)


frdiff<-function()
{
  Pa<-Reduce("+",alpha)/NUMI 
  Pvga<-Mu 
  Pv=colSums(t(Pvga)*Pa) # sums to 1
  oldsum<-colSums(Pvga^2/Pv)
  rwith=rwithout=rep(0,kLL)
  for (v in 1:kLL)
  {
    newPvga<-Pvga;newPvga[v,]=0
    #newPvga<-t(t(newPvga)/(1-Pvga[v,])) # renormalisation # share out Pvga[v,] in proportion to remaining groups
    newPvga<-t(t(newPvga)/colSums(newPvga)) # renormalisation # share out Pvga[v,] in proportion to remaining groups
    newPv=colSums(t(newPvga)*Pa) # sums to 1
    newsum<-colSums(newPvga[-v,]^2/newPv[-v]) # -v to avoid division by zero
    for (a in 1:L)
    {
      rwith[v]=rwith[v]+Pa[a]/(1-Pa[a])*(oldsum[a]-1) # all the same obviously
      rwithout[v]=rwithout[v]+Pa[a]/(1-Pa[a])*(newsum[a]-1) # all the same obviously
    }
  }
  rwith=rwith/L
  rwithout=rwithout/L
  rwith-rwithout
}



expected_info_diff<-frdiff
infoloss<-expected_info_diff();names(infoloss)<-panels[1:kLL];infoloss<-sort(infoloss,index=T) 
proposal_scores<-NULL
if (!is.nan(drop.threshold))
{
  if (any(infoloss$x<drop.threshold,na.rm=T))
    cat("drop unclear panels: proposing the following groups:\n")
  rg=0
  while (any(infoloss$x<drop.threshold,na.rm=T))
  {
    rg=rg+1
    action<-NULL
    prop.groups<-panels[infoloss$ix[infoloss$x<drop.threshold]]
    if (dropfast)
    {
      infoloss$x[-c(1:min(maxdrops,kLL-L))]<-1 # don't ever drop more than maxdrops at a time and always leave at least L
      rem.group=match(names(infoloss$x)[which(infoloss$x<drop.threshold)],panels) # drop without checking 
    }
    if (!dropfast)
      rem.group=which(panels==prop.groups[1]) # take first one; panels may be updated during this loop
    if (verbose) 
      cat(panels[rem.group], "... ") # given in order of removal
    if (target=="simulated")
      r2.before<-round(dip_fr2(g.true_anc, localanc),4)
    # keep all-in things below in case not actually dropping this group
    full.Mu<-Mu;full.label<-label;full.LL<-LL;full.kLL<-kLL;full.KNOWN<-KNOWN;full.NL<-NL;full.NUMP<-NUMP;full.NN<-NN
    if (HPC) 
    {
      full.donates=full.donatesl=full.donatesr=list()
      if (prethin)
        {full.prethin_donates=full.prethin_donatesl=full.prethin_donatesr=list()}
      for(ch in 1:nchrno) 
      {
        full.donates[[ch]]=full.donatesl[[ch]]=full.donatesr[[ch]]=list()
        if (prethin)
          {full.prethin_donates[[ch]]=full.prethin_donatesl[[ch]]=full.prethin_donatesr[[ch]]=list()}
	for(ind in 1:NUMI) 
	{
	  open(donates[[ch]][[ind]])
	  full.donates[[ch]][[ind]]=clone(donates[[ch]][[ind]],filename=paste0(ffpath,target,"_full_donates_",ch,"_",ind,".ff"),overwrite=T)
	  close(donates[[ch]][[ind]]);close(full.donates[[ch]][[ind]])
	  open(donatesl[[ch]][[ind]])
	  full.donatesl[[ch]][[ind]]=clone(donatesl[[ch]][[ind]],filename=paste0(ffpath,target,"_full_donatesl_",ch,"_",ind,".ff"),overwrite=T)
	  close(donatesl[[ch]][[ind]]);close(full.donatesl[[ch]][[ind]])
	  open(donatesr[[ch]][[ind]])
	  full.donatesr[[ch]][[ind]]=clone(donatesr[[ch]][[ind]],filename=paste0(ffpath,target,"_full_donatesr_",ch,"_",ind,".ff"),overwrite=T)
	  close(donatesr[[ch]][[ind]]);close(full.donatesr[[ch]][[ind]])
	  if (prethin)
	  {
	    open(prethin_donates[[ch]][[ind]])
	    full.prethin_donates[[ch]][[ind]]=clone(prethin_donates[[ch]][[ind]],filename=paste0(ffpath,target,"_full_prethin_donates_",ch,"_",ind,".ff"),overwrite=T)
	    close(prethin_donates[[ch]][[ind]]);close(full.prethin_donates[[ch]][[ind]])
	    open(prethin_donatesl[[ch]][[ind]])
	    full.prethin_donatesl[[ch]][[ind]]=clone(prethin_donatesl[[ch]][[ind]],filename=paste0(ffpath,target,"_full_prethin_donatesl_",ch,"_",ind,".ff"),overwrite=T)
	    close(prethin_donatesl[[ch]][[ind]]);close(full.prethin_donatesl[[ch]][[ind]])
	    open(prethin_donatesr[[ch]][[ind]])
	    full.prethin_donatesr[[ch]][[ind]]=clone(prethin_donatesr[[ch]][[ind]],filename=paste0(ffpath,target,"_full_prethin_donatesr_",ch,"_",ind,".ff"),overwrite=T)
	    close(prethin_donatesr[[ch]][[ind]]);close(full.prethin_donatesr[[ch]][[ind]])
	  }
	}
      }
    }
    if (!HPC)
    {
      full.donates=donates
      full.donatesl=donatesl
      full.donatesr=donatesr
      if (prethin)
      {
        full.prethin_donates=prethin_donates
        full.prethin_donatesl=prethin_donatesl
        full.prethin_donatesr=prethin_donatesr
      }
    }
    full.d.w=d.w
    full.panels<-panels
    full.max.donors<-max.donors;
    full.initProb<-initProb;full.ndonors<-ndonors;full.localanc<-localanc
    full.gfbs<-gfbs
    full.lambda=lambda;full.Q=Q;full.alpha=alpha;full.theta=theta;full.rho=rho
    source("drop_group_report.R") # drop it and recalculate donates, fb, localanc, and drop.anc.r2
    if (target=="simulated")
    {
      r2.after<-round(dip_fr2(g.true_anc, localanc),4)
      r2.diff<-round(r2.after-r2.before,4)
      # correct is: (increase and drop) or (decrease and keep) wrong is anything else
      if ((r2.diff>=0&(current.anc.r2-drop.anc.r2 )<=buf) | (r2.diff<0&(current.anc.r2-drop.anc.r2)>buf)) 
	count.correct=count.correct+1 else count.wrong=count.wrong+1
    }
    if ((current.anc.r2-drop.anc.r2)>buf) # buffer to bias away from dropping groups
    {
      action<-c(action, "keeping")
      proposal_scores<-rbind(proposal_scores,c(rg, "keeping", prop.groups[1], infoloss$x[1]))
      # need to make sure just kept one is set to largest and greater than threshold so that we don't get stuck in a loop
      if (!dropfast)
      {
	infoloss<-infoloss$x[order(infoloss$ix)] # reverse index lookup
	names(infoloss)<-panels[1:kLL]
	infoloss[rem.group]=drop.threshold+max(infoloss)+1
	infoloss<-sort(infoloss,index=T)
      }
      if (dropfast)
	infoloss$x[]<-drop.threshold-min(infoloss$x)+infoloss$x # min is now drop.threshold and order is the same
      if (verbose) 
      {
	if (target=="simulated")
	  cat("keeping: expected(r^2)=(", current.anc.r2, "->", drop.anc.r2, ") ", 
	  "actual(r^2)=(", r2.before, "->", r2.after, ") ", ifelse(r2.diff<0,"(+)","(-)"),sep="", "\n")
	if (target!="simulated")
	  cat("keeping: expected(r^2) =(", current.anc.r2, "->", drop.anc.r2, ")\n", sep="")
      }
      Mu<-full.Mu;label<-full.label;LL<-full.LL;kLL<-full.kLL;KNOWN<-full.KNOWN;NL<-full.NL;NUMP<-full.NUMP;NN<-full.NN # put all back in
      if (HPC)
      {
        donates=donatesl=donatesr=list()
	if (prethin)
          {prethin_donates=prethin_donatesl=prethin_donatesr=list()}
	for(ch in 1:nchrno) 
	{
          donates[[ch]]=donatesl[[ch]]=donatesr[[ch]]=list()
	  if (prethin)
            {prethin_donates[[ch]]=prethin_donatesl[[ch]]=prethin_donatesr[[ch]]=list()}
	  for(ind in 1:NUMI) 
	  {
	    open(full.donates[[ch]][[ind]])
  	    donates[[ch]][[ind]]=clone(full.donates[[ch]][[ind]],filename=paste0(ffpath,target,"_donates_",ch,"_",ind,".ff"),overwrite=T)
	    close(donates[[ch]][[ind]]);close(full.donates[[ch]][[ind]])
	    open(full.donatesl[[ch]][[ind]])
	    donatesl[[ch]][[ind]]=clone(full.donatesl[[ch]][[ind]],filename=paste0(ffpath,target,"_donatesl_",ch,"_",ind,".ff"),overwrite=T)
	    close(donatesl[[ch]][[ind]]);close(full.donatesl[[ch]][[ind]])
	    open(full.donatesr[[ch]][[ind]])
	    donatesr[[ch]][[ind]]=clone(full.donatesr[[ch]][[ind]],filename=paste0(ffpath,target,"_donatesr_",ch,"_",ind,".ff"),overwrite=T)
	    close(donatesr[[ch]][[ind]]);close(full.donatesr[[ch]][[ind]])
	    if (prethin)
	    {
	      open(full.prethin_donates[[ch]][[ind]])
  	      prethin_donates[[ch]][[ind]]=clone(full.prethin_donates[[ch]][[ind]],filename=paste0(ffpath,target,"_prethin_donates_",ch,"_",ind,".ff"),overwrite=T)
	      close(prethin_donates[[ch]][[ind]]);close(full.prethin_donates[[ch]][[ind]])
	      open(full.prethin_donatesl[[ch]][[ind]])
	      prethin_donatesl[[ch]][[ind]]=clone(full.prethin_donatesl[[ch]][[ind]],filename=paste0(ffpath,target,"_prethin_donatesl_",ch,"_",ind,".ff"),overwrite=T)
	      close(prethin_donatesl[[ch]][[ind]]);close(full.prethin_donatesl[[ch]][[ind]])
	      open(full.prethin_donatesr[[ch]][[ind]])
	      prethin_donatesr[[ch]][[ind]]=clone(full.prethin_donatesr[[ch]][[ind]],filename=paste0(ffpath,target,"_prethin_donatesr_",ch,"_",ind,".ff"),overwrite=T)
	      close(prethin_donatesr[[ch]][[ind]]);close(full.prethin_donatesr[[ch]][[ind]])
	    }
	  }
	}
      }
      if (!HPC)
      {
	donates=full.donates;donatesl=full.donatesl;donatesr=full.donatesr;
	if (prethin)
	  {prethin_donates=full.prethin_donates;prethin_donatesl=full.prethin_donatesl;prethin_donatesr=full.prethin_donatesr;}
      }
      panels<-full.panels
      d.w=full.d.w
      max.donors<-full.max.donors
      initProb<-full.initProb;ndonors<-full.ndonors;localanc<-full.localanc
      gfbs<-full.gfbs
      lambda=full.lambda;Q=full.Q;alpha=full.alpha;theta=full.theta;rho=full.rho
      mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
      for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,Q[[ind]],Mu,rho,NL)
      rem.group<-NULL
    } else 
    {
      action<-c(action,"removing")
      proposal_scores<-rbind(proposal_scores,c(rg, "removing", prop.groups[1], infoloss$x[1]))
      if (PLOT)
      {
	par(mfrow=c(1,2))
	ord.Mu<-plot_Mu(full.Mu,full.alpha,full.NL,MODE="joint"); abline(h=which(rownames(ord.Mu)==full.panels[rem.group]))
	ord.Mu<-plot_Mu(Mu,alpha,NL,MODE="joint")
      }
      if (verbose) 
      {
	if (target=="simulated")
	  cat("removing: expected(r^2) = (", current.anc.r2, "->", drop.anc.r2, ") ", 
	      "actual(r^2)=(", r2.before, "->", r2.after, ") ", ifelse(r2.diff<0,"(-)","(+)"),sep="", "\n")
	if (target!="simulated")
	  cat("removing: expected(r^2) = (", current.anc.r2, "->", drop.anc.r2, ")\n", sep="") 
      }
      current.anc.r2<-drop.anc.r2
      # recalculate infoloss for all remaining groups
      infoloss<-expected_info_diff();names(infoloss)<-panels[1:kLL];infoloss<-sort(infoloss,index=T) 
      # use remaining infoloss w/o recalculation
      #infoloss<-infoloss$x[order(infoloss$ix)];infoloss<-infoloss[-rem.group];infoloss<-sort(infoloss,index=T) # reverse index lookup
    }
    if (HPC) for (ch in 1:nchrno) for(ind in 1:NUMI) 
    {
      delete(full.donates[[ch]][[ind]])
      delete(full.donatesl[[ch]][[ind]])
      delete(full.donatesr[[ch]][[ind]])
      if (prethin)
      {
        delete(full.prethin_donates[[ch]][[ind]])
        delete(full.prethin_donatesl[[ch]][[ind]])
        delete(full.prethin_donatesr[[ch]][[ind]])
    }
    }
    rm(full.Mu,full.label,full.LL,full.kLL,full.KNOWN,full.NL,full.NUMP,full.NN,full.d.w,full.panels,
       full.max.donors,full.initProb,full.donates,full.donatesl,full.donatesr,full.ndonors,full.localanc,full.gfbs)
      if (prethin)
	rm(full.prethin_donates,full.prethin_donatesl,full.prethin_donatesr)
  }
  if (!dropfast)
  {
    proposal_scores<-as.data.frame(proposal_scores)
    colnames(proposal_scores)<-c("index","action","panel","infoloss")
    proposal_scores$index<-as.numeric(as.character(proposal_scores$index))
    proposal_scores$infoloss<-as.numeric(as.character(proposal_scores$infoloss))
    cat("Performance of proposal mechanism for group dropping: ")
    print(wilcox.test(proposal_scores$index~proposal_scores$action))
    #print(t.test(proposal_scores$infoloss~proposal_scores$action))
  }
  final.anc.r2<-dip_expected_fr2(localanc)
  if (verbose)
  {
    if ((target=="simulated") & (count.correct+count.wrong)>0)
    {
      final.r2<-round(dip_fr2(g.true_anc,localanc),4)
      #cat("made", count.correct, "correct calls and", count.wrong, "wrong calls based on r^2 with truth\n")
      cat("total change in expected r^2:", final.anc.r2-orig.anc.r2, "total change in r^2 with truth:", final.r2-orig.r2, "\n")
    } else cat("total change in expected r^2:", final.anc.r2-orig.anc.r2, "\n")
  }
}
if (kLL<init.kLL) source("create_logfile.R")
if (LOG)
{
  if (dropfast) source("klikelihood.R") # one forward pass to find new log-likelihood
  runtime<-as.numeric(Sys.time());diff.time<-runtime-old.runtime;old.runtime<-runtime;
  writelog("dropgroup")
}

