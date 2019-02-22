# function to simulate admixed genomes at a fixed number of generations ago from different donor genomes
admix_genomes=function(t.chrnos, t.ch, t.NUMA, t.NUMP, t.KNOWN, t.NN, t.multipanels, t.L, t.S, t.G, t.nl, t.kLL, t.NL, 
		       t.sim.alpha, t.sim.lambda, t.rates, t.g.map, t.dr, t.resultsdir, prop.missing=0,verbose=TRUE) {
  if (verbose) cat("creating admixed Chr ", t.chrnos[t.ch], "\n", sep="")

  d.w.ch=t.w.ch=list(u=list(),w=list())
  Y<-matrix(NA, t.NUMA, t.S[t.ch])
  true_anc.ch<-array(0,c(t.L,t.NUMA,t.S[t.ch]))

  # or simulate all breakpoints along all target chromosomes; then advance along positions, assigning last used donor to next breakpoint. 
  for (k in (t.NUMP+(1:t.NUMA))) # these are the admixed target haplotypes
  {
    tmpk=which({1:t.NN}[!t.KNOWN]==k)
    tmp=cumsum(c(1,ifelse((0:tmpk)%%2==0,1,3)))[tmpk]
    haps2=c(tmp,tmp+2) # the +2 ensures we don't use 2 haps from same donor ind
    ind=as.integer((tmpk+1)/2)
    tmps=0 # start at left of chromosome
    #tmp2k=tmpk # start with same index as the target
    tmp2k=haps2[1] # start with double the index of the target minus 1
    while (tmps[length(tmps)]<t.S[t.ch]) # while still in this chromosome
    {
      tmpia<-sample(1:t.L,1,prob=t.sim.alpha[[ind]]) # sample an ancestry 
      chunklengthM=rexp(1,t.sim.lambda[[ind]]) # in Morgans as per HapMix
      chunklengthM=round(chunklengthM/t.dr)*t.dr # to the nearest gridpoint
      RHS=which.min(abs(t.rates[tmps[length(tmps)]+1]+chunklengthM-t.rates)) # in units of the rates map; match to the genetic loci we have
      # make sure all SNPs later assigned to this gridpoint are switched together by taking rightmost SNP on this gridpoint
      RHS=max(which(t.g.map==t.g.map[RHS]))
      if (RHS==(t.S[t.ch]-1)) RHS=t.S[t.ch] # required as sometimes there's a gap of zero appended to the rates
      tmps=c(tmps, RHS)
      tmpil<-t.kLL+tmpia # use one panel in each anc
      l<-(tmps[length(tmps)-1]+1):(tmps[length(tmps)]) # vector of the markers in this ancestry window
      #cat(l[1],l[length(l)],"\n");readline()
      true_anc.ch[tmpia,tmpk,l]<-1 # write new true ancestry
      #tmp2k=(tmp2k+1-1)%%t.NUMA+2 # this adds 2 to tmp2k then shifts back to 1:t.NUMA; the 2 avoids haps of same donor ind 
      #tmp2k=(tmp2k+1-1)%%(t.NUMA*2)+2 # this adds one to tmp2k then shifts back to 1:(t.NUMA*2)
      tmp2k=haps2[haps2!=tmp2k] # use 2 donors per target; switch to one not last used
      Y[tmpk,l]<-as.integer(t.multipanels[[tmpil]][tmp2k,l]) 
    }
  }
  write.table(t(Y),file=paste0(t.resultsdir,"simulatedgenofile.",t.chrnos[t.ch],sep=""),row.names=F,col.names=F,sep="") # write out admixed individuals
  # add missing values
  if (prop.missing>0)
  {
    H=ifelse(t.NUMA>1, 2, 1)
    for (ind in 1:NUMI)
    {
      tmp=sample(1:t.S[t.ch],t.S[t.ch]*prop.missing) # both haps of an ind always missing together
      for (h in 1:H) 
      {
	k=(ind-1)*2+h
	Y[k,tmp]=9 # missing values indicated with a 9
      }
    }
  }

  k=1
  for (l in 1:t.kLL)
  {
    # do donors
    for (n in 1:t.NL[l])
    {
      tmpY<-as.integer(t.multipanels[[l]][n,]) 
      d.w.ch=cpp_unique_haps(tmpY,k,t.S[t.ch],t.G[t.ch],t.g.map-1,max(table(t.g.map)),d.w.ch)
      k<-k+1 # go to next one next
    }
    # now do targets
  }

  # now do targets
  l=t.kLL+1
  k=1
  for (n in 1:t.NL[l])
  {
    t.w.ch=cpp_unique_haps(Y[k,],k,t.S[t.ch],t.G[t.ch],t.g.map-1,max(table(t.g.map)),t.w.ch)
    k<-k+1 # go to next one next
  }
  return(list(d.w.ch=d.w.ch, t.w.ch=t.w.ch, true_anc.ch=true_anc.ch))
}

rdirichlet=function(n, t.alpha) {
  y=matrix(rgamma(n*length(t.alpha), t.alpha,1), ncol=length(t.alpha), byrow=TRUE)
  return(y/rowSums(y))
}

# some example simulations of admixture to try out using the HGDP dataset
create_sim=function(t.NUMA, t.L, t.o.lambda, mixers, t.panels, ratios=c(rdirichlet(1, rep(8,t.L))), fewer_ancs=NULL, verbose=TRUE) {
  sim.alpha<-sim.lambda<-list()
  NUMI<-max(1,t.NUMA/2)
  
  for (ind in 1:NUMI)
  {
    sim.alpha[[ind]]=c(rdirichlet(1, ratios*t.L*4)) # FLAG provide user control of this
    #sim.alpha[[ind]]=c(0.2,0.8) # leads to difficulties
    sim.lambda[[ind]]<-t.o.lambda
  }
  if (!is.null(fewer_ancs)) 
  {
    for (ind in fewer_ancs)
    {
      remanc=sample(1:t.L)[1:sample(1:(t.L-1),1)] # sample a subset of the t.L to remove from this individual, but never all of them
      sim.alpha[[ind]][remanc]=0
      if (t.L>2)
      {
	# never have ancs 1 and 2 both present
	if (sim.alpha[[ind]][1]) sim.alpha[[ind]][2]=0; if (sim.alpha[[ind]][2]) sim.alpha[[ind]][1]=0 
      }
      sim.alpha[[ind]]=sim.alpha[[ind]]/sum(sim.alpha[[ind]])
    } # remove first anc from first ind
    #for (ind in (1:NUMI)[-fewer_ancs]) # make those with this first anc only have a little
    #  {sim.alpha[[ind]][1]=0.1;sim.alpha[[ind]]=sim.alpha[[ind]]/sum(sim.alpha[[ind]])} # remove first anc from first ind
  }
  # will take in multipanels and simulates admixed individuals using the first panel in each anc
  # check that the supplied panels to mix are indeed members of panels 
  refs=t.panels[!(t.panels%in%mixers)]
  kLL=length(refs)
  if (verbose) cat("Admixing ", NUMI,  " individuals from ", paste(mixers, collapse=" and "),  " genomes ", sim.lambda, " generations ago\n", sep="")
  return(list(mixers=mixers, panels=refs, kLL=kLL, sim.alpha=sim.alpha, sim.lambda=sim.lambda))
}
