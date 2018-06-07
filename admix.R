# function to simulate admixed genomes at a fixed number of generations ago from different donor genomes
admix_genomes=function(t.chrnos, t.ch, t.NUMA, t.NUMP, t.KNOWN, t.NN, t.multipanels, t.L, t.S, t.G, t.nl, t.kLL, t.sim.alpha, t.sim.lambda,
		       prop.missing=0) {
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
      chunklengthM=round(chunklengthM/dr)*dr # to the nearest gridpoint
      RHS=which.min(abs(rates[tmps[length(tmps)]+1]+chunklengthM-rates)) # in units of the rates map; match to the genetic loci we have
      # make sure all SNPs later assigned to this gridpoint are switched together by taking rightmost SNP on this gridpoint
      RHS=max(which(g.map==g.map[RHS]))
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
  write.table(t(Y),file=paste0(resultsdir,"simulatedgenofile.",t.chrnos[t.ch],sep=""),row.names=F,col.names=F,sep="") # write out admixed individuals
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
    for (n in 1:NL[l])
    {
      tmpY<-as.integer(t.multipanels[[l]][n,]) 
      d.w.ch=cpp_unique_haps(tmpY,k,t.S[t.ch],t.G[t.ch],g.map-1,max(table(g.map)),d.w.ch)
      k<-k+1 # go to next one next
    }
    # now do targets
  }

  # now do targets
  l=t.kLL+1
  k=1
  for (n in 1:NL[l])
  {
    t.w.ch=cpp_unique_haps(Y[k,],k,t.S[t.ch],t.G[t.ch],g.map-1,max(table(g.map)),t.w.ch)
    k<-k+1 # go to next one next
  }
  return(list(d.w.ch=d.w.ch, t.w.ch=t.w.ch, true_anc.ch=true_anc.ch))
}


# some example simulations of admixture to try out using the HGDP dataset
require(gtools) # for rdirichlet
example_sims=function(t.NUMA, t.L, t.o.lambda) {
  sim.alpha<-sim.lambda<-list()
  NUMI<-max(1,t.NUMA/2)
  if (!exists("mean.sim.alpha")) mean.sim.alpha=c(rdirichlet(1, rep(8,t.L)))
  #mean.sim.alpha<-c(0.9,0.1) # manual choice for t.L=2
  #mean.sim.alpha<-c(0.1,0.1,0.8) # manual choice for t.L=3
  for (ind in 1:NUMI)
  {
    sim.alpha[[ind]]=c(rdirichlet(1, mean.sim.alpha*t.L*4))
    #sim.alpha[[ind]]=c(0.2,0.8) # leads to difficulties
    sim.lambda[[ind]]<-t.o.lambda
  }
  #sim.alpha[[1]]=c(1,0,0);sim.alpha[[2]]=c(0,1,0)#;sim.alpha[[3]]=c(0,0,1); # one pure from each anc=pop; useful when theta and rho not updated
  #sim.alpha[[1]]=c(1,1,1);sim.alpha[[2]]=c(1,1,1); # two that are equal parts each anc=pop; useful when theta and rho not updated
  if (!exists("fewer_ancs")) fewer_ancs=NULL
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
  # takes in multipanels and simulates admixed individuals using the first panel in each anc
  if (!is.null(ANC)) # if set to TRUE then use one of these examples, depending on choice of t.L
  {
    if (ANC[1]==T) # if just set to TRUE but mixing panels not supplied directly
    {
      if (t.L==2)
      {
	#pops<-list(Eu=c("Iranian", "Lezgin", "Armenian", "Sindhi", "Brahui", "Georgian"),
	#	   As=c("Hezhen", "HanNchina", "Han", "Tu", "Oroqen", "Daur", "Xibo", "Tujia", "Yakut"))
	#mixers=c("Pathan","Mongola")
	pops<-list(Eu=c("Spanish","Norwegian","Polish","English","Welsh","Ireland","Scottish"),
		   Af=c("Mandenka","Hadza","BantuKenya","BantuSouthAfrica","Sandawe","BiakaPygmy","MbutiPygmy"))
	mixers=c("French","Yoruba")
      }
      if (t.L==3)
      {
	pops<-list(Eu=c("English", "Scottish", "Spanish", "Polish", "Melanesian", "Turkish", "Norwegian", "Hungarian"),
		   #As=c("Han", "Yi", "Daur", "She", "Maya"),# "HanNchina", "Hezhen", "Oroqen", "Tu", # easier
		   SEA=c("Sindhi", "Indian", "Makrani", "IndianJew", "Yemeni"), # harder 
		   Af=c("Mandenka", "BantuKenya", "Sandawe", "BantuSouthAfrica", "Pima", "Surui", "MbutiPygmy", "Hadza")
		   )
	mixers=c("French", "Pathan", "Yoruba") # harder
	#mixers=c("French", "Japanese", "Yoruba") # easier
      }
      if (t.L==4)
      {
	pops<-list(Eu=c("English", "Scottish", "Polish", "Spanish"),
		   As=c("Han", "Yi", "Daur", "She"),
		   Af=c("Mandenka","BantuKenya", "Sandawe"),
		   Sa=c("Maya", "Surui"))
	mixers=c("French", "Japanese", "Yoruba", "Pima")
      }
    }
    if (ANC[1]!=T) # something supplied
    { 
      mixers=ANC
      # check that the supplied panels to mix are indeed members of panels and that length(ANC)==t.L
      refs=sample(panels[!(panels%in%mixers)]);tmp=sample(1:t.L, length(refs), replace=T); # random order of unused panels, then split into approx equal groups
      pops=list()
      for (l in 1:t.L) pops[[l]]=c(ANC[l],refs[tmp==l]) # first panel is always the panel used to simulate the admixed individuals, remaining panels taken at random
    }
  }
  if (is.null(ANC))
  {
    ANC=mixers=sample(panels,t.L)
    if (verbose)
      cat("no ancestral panels specified so using:", ANC, "\n")
    refs=sample(panels[!(panels%in%mixers)]);tmp=sample(1:t.L, length(refs), replace=T); # random order of unused panels, then split into approx equal groups
    pops=list()
    for (l in 1:t.L) pops[[l]]=c(ANC[l],refs[tmp==l]) # # first panel is always the panel used to simulate the admixed individuals, remaining panels taken at random
  }
  pops[[t.L+1]]="simulated";names(pops)[t.L+1]="target"
  panels<-unlist(pops)
  kLL=length(panels)-1 # remove simulated as there is no such panel
  panels=panels[1:kLL]
  if (verbose) cat("Creating ", NUMI, " simulated ", t.L, "-way admixed target individuals with ", kLL, " panels\n", sep="")
  ANC<-NULL
  for (i in 1:t.L)
    ANC[i]=mixers[i] 
  return(list(ANC=ANC, mixers=mixers, panels=panels, kLL=kLL, sim.alpha=sim.alpha, sim.lambda=sim.lambda))
}
