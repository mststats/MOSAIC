# function that reads in the data and lays on a grid along recombination rates map
read_panels=function(datasource, t.nchrno, t.nl, t.FLAT, dr, t.o.lambda, t.resultsdir, mask=NULL, S=rep(NaN,t.nchrno),
		     firstind=1) {
  panels<-read.table(paste(datasource,"sample.names",sep=""), header=F);panels<-as.character(unique(panels[,1]))
  gobs=g.loc=list()
  maxmatch=maxmiss=0
  if (!is.null(mask))
  {
    maskpanels=NULL;for (maskpanel in mask) {tmp2=match(maskpanel,panels);if (length(tmp2)>0) maskpanels=c(maskpanels,tmp2)}
    panels<-panels[-maskpanels] # remove masked groups 
  }
  tmp=match(target,panels); if (!is.na(tmp)) panels=panels[-tmp] # remove target panel
  kLL=length(panels)
  if (is.null(ANC))
  {
    tmp<-panels
    tmp<-tmp[tmp!=target] # take everything else except the target as a potential donor panel
    panels=as.character(tmp)
    kLL=length(panels)
    if (target!="simulated") panels[kLL+1]=target
  }
  if (target=="simulated") true_anc<-g.true_anc<-list()
  if (!is.null(ANC) | target=="simulated") # always call this if looking at a simulation and / or if ANC=T
  {
    # example simulations
    tmp=example_sims(NUMA, L, t.o.lambda) # note that this may use reduced set of panels
    ANC=tmp$ANC;mixers=tmp$mixers;panels=tmp$panels;kLL=tmp$kLL;sim.alpha=tmp$sim.alpha;sim.lambda=tmp$sim.lambda
  }
  d.w=list() # map to unique donor  haps at each gridpoint
  t.w=list() # map to unique target haps at each gridpoint
  umatch=list() # lookup for donor,target to #matches using t.w and d.w
  maxmatchsize=NULL # maximum size of the umatch matrix at any gridpoint on each chromosome
  G=NULL
  for (ch in 1:t.nchrno)
  {
    multipanels<-list() # shouldn't read all of these in at once!
    for (i in 1:length(panels))
    {
      tmp<-scan(paste0(datasource,panels[i],"genofile.",chrnos[ch]),what="character",quiet=T)
      tmp<-strsplit(tmp,"")
      allS<-length(tmp)
      N2<-length(tmp[[1]])
      multipanels[[i]]<-matrix(sapply(tmp, as.double), N2, allS)
    }
    if (target=="simulated")
      for (i in 1:L)
      {
	tmp<-scan(paste0(datasource,mixers[i],"genofile.",chrnos[ch]),what="character",quiet=T)
	tmp<-strsplit(tmp,"")
	allS<-length(tmp)
	N2<-length(tmp[[1]])
	multipanels[[kLL+i]]<-matrix(sapply(tmp, as.double), N2, allS)
      }
    if (target!="simulated") names(multipanels)<-panels
    if (target=="simulated") names(multipanels)<-c(panels,mixers)
    snps<-read.table(paste0(datasource,"snpfile.",chrnos[ch])) 
    #### inputs #################################
    if (is.nan(S[ch])) S[ch]<-nrow(snps)
    #if (datasource=="HGDP_PEL/") snps<-cbind(snps[,1],snps[,2],snps[,3]/max(snps[,3])/2,snps[,3],snps[,4],snps[,5])
    if (S[ch]<nrow(snps))
    {
      # only look at S loci
      #locs<-sort(sample(1:nrow(snps),S)) # at random
      locs<-seq(1,nrow(snps),by=as.integer(nrow(snps)/S[ch])) # take a subset
      #locs<-1:S[ch] # first S[ch] only
      S[ch]<-length(locs) 
    }
    if (S[ch]==nrow(snps)) locs<-1:S[ch]
    snps<-snps[locs,]
    for (i in 1:length(multipanels)) multipanels[[i]]<-multipanels[[i]][,locs]
    
    # start at firstind i.e. remove all before this in target pop. 
    if (firstind!=1)
    {
      i=length(multipanels) # Admixed target is always stored last
      multipanels[[i]]<-multipanels[[i]][-(1:(2*(firstind-1))),]
    }
    NL<-c(rep(t.nl,kLL),NUMA)
    LL<-kLL+1
    for (i in 1:kLL)
      NL[i]<-min(NL[i],nrow(multipanels[[i]])) # make sure none are asked for more than they have
    if (target!="simulated")
      NL[LL]<-min(NL[LL],nrow(multipanels[[LL]])) # make sure none are asked for more than they have
    NN<-sum(NL)
    NUMA=NL[LL] # to make sure it doesn't go over available target haps
    NUMI=max(1,NUMA/2)
    label=rep(NaN,NN)
    tmp<-c(0,cumsum(NL))
    for (ll in 1:LL)
      label[(tmp[ll]+1):tmp[ll+1]]<-ll
    # KNOWN are of known ancestry
    KNOWN<-rep(F,NN)
    KNOWN[label<LL]<-T # last one / group only not known
    NUMP<-sum(KNOWN) # i.e. number in panels
    NUMA<-sum(!KNOWN) #i.e. number of targets / admixed

    all_rates<-matrix(scan(paste0(datasource,"rates.",chrnos[ch]),skip=1,quiet=T),ncol=2)
    locs<-as.integer(snps[,4])
    tmp=match(locs, all_rates[,1])
    rates=all_rates[tmp,2] # use ones with hap data; some may be missing if in snps file but not in rates file
    # rates are flat in sections so use rate to the left if missing
    tmp=which(is.na(tmp))
    rates[tmp]=all_rates[vapply(tmp,function(l) which.max(all_rates[all_rates[,1]<locs[l],1]-locs[l]),0L),2]
    rates<-rates/100 # /100 to move to morgans from centimorgans 
    if (t.FLAT) 
    {
      rates<-seq(rates[1],2*rates[S[ch]],l=S[ch])
      warning("using flat recombination rates map",immediate.=T)
    }
    rm(all_rates)
    # G will be correct / consistent w/in 0.5 and is large so almost exactly the same dr across chromosomes
    G[ch]<-as.integer((rates[S[ch]]-rates[1])/dr+1)
    g.rates<-seq(rates[1],rates[S[ch]],l=G[ch])
    g.map<-vapply(1:S[ch], function(s) which.min((rates[s]-g.rates)^2),0L) # create map from rates to grid
    if (target!="simulated")
    {
      d.w[[ch]]=t.w[[ch]]=list(u=list(),w=list())
      k=1
      for (l in 1:kLL)
      {
	# do mixers
	for (n in 1:NL[l])
	{
	  Y<-as.integer(multipanels[[l]][n,]) # take whole of the nth haplotype here
	  d.w[[ch]]=cpp_unique_haps(Y,k,S[ch],G[ch],g.map-1,max(table(g.map)),d.w[[ch]])
	  k<-k+1 # go to next one next
	}
	# now do targets
      }
      l=LL
      k=1
      for (n in 1:NL[l])
      {
	Y<-as.integer(multipanels[[l]][n,]) # take whole of the nth haplotype here
	t.w[[ch]]=cpp_unique_haps(Y,k,S[ch],G[ch],g.map-1,max(table(g.map)),t.w[[ch]])
	k<-k+1 # go to next one next
      }
    }
    if (target=="simulated")
    {
      tmp=admix_genomes(chrnos, ch, NUMA, NUMP, KNOWN, NN, multipanels, L, S, G, t.nl, kLL, NL, sim.alpha, sim.lambda, rates, g.map, 
			dr, t.resultsdir, panels)
      d.w[[ch]]=tmp$d.w.ch
      t.w[[ch]]=tmp$t.w.ch
      true_anc[[ch]]=tmp$true_anc.ch
    }
    umatch[[ch]]=create_umatch(d.w[[ch]],t.w[[ch]],g.map,G[ch])
    maxmatchsize[ch]=max(sapply(umatch[[ch]],function(x) prod(dim(x))))
    # don't need to keep the haps any more
    d.w[[ch]]$u=NULL;d.w[[ch]]$du=NULL
    t.w[[ch]]$u=NULL;t.w[[ch]]$du=NULL
    rm(multipanels)
    tmp=create_grid(G[ch],S[ch],g.map, chrnos[ch], NUMA, L, umatch[[ch]], t.w[[ch]], true_anc[[ch]], locs, NUMI)
    g.loc[[ch]]=tmp$g.loc;gobs[[ch]]=tmp$gobs;maxmatch_chr=tmp$maxmatch_chr;maxmiss_chr=tmp$maxmiss_chr 
    if (target=="simulated") g.true_anc[[ch]]=tmp$g.true_anc_chr
    maxmatch<-max(maxmatch,maxmatch_chr)
    maxmiss<-max(maxmiss,maxmiss_chr)
    rm(snps,rates) # leave in if planning to re-grid
  }

  if (verbose) cat("Fitting model to ", NUMI, " ", target, " ", L, "-way admixed target individuals using ", kLL, " panels\n", sep="")
  ans=list(maxmatch=maxmatch, maxmiss=maxmiss,g.loc=g.loc,gobs=gobs,d.w=d.w,t.w=t.w,umatch=umatch,NUMP=NUMP,LL=LL,NUMI=NUMI,
	   label=label,KNOWN=KNOWN,NN=NN,kLL=kLL,NL=NL,G=G,maxmatchsize=maxmatchsize,panels=panels)
  if (target=="simulated") 
    ans$g.true_anc=g.true_anc
    return(ans)
}
