# script that reads in the data and lays on a grid along recombination rates map
source("compressed_grid.R")
panels<-read.table(paste(datasource,"sample.names",sep=""), header=F);panels<-as.character(unique(panels[,1]))
if (!exists("mask")) mask=NULL
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
  pops<-list(panels=tmp, AD=target)
  panels=as.character(tmp)
  kLL=length(panels)
  if (target!="simulated") panels[kLL+1]=target
}
if (target=="simulated") true_anc<-g.true_anc<-list()
if (!is.null(ANC) | target=="simulated") # always call this if looking at a simulation and / or if ANC=T
  source("examples.R") # example simulations and real admixed panels
d.w=list() # map to unique donor  haps at each gridpoint
t.w=list() # map to unique target haps at each gridpoint
umatch=list() # lookup for donor,target to #matches using t.w and d.w
maxmatchsize=NULL # maximum size of the umatch matrix at any gridpoint on each chromosome
for (ch in 1:nchrno)
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
  if (!exists("firstind")) firstind=1
  # start at firstind i.e. remove all before this in target pop. 
  if (firstind!=1)
  {
    i=length(multipanels) # Admixed target is always stored last
    multipanels[[i]]<-multipanels[[i]][-(1:(2*(firstind-1))),]
  }
  NL<-c(rep(nl,kLL),NUMA)
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
  if (FLAT) 
  {
    rates<-seq(rates[1],2*rates[S[ch]],l=S[ch])
    warning("using flat recombination rates map",immediate.=T)
  }
  rm(all_rates)
  G[ch]<-as.integer((rates[S[ch]]-rates[1])/dr+1)
  g.rates<-seq(rates[1],rates[S[ch]],l=G[ch])
  g.map<-vapply(1:S[ch], function(s) which.min((rates[s]-g.rates)^2),0L) # create map from rates to grid
  d.w[[ch]]=t.w[[ch]]=list(u=list(),w=list())
  if (target!="simulated")
  {
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
    source("admix.R")
  umatch[[ch]]=create_umatch(d.w[[ch]],t.w[[ch]],g.map,G[ch])
  maxmatchsize[ch]=max(sapply(umatch[[ch]],function(x) prod(dim(x))))
  # don't need to keep the haps any more
  d.w[[ch]]$u=NULL;d.w[[ch]]$du=NULL
  t.w[[ch]]$u=NULL;t.w[[ch]]$du=NULL
  rm(multipanels)
  source("grid.R")
  rm(snps,rates) # leave in if planning to re-grid
}

if (verbose) cat("Fitting model to ", NUMI, " ", target, " ", L, "-way admixed target individuals using ", kLL, " panels\n", sep="")
