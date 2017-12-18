if (verbose) cat("creating admixed Chr ", chrnos[ch], "\n", sep="")
NL<-c(rep(nl,kLL),NUMA) # number of individuals in each group
for (k in 1:kLL)
  NL[k]<-min(NL[k], nrow(multipanels[[k]])) # make sure not to try to put too many...
LL=length(NL)
kLL=LL-1
NN<-sum(NL)
NUMP=NN-NUMA
label=rep(NaN,NN)
tmp<-c(0,cumsum(NL))
for (ll in 1:LL)
  label[(tmp[ll]+1):tmp[ll+1]]<-ll

# KNOWN are of known ancestry
KNOWN<-rep(F,NN)
KNOWN[label<LL]<-T # last one / group only not known
Y<-matrix(NA, NUMA, S[ch])
true_anc[[ch]]<-array(0,c(L,NUMA,S[ch]))

# FLAG or simulate all breakpoints along all target chromosomes; then advance along positions, assigning last used donor to next breakpoint. 
for (k in (NUMP+(1:NUMA))) # these are the admixed target haplotypes
{
  tmpk=which({1:NN}[!KNOWN]==k)
  tmp=cumsum(c(1,ifelse((0:tmpk)%%2==0,1,3)))[tmpk]
  haps2=c(tmp,tmp+2) # the +2 ensures we don't use 2 haps from same donor ind
  ind=as.integer((tmpk+1)/2)
  tmps=0 # start at left of chromosome
  #tmp2k=tmpk # start with same index as the target
  tmp2k=haps2[1] # start with double the index of the target minus 1
  while (tmps[length(tmps)]<S[ch]) # while still in this chromosome
  {
    tmpia<-sample(1:L,1,prob=sim.alpha[[ind]]) # sample an ancestry 
    chunklengthM=rexp(1,sim.lambda[[ind]]) # in Morgans as per HapMix
    chunklengthM=round(chunklengthM/dr)*dr # to the nearest gridpoint
    RHS=which.min(abs(rates[tmps[length(tmps)]+1]+chunklengthM-rates)) # in units of the rates map; match to the genetic loci we have
    # make sure all SNPs later assigned to this gridpoint are switched together by taking rightmost SNP on this gridpoint
    RHS=max(which(g.map==g.map[RHS]))
    if (RHS==(S[ch]-1)) RHS=S[ch] # required as sometimes there's a gap of zero appended to the rates
    tmps=c(tmps, RHS)
    tmpil<-kLL+tmpia # use one panel in each anc
    l<-(tmps[length(tmps)-1]+1):(tmps[length(tmps)]) # vector of the markers in this ancestry window
    #cat(l[1],l[length(l)],"\n");readline()
    true_anc[[ch]][tmpia,tmpk,l]<-1 # write new true ancestry
    #tmp2k=(tmp2k+1-1)%%NUMA+2 # this adds 2 to tmp2k then shifts back to 1:NUMA; the 2 avoids haps of same donor ind 
    #tmp2k=(tmp2k+1-1)%%(NUMA*2)+2 # this adds one to tmp2k then shifts back to 1:(NUMA*2)
    tmp2k=haps2[haps2!=tmp2k] # use 2 donors per target; switch to one not last used
    Y[tmpk,l]<-as.integer(multipanels[[tmpil]][tmp2k,l]) 
    if (OUTHAPMIX) hY[k,l]=Y[tmpk,l]
  }
}
write.table(t(Y),file=paste0(resultsdir,"simulatedgenofile.",chrnos[ch],sep=""),row.names=F,col.names=F,sep="") # write out admixed individuals
# add missing values
if (prop.missing>0)
{
  H=ifelse(NUMA>1, 2, 1)
  for (ind in 1:NUMI)
  {
    tmp=sample(1:S[ch],S[ch]*prop.missing) # both haps of an ind always missing together
    for (h in 1:H) 
    {
      k=(ind-1)*2+h
      Y[k,tmp]=9 # missing values indicated with a 9
    }
  }
}

k=1
for (l in 1:kLL)
{
  # do donors
  for (n in 1:NL[l])
  {
    tmpY<-as.integer(multipanels[[l]][n,]) 
    if (OUTHAPMIX) hY[k,]=tmpY
    d.w[[ch]]=cpp_unique_haps(tmpY,k,S[ch],G[ch],g.map-1,max(table(g.map)),d.w[[ch]])
    k<-k+1 # go to next one next
  }
  # now do targets
}

# now do targets
l=LL
k=1
for (n in 1:NL[l])
{
  t.w[[ch]]=cpp_unique_haps(Y[k,],k,S[ch],G[ch],g.map-1,max(table(g.map)),t.w[[ch]])
  k<-k+1 # go to next one next
}

