# if we have a version impervious to phasing it will last longer without changing during thin / phase / EM
require(compiler)
# code to select enough donors to capture prop.don of the copying 
#sourceCpp("switches.cpp") # gives cppswitches
create_donates<-function(getswitches,ch,ind,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.flips,t.kLL,t.Mu,t.rho,t.theta,HPC,prethin=F)
{
  THIN=F
  if (NUMA==1) H=1 else H=2
  hap<-c(ind*2-1,ind*2)
  # we need these (temporarily) to calculate E[switches] and (bizarrely) to calculate the thinned version of themselves
  t.donates<-0:(NUMP-1) # Holds vectors to indicate which donors are copied from at each gridpoint. 
  t.donatesl<-t.donatesr<-t.donates # left and right one locations of aligned indices
  t.ndonors<-rep(NUMP,G[ch]) # too large; should limit this really
  if (max.donors==NUMP & !getswitches) 
    return(list(ndonors=t.ndonors,donates=t.donates,donatesl=t.donatesl,donatesr=t.donatesr)) # these are vectors; dealt with automatically in cpp calls
  # note that we always need this when thinning; see probmass below
  probmass<-matrix(0,G[ch],NUMP)
  switches<-list() # returns empty list if not needed; serves as placeholder
  # fit a HMM with no latent ancestry
  L<-1;tmpPI<-matrix(0,1,1);
  noanc_initProb<-matrix(0,H,kLL); for (h in 1:H) noanc_initProb[h,]<-t.Mu/NL[1:kLL]
  noanc_mutmat<-fmutmat(t.theta, L, maxmiss, maxmatch)
  noanc_transitions<-s_trans(L,t.kLL,tmpPI,t.Mu,t.rho,NL)
  noanc_fors<-noanc_sumfors<-noanc_backs<-noanc_scalefactor<-noanc_scalefactorb<-list()
  for (h in 1:H) { 
    noanc_fors[[h]]<-rep(0,G[ch]*NUMP)
    noanc_sumfors[[h]]<-matrix(0,G[ch],L)
    noanc_backs[[h]]<-rep(0,G[ch]*NUMP)
    noanc_scalefactor[[h]]<-rep(0,G[ch])
    noanc_scalefactorb[[h]]<-rep(0,G[ch])
  }
  if (prethin & !getswitches) # only need to use this if pre-thinning and not finding switch counts based on all donors
  {
    THIN=T
    t.ndonors=prethin_ndonors[[ch]][[ind]]
    if (HPC)
    {
      open(prethin_donates[[ch]][[ind]]);open(prethin_donatesl[[ch]][[ind]]);open(prethin_donatesr[[ch]][[ind]]);
      t.donates=prethin_donates[[ch]][[ind]][];t.donatesl=prethin_donatesl[[ch]][[ind]][];t.donatesr=prethin_donatesr[[ch]][[ind]][]
      close(prethin_donates[[ch]][[ind]]);close(prethin_donatesl[[ch]][[ind]]);close(prethin_donatesr[[ch]][[ind]]);
    }
    if (!HPC)
    {
      t.donates=prethin_donates[[ch]][[ind]];t.donatesl=prethin_donatesl[[ch]][[ind]];t.donatesr=prethin_donatesr[[ch]][[ind]]
    }
  }
  for (h in 1:H)
  {
    k=hap[h]
    cppforward(k,NUMA,NUMP,THIN,NUMP,kLL,L,0,G[ch],G[ch],noanc_transitions,t.umatch,t.maxmatchsize,d.w[[ch]],t.w[[ch]],t.gobs,noanc_mutmat,maxmiss,noanc_initProb[h,],
	       label,t.ndonors,t.donates,t.donatesl,t.flips,noanc_fors[[h]],noanc_sumfors[[h]],noanc_scalefactor[[h]])
    cppbackward(k,NUMA,NUMP,THIN,NUMP,L,0,G[ch],G[ch],noanc_transitions,t.umatch,t.maxmatchsize,d.w[[ch]],t.w[[ch]],t.gobs,noanc_mutmat,maxmiss,
		label,t.ndonors,t.donates,t.donatesr,t.flips,noanc_backs[[h]],noanc_scalefactorb[[h]])
    probmass<-probmass+cppforback(NUMP,THIN,NUMP,L,G[ch],t.ndonors,t.donates,noanc_fors[[h]],noanc_backs[[h]]) # 1 as first argument b/c using all now
    if (getswitches) 
      switches[[h]]<-t(matrix(cppswitches(h,NUMA,NUMP,THIN,NUMP,G[ch],NL,label,noanc_sumfors[[h]],noanc_backs[[h]],
    				 noanc_transitions,t.flips,noanc_mutmat,maxmiss,t.umatch,t.maxmatchsize,t.dw,t.tw,t.gobs,t.ndonors,t.donates)$switches,NUMP))
  }
  if (max.donors==NUMP & getswitches) 
    return(list(ndonors=t.ndonors,donates=t.donates,donatesl=t.donatesl,donatesr=t.donatesr,switches=switches)) 
  #print(range(t.ndonors));print(range(noanc_fors[[1]]));readline()
  # next lines tested and absolutely required
  if (NUMA>1)
    for (h in 1:2) # both haps fors and backs must be calculated before the next line is called
      probmass<-probmass+cppforback(NUMP,THIN,NUMP,L,G[ch],t.ndonors,t.donates,noanc_fors[[h]],noanc_backs[[h+ifelse(h%%2,1,-1)]]) # 1 as first argument b/c using all now
  # if NUMA==1 probmass is now f[1]*b[1] 
  # if NUMA>1 probmass is now f[1]*b[1], f[2]*b[2], f[1]*b[2], f[2]*b[1]
  probmass<-t(t(probmass)/rowSums(probmass))
  f<-function(x)
  {
     tmp<-order(x,decreasing=T)[1:max.donors] # descending order of donors; can't use partial sorting as index can't be returned
     quants<-cumsum(x[tmp]);quants[max.donors]=1
     donors<-tmp
     cutoff<-max(which(quants>prop.don)[1], min.donors) # index of first one to go over prop.don but take at least min.donors
     donors<-donors[1:max.donors] # easier to return also the not needed ones
     #donors=1:max.donors;cutoff=max.donors;# debugging
    return(list(donors=donors, ndonors=cutoff))
  }
  g.donors<-apply(probmass,1,f) 
  t.ndonors<-vapply(g.donors,function(x) x$ndonors, FUN.VALUE=0)
  t.donates<-matrix(vapply(g.donors,function(x) x$donors,FUN.VALUE=rep(0,max.donors)),max.donors)
  # now calculate the right shifted and left shifted locations of the matching indices 
  t.donatesl<-t.donatesr<-matrix(0L,max.donors,G[ch]) # re-sizing and 0s first
  t.donatesl[,2:G[ch]]<-matrix(unlist(tapply(2:G[ch], 2:G[ch], function(g) match(t.donates[,g], t.donates[,g-1])),use.names=F),max.donors)
  t.donatesr[,1:(G[ch]-1)]<-matrix(unlist(tapply(1:(G[ch]-1), 1:(G[ch]-1), function(g) match(t.donates[,g], t.donates[,g+1])),use.names=F),max.donors)
  t.donatesl[is.na(t.donatesl)]<-0 # match() returns NA where no match occurs, replace with 0s
  t.donatesr[is.na(t.donatesr)]<-0 # match() returns NA where no match occurs, replace with 0s
  return(list(ndonors=t.ndonors,donates=t.donates-1,donatesl=t.donatesl-1,donatesr=t.donatesr-1,switches=switches)) # -1 to convert to cpp indexing
}
