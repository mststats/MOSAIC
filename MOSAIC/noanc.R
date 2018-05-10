# need to run this otherwise so that have approx. correct parameters to start, else donates usage is very poor
# just fit noanc on a couple of targets and a couple of chromosomes
if (length(samp_chrnos)>nchrno) samp_chrnos=chrnos # use all if try to use too many
o.nchrno=nchrno;o.chrnos=chrnos;chrnos=samp_chrnos;nchrno=length(samp_chrnos);
o.NUMA=NUMA;o.NUMI=NUMI;NUMA=min(o.NUMA,subNUMA);NUMI=max(NUMA/2,1)
o.NUMP<-NUMP;o.label<-label;o.KNOWN<-KNOWN;o.NL<-NL;o.NN<-NN
dons<-NULL;for (k in 1:kLL) dons<-c(dons,sort(sample(which(label==k),min(subNL,sum(label==k)))))
NUMP<-length(dons);label<-c(label[dons],label[!KNOWN]);NL<-c(table(label));NN<-sum(NL);KNOWN<-c(KNOWN[dons],KNOWN[!KNOWN])
# if required, use subset of targets and subset of donors on subset of chromosomes
if (nchrno!=o.nchrno | NUMA!=o.NUMA | NUMP!=o.NUMP) 
{
  o.umatch=umatch
  o.d.w=d.w
  o.t.w=t.w
  o.maxmatchsize=maxmatchsize
  o.G=G;o.flips=flips;o.gobs=gobs;
  G=G[match(samp_chrnos, o.chrnos)]; flips=gobs=list()
  for (ch in 1:nchrno)
  {
    umatch[[ch]]=umatch[[which(samp_chrnos[ch]==o.chrnos)]]
    maxmatchsize[ch]=o.maxmatchsize[which(samp_chrnos[ch]==o.chrnos)]
    d.w[[ch]]=o.d.w[[which(samp_chrnos[ch]==o.chrnos)]]
    t.w[[ch]]=o.t.w[[which(samp_chrnos[ch]==o.chrnos)]]
    if (NUMP!=o.NUMP)
      for (g in 1:G[ch])
	d.w[[ch]]$w[[g]]=d.w[[ch]]$w[[g]][dons] # subset of the donors to fit the no latent ancestry model parameters 
    if (NUMA!=o.NUMA)
      for (g in 1:G[ch])
	t.w[[ch]]$w[[g]]=t.w[[ch]]$w[[g]][1:NUMA] # subset of the targets to fit the no latent ancestry model parameters 
  }
  for (ind in 1:NUMI) {flips[[ind]]=list(); for (ch in 1:nchrno) flips[[ind]][[ch]]=o.flips[[ind]][[which(samp_chrnos[ch]==o.chrnos)]]}
  for (ch in 1:nchrno) {gobs[[ch]]=list(); for (ind in 1:NUMI) gobs[[ch]][[ind]]=o.gobs[[which(samp_chrnos[ch]==o.chrnos)]][[ind]]}
}
###############################  fit no anc model and EM a little #######################
o.L<-L;L<-1
noanc.rho=rho<-mean(rho)
noanc.theta=theta<-mean(theta)
# first compute the no-ancestry equivalent parameters Mu, rho, and theta. One for each ind.
ind.Mu=list()
for (ind in 1:NUMI) 
{
  # use current ind specific parameters from ancestry aware model
  ind.Mu[[ind]]=matrix(rowSums(t(t(Mu)*alpha[[ind]])),kLL) # p(g) = sum_a(p(g|a)p(a))
}
noanc.Mu=Mu=matrix(rowSums(Mu%*%Reduce("+",alpha)/NUMI),kLL) 
o.doMu<-doMu;doMu<-T;o.dotheta<-dotheta;dotheta<-T;o.dorho<-dorho;dorho<-T;
o.doPI<-doPI;doPI<-F;o.PI<-PI;PI=list();for (ind in 1:NUMI) PI[[ind]]<-matrix(0,1,1)
o.alpha=alpha;alpha=list();for (ind in 1:NUMI) alpha[[ind]]=1
o.lambda=lambda;lambda=list();for (ind in 1:NUMI) lambda[[ind]]=0
o.prop.don<-prop.don;o.max.donors<-max.donors
prop.don<-1;max.donors<-NUMP # use all donor haplotypes here 
for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,PI[[ind]],ind.Mu[[ind]],rho,NL)
mutmat<-fmutmat(theta, L, maxmiss, maxmatch) # possibly overkill / some redundancy as maxmiss and maxmatch may have fallen for this subset
LOG=F;source("all_donates.R") # dummy run; this will return all donors at all gridpoints and is not affected by parameter values
cloglike=NaN
LOG=o.LOG
initProb=initprobs(T,NUMA,L,NUMP,kLL,PI,Mu,rho,alpha,label,NL)

if(verbose) 
  cat("Fitting no-ancestry model\n") 
if (LOG) 
{
  tmp=create_logfile(resultsdir,target,kLL,L,NUMI,firstind,chrnos,nchrno,NN,GpcM)
  runtime=old.runtime=tmp$rtime;diff.time=0;len=tmp$len
  noancEMlogfile=EMlogfile=tmp$logfile
}
total=50 # only estimating some of the parameters, not required to be super accurate
#stop("wait")
if (EM) source("mosaic.R") # no anc fit and all donors included; should remove EM output
#stop("finished noanc part")
if (!exists("getnoancgfbs")) getnoancgfbs=F 
if (getnoancgfbs)
  noanc_gfbs=get_gfbs()
L<-o.L
# return parameters, etc to correct sizes
doMu<-o.doMu;dorho=o.dorho;dotheta=o.dotheta
doPI<-o.doPI;PI<-o.PI;alpha=o.alpha;lambda=o.lambda
noanc.rho=rho;rho<-rep(noanc.rho,L) # note that this includes all the latent ancestry switches
noanc.Mu<-Mu;Mu<-NULL; for (l in 1:L) Mu<-cbind(Mu,noanc.Mu)
noanc.theta=theta;theta<-rep(noanc.theta,L)
# next line gets called if some groups dropped but it's fast so potential redundancy is ok
for (ind in 1:NUMI) transitions[[ind]]<-s_trans(L,kLL,PI[[ind]],Mu,rho,NL)
mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
###################### return to original size of problem ###############################
if (nchrno!=o.nchrno | NUMA!=o.NUMA | NN!=o.NN)
{
  NUMP<-o.NUMP;label<-o.label;KNOWN<-o.KNOWN;NN<-o.NN;NL<-o.NL
  maxmatchsize=o.maxmatchsize
  for(ch in 1:nchrno)
  {
    umatch[[ch]]=o.umatch[[ch]]
    d.w[[ch]]=o.d.w[[ch]]
    t.w[[ch]]=o.t.w[[ch]]
  }
  rm(o.umatch,o.d.w,o.t.w)
  G=o.G;flips=o.flips;gobs=o.gobs;nchrno=o.nchrno;chrnos=o.chrnos;NUMA=o.NUMA;NUMI=o.NUMI # return to full size
  rm(o.flips,o.gobs,o.G)
}
prop.don<-o.prop.don;max.donors<-o.max.donors

