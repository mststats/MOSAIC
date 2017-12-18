HapMixpath="../../HapmixReleasev2/"
HapMixrunpath=paste0(HapMixpath,"RAINBOW/RUN/")
source("calc_r2.R")
source("plot_funcs.R")
if (!exists("showunphased")) showunphased=F
hap_hmixlocalanc=hmixlocalanc=list()
hmixloc=list()
gs=list()
oldwd=getwd()
setwd(HapMixpath)

if (HapMixMODE=="HAP") 
{
  LIM=NUMA
  NCHR=1;MODE="PROB"
}
if (HapMixMODE=="DIP" | HapMixMODE=="BOTH") 
{
  LIM=NUMA/2
  NCHR=2;MODE="LOCALANC";
}

for (chrno in chrnos) source("read_results.R")
setwd(oldwd)
for (ch in 1:nchrno)
{
  if (HapMixMODE=="DIP" | HapMixMODE=="BOTH") 
    load(paste0(HapMixrunpath, "HGDP.DIP_hapmix_anc_",chrnos[ch],".rdata"))
  if (HapMixMODE=="HAP") 
    load(paste0(HapMixrunpath, "HGDP.HAP_hapmix_anc_",chrnos[ch],".rdata"))
  hmixlocalanc[[ch]]=hlocalanc
  # hap_hmixlocalanc only used after conversion to diploid; convenient format for r2 calcs
  if (HapMixMODE=="DIP" | HapMixMODE=="BOTH") 
  {
    hap_hmixlocalanc[[ch]]=array(0,c(dim(hlocalanc)[1],dim(hlocalanc)[2]*2,dim(hlocalanc)[3]))
    hap_hmixlocalanc[[ch]][,seq(1,NUMA,2),]=hap_hmixlocalanc[[ch]][,seq(2,NUMA,2),]=hlocalanc/2 # split prob mass evenly across haps
  }
  if (HapMixMODE=="HAP")
  {
    hap_hmixlocalanc[[ch]]=hlocalanc 
  }
  hmixloc[[ch]]=hloc
  #gs[[ch]]=1:G[ch] # all gridpoints
  #gs[[ch]]=which(apply(apply((g.true_anc[[ch]]==0 | g.true_anc[[ch]]==1), 2:3, all),2,all)) # only look where g.true defined as 0 or 1
  gs[[ch]]=which(gobs[[ch]][[1]]>0) # only look where there are SNPs assigned; assume all have obs at same gridpoints
  # only look where there are SNPs assigned and truth is 0 or 1
  #tmp=match(which(apply(apply((g.true_anc[[ch]]==0 | g.true_anc[[ch]]==1), 2:3, all),2,all)),gs[[ch]])
  #tmp=tmp[!is.na(tmp)]
  #gs[[ch]]=gs[[ch]][tmp]
  gs[[ch]]=unique(gs[[ch]])
}
# subset to look at gridpoints where we know truth only
g.true_anc_gs=localanc_gs=noanc_unphased_localanc_gs=list()
for (ch in 1:nchrno)
{
  g.true_anc_gs[[ch]]=g.true_anc[[ch]][,,gs[[ch]]]
  localanc_gs[[ch]]=localanc[[ch]][,,gs[[ch]]]
  noanc_unphased_localanc_gs[[ch]]=noanc_unphased_localanc[[ch]][,,gs[[ch]]]
}

# sometimes HapMix needs re-r=ordering of ancestries
for (ch in 1:nchrno) 
{
  for (k in 1:NUMI)
    {
    haps=c(k*2-1,k*2)
    reorderhlocal=hap_hmixlocalanc[[ch]][2:1,haps,]
    if (sum(((hap_hmixlocalanc[[ch]][,haps[1],]+hap_hmixlocalanc[[ch]][,haps[2],])-(true_anc[[ch]][,haps[1],]+true_anc[[ch]][,haps[2],]))^2) >  
	sum(((reorderhlocal[,1,]+reorderhlocal[,2,])-(true_anc[[ch]][,haps[1],]+true_anc[[ch]][,haps[2],]))^2))
      hap_hmixlocalanc[[ch]][,haps,]=reorderhlocal
    }
}
