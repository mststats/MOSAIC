# various contrived special cases to test how the model reacts

# create an extra latent ancestry
if (!exists("extra_ancs")) extra_ancs=F
if (extra_ancs)
{
  o.L=L;L=L+extra_ancs
  for (ind in 1:NUMI) sim.alpha[[ind]]=c(sim.alpha[[ind]],rep(0,extra_ancs))
  tmpg=g.true_anc;for (ch in 1:nchrno) {g.true_anc[[ch]]=array(0,c(L,NUMA,G[ch]));g.true_anc[[ch]][1:o.L,,]=tmpg[[ch]]};rm(tmpg)
}

if (!exists("singleind")) singleind=F # this is purely for sanity testing the code
if (singleind) # stitch all genomes together as a single individual
{
  if (HPC) stop("can't do this under HPC mode without hassle")
  single.gobs=single.MATCH=single.g.true_anc=single.g.loc=list()
  for (ch in 1:nchrno)
  {
    for (ind in 1:NUMI) 
    {
      single.gobs[[(ch-1)*NUMI+ind]]=list(gobs[[ch]][[ind]])
      single.MATCH[[(ch-1)*NUMI+ind]]=list(MATCH[[ch]][[ind]])
      hap=c((ind-1)*2+1,ind*2)
      single.g.true_anc[[(ch-1)*NUMI+ind]]=g.true_anc[[ch]][,hap,]
      single.g.loc[[(ch-1)*NUMI+ind]]=g.loc[[ch]]
    }
  }
  gobs=single.gobs
  MATCH=single.MATCH
  g.true_anc=single.g.true_anc
  g.loc=single.g.loc
  G=c(t(matrix(rep(G,NUMI),nchrno)))  
  sim.alpha=list(Reduce("+",sim.alpha)/NUMI)
  nchrno=NUMI*nchrno
  chrnos=1:nchrno
  NUMI=1;NUMA=2 
}

if (singleQ)
{
  tmpalpha=Reduce("+",sim.alpha)/NUMI
  for (ind in 1:NUMI)
    sim.alpha[[ind]]=tmpalpha
}
