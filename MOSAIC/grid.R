# script to compute mapping from markers to evenly (in genetic distace) spaced gridpoints
# S is number of observed loci, G is number of gridpoints
if (verbose) cat("mapping chr", chrnos[ch], "to a grid...\n")
# need to evenly spread total gridpoints across chromosomes s.t. each gap is the same #morgans
# note that rates is chrno dependent but we are in a loop over chrno in this script
# G will be correct / consistent w/in 0.5 and is large so almost exactly the same dr across chromosomes
G[ch]<-as.integer((rates[S[ch]]-rates[1])/dr+1)
if (FLAT) {G[ch]=as.integer(S[ch]);dr=(rates[S[ch]]-rates[1])/(G[ch]-1);}
g.rates<-seq(rates[1],rates[S[ch]],l=G[ch]) # even grid across recombination rates
# create a map of observed loci to gridded loci; grid is even on rates not distances
if (verbose) cat("Finding new positions on chr", chrnos[ch], "...\n")
g.map<-vapply(1:S[ch], function(s) which.min((rates[s]-g.rates)^2),0L) # create map from rates to grid
# the above is lazy. Should use all_rates rather than thinned to SNPs rates. 
rm(g.rates)
if (verbose) cat("Finding number at each location on chr", chrnos[ch], "...\n")
#g.loc gives physical locus of each gridpoint; average if more than 1 or more obs; else average of nearest two obs
g.loc[[ch]]<-rep(0L,G[ch]); for (s in 1:S[ch]) g.loc[[ch]][g.map[s]]<-g.loc[[ch]][g.map[s]]+locs[s]/sum(g.map==g.map[s]) 
if (target=="simulated")
{
  if (verbose) cat("Mapping true ancestry array for chr", chrnos[ch], "to the grid\n")
  g.true_anc[[ch]]<-array(0,c(L,NUMA,G[ch]))
  for (i in 1:L)
    for (k in 1:NUMA)
      for (s in 1:S[ch])
	g.true_anc[[ch]][i,k,g.map[s]]<-g.true_anc[[ch]][i,k,g.map[s]]+true_anc[[ch]][i,k,s]/sum(g.map==g.map[s])
}
emptyg<-which(g.loc[[ch]]==0)
if (length(emptyg)>0)
{
  if (emptyg[1]==1) # fix first if empty
  {
    tmpg<-min(which(g.loc[[ch]]!=0 & (1:G[ch])>1)) # nearest one up
    if (target=="simulated") g.true_anc[[ch]][i,k,1]<-g.true_anc[[ch]][i,k,tmpg]
    g.loc[[ch]][1]<-g.loc[[ch]][tmpg]
    emptyg<-emptyg[-1] # mark as none empty
  }
  maxg<-length(emptyg)
  if (emptyg[maxg]==G[ch]) # fix last if empty
  {
    tmpg<-max(which(g.loc[[ch]]!=0 & (1:G[ch])<G[ch])) # nearest one down
    if (target=="simulated") g.true_anc[[ch]][i,k,G[ch]]<-g.true_anc[[ch]][i,k,tmpg]
    g.loc[[ch]][G[ch]]<-g.loc[[ch]][tmpg]
    emptyg<-emptyg[-maxg] # mark as non-empty
  }
  # now fill in other empty gridpoints
  for (eg in 1:length(emptyg))
  {
    g=emptyg[eg]
    a<-max(which(g.loc[[ch]]!=0 & (1:G[ch])<g))
    b<-min(which(g.loc[[ch]]!=0 & (1:G[ch])>g))
    wa=1/(g-a)
    wb=1/(b-g)
    ws=wa+wb
    wa=wa/ws
    wb=wb/ws
    #a=b # use to the right as per admix.R 
    #b=a # use to the right as per admix.R 
    if (target=="simulated") 
      for (i in 1:L)
	for (k in 1:NUMA)
	  g.true_anc[[ch]][i,k,g]<-wa*g.true_anc[[ch]][i,k,a]+wb*g.true_anc[[ch]][i,k,b] # linear interpolation; sample jump would be more realistic, but no data to see locus
    g.loc[[ch]][g]<-wa*g.loc[[ch]][a]+wb*g.loc[[ch]][b] # weighted average of physical locus of the two nearest gridpoint w/ data
  }
}
rm(emptyg)
gobs[[ch]]<-list()
for (ind in 1:NUMI) 
{
  haps=c(ind*2-1,ind*2)
  maxmatch<-max(maxmatch,max(unlist(umatch[[ch]])))
  gobs[[ch]][[ind]]<-sapply(1:G[ch],function(g) max(umatch[[ch]][[g]][,t.w[[ch]]$w[[g]][haps]+1]))
  gobs[[ch]][[ind]][is.na(gobs[[ch]][[ind]])]<-0
  maxmiss<-max(maxmiss,gobs[[ch]][[ind]]-sapply(1:G[ch], function(g) min(umatch[[ch]][[g]][,t.w[[ch]]$w[[g]][haps]+1])),na.rm=T)
}
