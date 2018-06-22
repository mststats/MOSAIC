# function to compute mapping from markers to evenly spaced (in genetic distance) gridpoints
# t.S_chr is number of observed loci, t.G_chr is number of gridpoints; note this is called for a particular chromosome
create_grid=function(t.G_chr,t.S_chr,t.g.map, t.chrno, t.NUMA, t.L, t.umatch_chr, t.t.w_chr, t.true_anc_chr, t.locs, t.NUMI){
  if (verbose) cat("mapping chr", t.chrno, "to a grid...\n")
  # need to evenly spread total gridpoints across chromosomes s.t. each gap is the same #morgans
  # create a map of observed loci to gridded loci; grid is even on t.rates not distances
  if (verbose) cat("Finding new positions on chr", t.chrno, "...\n")
  # the above is lazy. Should use all_rates rather than thinned to SNPs rates. 
  if (verbose) cat("Finding number at each location on chr", t.chrno, "...\n")
  #g.loc_chr gives physical locus of each gridpoint; average if more than 1 or more obs; else average of nearest two obs
  g.loc_chr<-rep(0L,t.G_chr); for (s in 1:t.S_chr) g.loc_chr[t.g.map[s]]<-g.loc_chr[t.g.map[s]]+t.locs[s]/sum(t.g.map==t.g.map[s]) 
  if (target=="simulated")
  {
    if (verbose) cat("Mapping true ancestry array for chr", t.chrno, "to the grid\n")
    g.true_anc_chr<-array(0,c(t.L,t.NUMA,t.G_chr))
    for (i in 1:t.L)
      for (k in 1:t.NUMA)
	for (s in 1:t.S_chr)
	  g.true_anc_chr[i,k,t.g.map[s]]<-g.true_anc_chr[i,k,t.g.map[s]]+t.true_anc_chr[i,k,s]/sum(t.g.map==t.g.map[s])
  }
  emptyg<-which(g.loc_chr==0)
  if (length(emptyg)>0)
  {
    if (emptyg[1]==1) # fix first if empty
    {
      tmpg<-min(which(g.loc_chr!=0 & (1:t.G_chr)>1)) # nearest one up
      if (target=="simulated") g.true_anc_chr[i,k,1]<-g.true_anc_chr[i,k,tmpg]
      g.loc_chr[1]<-g.loc_chr[tmpg]
      emptyg<-emptyg[-1] # mark as none empty
    }
    maxg<-length(emptyg)
    if (emptyg[maxg]==t.G_chr) # fix last if empty
    {
      tmpg<-max(which(g.loc_chr!=0 & (1:t.G_chr)<t.G_chr)) # nearest one down
      if (target=="simulated") g.true_anc_chr[i,k,t.G_chr]<-g.true_anc_chr[i,k,tmpg]
      g.loc_chr[t.G_chr]<-g.loc_chr[tmpg]
      emptyg<-emptyg[-maxg] # mark as non-empty
    }
    # now fill in other empty gridpoints
    for (eg in 1:length(emptyg))
    {
      g=emptyg[eg]
      a<-max(which(g.loc_chr!=0 & (1:t.G_chr)<g))
      b<-min(which(g.loc_chr!=0 & (1:t.G_chr)>g))
      wa=1/(g-a)
      wb=1/(b-g)
      ws=wa+wb
      wa=wa/ws
      wb=wb/ws
      #a=b # use to the right as per admix.R 
      #b=a # use to the right as per admix.R 
      if (target=="simulated") 
	for (i in 1:t.L)
	  for (k in 1:t.NUMA)
	    g.true_anc_chr[i,k,g]<-wa*g.true_anc_chr[i,k,a]+wb*g.true_anc_chr[i,k,b] # linear interpolation; sample jump would be more realistic, but no data to see locus
      g.loc_chr[g]<-wa*g.loc_chr[a]+wb*g.loc_chr[b] # weighted average of physical locus of the two nearest gridpoint w/ data
    }
  }
  rm(emptyg)
  gobs_chr<-list()
  for (ind in 1:t.NUMI) 
  {
    haps=c(ind*2-1,ind*2)
    gobs_chr[[ind]]<-sapply(1:t.G_chr,function(g) max(t.umatch_chr[[g]][,t.t.w_chr$w[[g]][haps]+1]))
    gobs_chr[[ind]][is.na(gobs_chr[[ind]])]<-0
    # return maxmatch and maxmiss for this chromosome only and find max across all afterwards
    maxmatch_chr<-max(unlist(t.umatch_chr))
    maxmiss_chr<-max(gobs_chr[[ind]]-sapply(1:t.G_chr, function(g) min(t.umatch_chr[[g]][,t.t.w_chr$w[[g]][haps]+1])),na.rm=T)
  }
  if (target=="simulated") return(list(g.loc_chr=g.loc_chr,gobs_chr=gobs_chr,maxmatch_chr=maxmatch_chr,maxmiss_chr=maxmiss_chr,
				       g.true_anc_chr=g.true_anc_chr))
  if (target!="simulated") return(list(g.loc_chr=g.loc_chr,gobs_chr=gobs_chr,maxmatch_chr=maxmatch_chr,maxmiss_chr=maxmiss_chr))
}
