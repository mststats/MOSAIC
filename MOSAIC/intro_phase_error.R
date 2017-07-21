intro_phase_error<-function(ch, ind, t.flips, verbose) {
  hap<-c(ind*2-1,ind*2) # ind indexes over genotypes, hap the two haplotypes
  het.g<-rep(T,G[ch]) # use all gridpoints; change to phase effective everywhere
  het.g[1:3]<-F # het gridpoints not 1 or G[ch]
  if (!exists("g.w")) g.w=1e2
  # take two admixed haps and introduce NE phase flip errors 
  NE=as.integer(RPE*sum(het.g)/100) # %RPE phase flip errors
  NE=max(1,NE) # do at least one
  ind.phase.error.locs<-NULL
  if (verbose) cat("introducing ", NE*2, "phase errors on chromosome", chrnos[ch], "of ind", ind, "\n")
    if (NE>0)
  for (ne in 1:NE) {
    # make t.flips that include a het site
    g<-sample(1:G[ch], 1, prob=het.g*1) # this is the het site. g.l and g.u will be chosen around it w/in g.w either side
    #ind.phase.error.locs<-c(ind.phase.error.locs,g)
    #t.flips[g:G[ch]]<-!t.flips[g:G[ch]] # could have been flipped before
    g.l<-sample(max(1,g-g.w):g,1)
    g.u<-sample((g+1):min(G[ch],g+g.w),1)
    ind.phase.error.locs<-c(ind.phase.error.locs,g.l,g.u)
    t.flips[g.l:G[ch]]<-!t.flips[g.l:G[ch]] # could have been flipped before
    t.flips[g.u:G[ch]]<-!t.flips[g.u:G[ch]] # could have been flipped before
    }
  return(list(ind.phase.error.locs=ind.phase.error.locs, err.flips=t.flips))
}
