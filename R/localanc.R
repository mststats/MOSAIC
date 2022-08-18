# function to calculate the local ancestry estimates based on temporarily constructed parts of the gridded forward-backward probabilities 
get_localanc=function(t.NUMP, t.nchrno, t.max.donors, t.donates, t.donatesl, t.donatesr, t.transitions, t.umatch, t.maxmatchsize, t.d.w, t.t.w, 
		      t.gobs, t.mutmat, t.maxmiss, t.initProb, t.label, t.ndonors, t.flips, t.HPC,  
		      t.G,t.A,t.kLL,t.NUMA,t.NUMI,tol=1e-8,t.g.true_anc=NULL) {
  localanc<-list()
  THIN=ifelse(t.max.donors==t.NUMP, FALSE, TRUE)
  t.NUMI=t.NUMA/2
  for (ch in 1:t.nchrno)
  {
    localanc[[ch]]<-array(0,c(t.A,t.NUMA,t.G[ch]))
    if (t.HPC==1)
    {
      donates_chr=getdonates(t.donates[[ch]],t.NUMI)
      donatesl_chr=getdonates(t.donatesl[[ch]],t.NUMI)
      donatesr_chr=getdonates(t.donatesr[[ch]],t.NUMI)
      tmplocalanc<-foreach(k=1:t.NUMA) %dopar% 
      {
	ind=as.integer((k+1)/2)
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.A);t.sumfors<-matrix(0,t.G[ch],t.A);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesl_chr[[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	t.backs<-rep(0,t.G[ch]*t.max.donors*t.A);t.scalefactorb<-rep(0,t.G[ch]);
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		    t.mutmat,t.maxmiss,t.label,t.ndonors[[ch]][[ind]],donates_chr[[ind]],donatesr_chr[[ind]],t.flips[[ind]][[ch]],t.backs,t.scalefactorb)
	tmp2=cppgforback(t.max.donors,THIN,t.kLL,t.NUMP,t.label,t.A,t.G[ch],t.ndonors[[ch]][[ind]],donates_chr[[ind]],t.fors,t.backs)
	ans=matrix(0,t.A,t.G[ch])
	for (j in 1:t.A) 
	  for (jk in 1:t.kLL) 
	    ans[j,]=ans[j,]+tmp2[,(j-1)*t.kLL+jk]
	ans
      }
    }
    if (t.HPC==2)
    {
      tmplocalanc<-foreach(k=1:t.NUMA) %dopar% 
      {
	ind=as.integer((k+1)/2)
	donates_chr_ind=getdonates_ind(t.donates[[ch]][[ind]])
	donatesl_chr_ind=getdonates_ind(t.donatesl[[ch]][[ind]])
	donatesr_chr_ind=getdonates_ind(t.donatesr[[ch]][[ind]])
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.A);t.sumfors<-matrix(0,t.G[ch],t.A);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesl_chr_ind,t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	t.backs<-rep(0,t.G[ch]*t.max.donors*t.A);t.scalefactorb<-rep(0,t.G[ch]);
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		    t.mutmat,t.maxmiss,t.label,t.ndonors[[ch]][[ind]],donates_chr_ind,donatesr_chr_ind,t.flips[[ind]][[ch]],t.backs,t.scalefactorb)
	tmp2=cppgforback(t.max.donors,THIN,t.kLL,t.NUMP,t.label,t.A,t.G[ch],t.ndonors[[ch]][[ind]],donates_chr_ind,t.fors,t.backs)
	rm(donates_chr_ind,donatesl_chr_ind,donatesr_chr_ind)
	ans=matrix(0,t.A,t.G[ch])
	for (j in 1:t.A) 
	  for (jk in 1:t.kLL) 
	    ans[j,]=ans[j,]+tmp2[,(j-1)*t.kLL+jk]
	ans
      }
    }
    if (!t.HPC)
    {
      tmplocalanc<-foreach(k=1:t.NUMA) %dopar% 
      {
	ind=as.integer((k+1)/2)
	t.fors<-rep(0,t.G[ch]*t.max.donors*t.A);t.sumfors<-matrix(0,t.G[ch],t.A);t.scalefactor<-rep(0,t.G[ch]);
	cppforward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.kLL,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.t.w[[ch]],t.gobs[[ch]][[ind]],
		   t.mutmat,t.maxmiss,t.initProb[k,],t.label,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesl[[ch]][[ind]],t.flips[[ind]][[ch]],t.fors,t.sumfors,t.scalefactor)
	t.backs<-rep(0,t.G[ch]*t.max.donors*t.A);t.scalefactorb<-rep(0,t.G[ch]);
	cppbackward(k,t.NUMA,t.max.donors,THIN,t.NUMP,t.A,0,t.G[ch],t.G[ch],t.transitions[[ind]],t.umatch[[ch]],t.maxmatchsize[ch],t.d.w[[ch]],t.w[[ch]],t.gobs[[ch]][[ind]],
		    t.mutmat,t.maxmiss,t.label,t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.donatesr[[ch]][[ind]],t.flips[[ind]][[ch]],t.backs,t.scalefactorb)
	tmp2=cppgforback(t.max.donors,THIN,t.kLL,t.NUMP,t.label,t.A,t.G[ch],t.ndonors[[ch]][[ind]],t.donates[[ch]][[ind]],t.fors,t.backs)
	ans=matrix(0,t.A,t.G[ch])
	for (j in 1:t.A) 
	  for (jk in 1:t.kLL) 
	    ans[j,]=ans[j,]+tmp2[,(j-1)*t.kLL+jk]
	ans
      }
    }
    for (k in 1:t.NUMA)
    {
      localanc[[ch]][,k,]=tmplocalanc[[k]]
      #localanc[[ch]][,k,]<-round(localanc[[ch]][,k,],3) # Hapmix does this!
      localanc[[ch]][,k,][localanc[[ch]][,k,]<tol]<-tol;localanc[[ch]][,k,][localanc[[ch]][,k,]>(1-tol)]<-1-tol
      localanc[[ch]][,k,]<-t(t(localanc[[ch]][,k,])/apply(localanc[[ch]][,k,],2,sum))
    }
  }

  lans=list(localanc=localanc)
  # re-order truth to be same labels
  if (!is.null(t.g.true_anc))
  {
    r.ord<-function(ord) 
    {
      tmp.anc<-t.g.true_anc
      for (ch in 1:t.nchrno) tmp.anc[[ch]][,,]=t.g.true_anc[[ch]][ord,,]
      return(dip_fr(tmp.anc,localanc))
    }
    all.ord<-permn(t.A)
    best.r<-0
    for (i in 1:length(all.ord))
      if (r.ord(all.ord[[i]])>best.r) 
      {
	best.r=r.ord(all.ord[[i]])
	best.ord<-all.ord[[i]]
      }
    for (ch in 1:t.nchrno) 
      t.g.true_anc[[ch]][,,]<-t.g.true_anc[[ch]][best.ord,,]
    lans$g.true_anc=t.g.true_anc
  }
  return(lans)
}

# function to calculate the local ancestry estimates based on the gridded forward-backward probabilities gfbs
get_localanc_using_gfbs=function(t.NUMP, t.nchrno, t.max.donors, t.donates, t.donatesl, t.donatesr, t.transitions, t.umatch, t.maxmatchsize, t.d.w, t.t.w, 
		      t.gobs, t.mutmat, t.maxmiss, t.initProb, t.label, t.ndonors, t.flips, t.HPC,  
		      t.G,t.A,t.kLL,t.NUMA,t.NUMI,tol=1e-8,t.g.true_anc=NULL) {
  t.gfbs=get_gfbs(t.NUMP, t.nchrno, t.max.donors, t.donates, t.donatesl, t.donatesr, t.NUMA, t.A, t.G, t.kLL, t.transitions, t.umatch, t.maxmatchsize, 
		  t.d.w, t.t.w, t.gobs, t.mutmat, t.maxmiss, t.initProb, t.label, t.ndonors, t.flips, t.HPC)
  localanc<-list()
  for (ch in 1:t.nchrno)
  {
    localanc[[ch]]<-array(1/t.A,c(t.A,t.NUMA,t.G[ch]))
    for (k in 1:t.NUMA)
    {
      localanc[[ch]][,k,]<-0
      for (j in 1:t.A) 
	for (jk in 1:t.kLL) 
	  localanc[[ch]][j,k,]<-localanc[[ch]][j,k,] + t.gfbs[[ch]][[k]][,(j-1)*t.kLL+jk]
      #localanc[[ch]][,k,]<-round(localanc[[ch]][,k,],3) # Hapmix does this!
      localanc[[ch]][,k,][localanc[[ch]][,k,]<tol]<-tol;localanc[[ch]][,k,][localanc[[ch]][,k,]>(1-tol)]<-1-tol
      localanc[[ch]][,k,]<-t(t(localanc[[ch]][,k,])/apply(localanc[[ch]][,k,],2,sum))
    }
  }
  ans=list(localanc=localanc)
  # re-order truth to be same labels
  if (!is.null(t.g.true_anc))
  {
    r.ord<-function(ord) 
    {
      tmp.anc<-t.g.true_anc
      for (ch in 1:t.nchrno) tmp.anc[[ch]][,,]=t.g.true_anc[[ch]][ord,,]
      return(dip_fr(tmp.anc,localanc))
    }
    all.ord<-permn(t.A)
    best.r<-0
    for (i in 1:length(all.ord))
      if (r.ord(all.ord[[i]])>best.r) 
      {
	best.r=r.ord(all.ord[[i]])
	best.ord<-all.ord[[i]]
      }
    for (ch in 1:t.nchrno) 
      t.g.true_anc[[ch]][,,]<-t.g.true_anc[[ch]][best.ord,,]
    ans$g.true_anc=t.g.true_anc
  }
  return(ans)
}

# function to find snp positions values of x based on gridded values of x
grid_to_pos_chr=function(x,pos,glocs) { # arguments are thing-to-map, SNP positions, grid locations
  S=length(pos)
  g.map<-vapply(1:S, function(s) which.min((pos[s]-glocs)^2),0L) # create map from rates to grid
  if (length(dim(x))==2)
    ans=x[,g.map]
  if (length(dim(x))==3)
    ans=x[,,g.map]
  return(ans)
}

grid_to_pos=function(x,pathin,glocs,chrnos) { # arguments are thing-to-map, SNP positions, grid locations
  ans=list()
  for (ch in 1:length(chrnos)) {
    pos=read.table(paste0(pathin,"snpfile.",chrnos[ch]))[,4]
    S=length(pos)
    g.map<-vapply(1:S, function(s) which.min((pos[s]-glocs[[ch]])^2),0L) # create map from rates to grid
    if (length(dim(x[[ch]]))==2)
      ans[[ch]]=x[[ch]][,g.map]
    if (length(dim(x[[ch]]))==3)
      ans[[ch]]=x[[ch]][,,g.map]
  }
  return(ans)
}

# calculate proportions from local ancestry information
alpha_from_local=function(y) {ans=apply(sapply(y, function(x) apply(x,1,sum)),1,sum);return(ans/sum(ans))}
alpha_from_local_ind=function(y,ind) {haps=c(ind*2-1,ind*2);ans=apply(sapply(y, function(x) apply(x[,haps,],1,sum)),1,sum);return(ans/sum(ans))}
