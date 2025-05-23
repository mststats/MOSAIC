# functions to estimate Fst between genomes and also Rst (see MOSAIC paper)
# use Equation 6 of Bhatia et al Genome Research (2013), but sum numerator and denominator over loci
r_wc_ests=function(p,n) # best comparison b/w Mu and 1/Fst values
{
  sumn=sum(n) #-1
  a=prod(n)/(sumn)
  bcc=sum(n*p*(1-p))/(sumn-2)#+1)
  numer=2*a*bcc
  denom=a*((p[1]-p[2])^2)+(2*a-1)*bcc
  #denom=a*((p[1]-p[2])^2)+(2*a)*bcc
  return(c(numer,denom))
}

r_calc_freqs=function(s,t.A,t.populations,t.y) 
{
  n=p=rep(0,t.A)
  for (l in 1:t.A) 
  {
    tmppops=which(t.populations[,s]==l)
    n[l]=sum(!is.nan(tmppops))
    p[l]=mean(t.y[s,tmppops],na.rm=TRUE)
  }
  return(list(p,n)) # returns two vectors of length t.A; freq and number
}

wc_ests=cmpfun(r_wc_ests,list(optimize=3))
calc_freqs=cmpfun(r_calc_freqs,list(optimize=3))
which.max.thresh=function(x,thresh) {ans=which.max(x); if (x[ans]<thresh) ans=NaN; ans}
r_maximal_alleles=function(t.target,chrnos,glocs,t.localanc,pathin1,pathin2,thresh=0.8) # assign each hap locally to an anc and return "ancestral" haps
{
  t.A=dim(t.localanc[[1]])[1]
  G=sapply(t.localanc,function(x) dim(x)[3])
  NUMA=dim(t.localanc[[1]])[2]
  allp=alln=list()
  for (l in 1:t.A)
    allp[[l]]=alln[[l]]=list()
  for (ch in 1:length(chrnos))
  {
    snps=read.table(file.path(pathin1,paste0("snpfile.",chrnos[ch])))
    S=nrow(snps)
    tmpfilename=file.path(pathin2,paste0(t.target,"genofile.",chrnos[ch]))
    tmp<-scan(tmpfilename,what="character",quiet=TRUE,nlines=1)
    N2<-nchar(tmp)
    y=matrix(suppressWarnings(as.integer(as.matrix(laf_open_fwf(tmpfilename, column_widths=rep(1,N2),column_types=rep("character",N2))[,]))),ncol=N2)
    populations=matrix(NaN,NUMA,S)
    tmp=grid_to_pos_chr(t.localanc[[ch]],snps[,4],glocs[[ch]])
    k=0
    for (ind in 1:(NUMA/2))
      for (h in 1:2)
      {
	k=(ind-1)*2+h
	populations[k,]=apply(tmp[,k,],2,which.max.thresh,thresh)
      }
    tmp=lapply(1:S,calc_freqs,t.A,populations,y)
    freqs_mat=matrix(sapply(tmp,function(x) x[[1]]),t.A) 
    for (l in 1:t.A) allp[[l]][[ch]]=freqs_mat[l,] # allele freqs for ancs (rows) along chromosome (cols)
    counts_mat=matrix(sapply(tmp,function(x) x[[2]]),t.A) 
    for (l in 1:t.A) alln[[l]][[ch]]=counts_mat[l,] # counts for ancs (rows) along chromosome (cols)
  }
  return(list("freqs"=allp,"counts"=alln))
}
maximal_alleles=cmpfun(r_maximal_alleles,list(optimize=3))

r_wc_fst=function(freqs1,counts1,freqs2,counts2) 
{
  Hs=Ht=NULL
  nchrno=length(freqs1)
  if ((length(counts1[[1]])>1) & (length(counts2[[1]])>1))
    f_wc=function(s) wc_ests(c(freqs1[[ch]][s],freqs2[[ch]][s]),c(counts1[[ch]][s],counts2[[ch]][s]))
  if ((length(counts1[[1]])>1) & (length(counts2[[1]])==1))
    f_wc=function(s) wc_ests(c(freqs1[[ch]][s],freqs2[[ch]][s]),c(counts1[[ch]][s],counts2[[ch]]))
  if ((length(counts1[[1]])==1) & (length(counts2[[1]])>1))
    f_wc=function(s) wc_ests(c(freqs1[[ch]][s],freqs2[[ch]][s]),c(counts1[[ch]],counts2[[ch]][s]))
  if ((length(counts1[[1]])==1) & (length(counts2[[1]])==1))
    f_wc=function(s) wc_ests(c(freqs1[[ch]][s],freqs2[[ch]][s]),c(counts1[[ch]],counts2[[ch]]))
  for (ch in 1:nchrno)
  {
    S=length(freqs1[[ch]])
    tmp=sapply(1:S,f_wc)
    ch_Hs=tmp[1,]
    ch_Ht=tmp[2,]
    Hs=c(Hs,ch_Hs)
    Ht=c(Ht,ch_Ht)
  }
  return((1-sum(Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE)))
}

v_wc_fst=function(freqs1,counts1,freqs2,counts2) 
{
  Hs=Ht=NULL
  nchrno=length(freqs1)
  for (ch in 1:nchrno)
  {
    sumn=counts1[[ch]]+counts2[[ch]]
    a=(counts1[[ch]]*counts2[[ch]])/(sumn)
    bcc=(counts1[[ch]]*freqs1[[ch]]*(1-freqs1[[ch]])+counts2[[ch]]*freqs2[[ch]]*(1-freqs2[[ch]]))/(sumn-2) 
    numer=2*a*bcc
    denom=a*((freqs1[[ch]]-freqs2[[ch]])^2)+(2*a-1)*bcc
    Hs=c(Hs,numer)
    Ht=c(Ht,denom)
  }
  return((1-sum(Hs,na.rm=TRUE)/sum(Ht,na.rm=TRUE)))
}

#wc_fst=cmpfun(r_wc_fst,list(optimize=3))
wc_fst=cmpfun(v_wc_fst,list(optimize=3)) # gives the same as the above

R_Fst=function(x) 
{ # difference quotient of Fst i.e. (p-q)^2/(0.5*(p+q))
  combn_aa=utils::combn(ncol(x),2)
  laa=ncol(combn_aa)
  d_fst=rep(NaN,laa)
  for (aa in 1:laa) {a1=combn_aa[1,aa]
  a2=combn_aa[2,aa]
  d_fst[aa]=mean(((x[,a1]-x[,a2])^2)/(0.5*(x[,a1]+x[,a2])),na.rm=TRUE)} # squared diff over average
  return(d_fst)
}

# function to calculate Fst between all pairs of latent ancestries and between each latent ancestry and each donor panel
Fst_combos=function(target, A, NN, panels, pathin="MOSAIC_RESULTS/FREQS") {
  pdata=ancestral_freqs=NULL # placeholder, overwritten by next line
  load(file=file.path(pathin, paste0(target, "_", A, "way_", NN, "_freqs.rdata"))) # use pre-calculated freq / count pairs from running write_admixed_summary
  anc_fst=rep(NaN,choose(A,2))
  Ls=utils::combn(A,2)
  for (l in 1:ncol(Ls))
  {
    l1=Ls[1,l];l2=Ls[2,l]
    anc_fst[l]=wc_fst(ancestral_freqs$freqs[[l1]],ancestral_freqs$counts[[l1]],ancestral_freqs$freqs[[l2]],ancestral_freqs$counts[[l2]]) 
    names(anc_fst)[l]=paste0("anc",l1,"x","anc",l2)
  }
  panels_fst=matrix(NaN,length(panels),A)
  for (l1 in 1:length(panels))
  {
    load(file.path(pathin, paste0(panels[l1],"_freqs.rdata")))
    for (l2 in 1:A) 
      panels_fst[l1,l2]=wc_fst(ancestral_freqs$freqs[[l2]],ancestral_freqs$counts[[l2]],pdata$freqs,pdata$counts)
  }
  rownames(panels_fst)=panels
  colnames(panels_fst)=paste("anc",1:A)
  Rst=R_Fst(panels_fst);names(Rst)=names(anc_fst)
  return(list("ancs"=anc_fst, "Rst"=Rst, "panels"=panels_fst))
}

Fst_panels=function(panel1,panel2, pathin="MOSAIC_RESULTS/FREQS") {
  load(file.path(pathin, paste0(panel1,"_freqs.rdata")))
  tmp1=pdata
  load(file.path(pathin, paste0(panel2,"_freqs.rdata")))
  tmp2=pdata
  return(wc_fst(tmp1$freqs,tmp1$counts,tmp2$freqs,tmp2$counts))
}

