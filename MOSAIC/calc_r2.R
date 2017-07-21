dip_expected_fr2_chr_ind<-function(x,ch,ind)
{
  tmpG=dim(x[[ch]])[3]
  hap<-c(ind*2-1,ind*2)
  ans=0
  for (a in 1:L)
  {
    p0<-(1-x[[ch]][a,hap[1],])*(1-x[[ch]][a,hap[2],]) # prob hom of not anc a
    p1<-x[[ch]][a,hap[1],]*(1-x[[ch]][a,hap[2],])+(1-x[[ch]][a,hap[1],])*x[[ch]][a,hap[2],] # prob het anc a
    p2<-x[[ch]][a,hap[1],]*x[[ch]][a,hap[2],] # prob hom of anc a
    px=p1+2*p2 # expected number of a alleles out of 2
    varp=sum((px)^2)/(tmpG)-sum(px/tmpG)^2
    varp=ifelse(varp<0,0,varp) # effectively ignore if no contribution due to no or all a ancestry here
    #varx=(sum(p1*(1-p1)+4*p0*p2)/(tmpG)+varp)
    #ans=varp/varx # not good when tested 
    #acov=sum(px^2)/(tmpG-1)-sum(px)^2/tmpG/(tmpG-1)
    #avarx=sum(p1+4*p2)/(tmpG-1)-sum(px)^2/tmpG/(tmpG-1) # unstable
    avarx = (sum(p1*(1-p1)+4*p0*p2)/(tmpG)+varp)
    #acor=sqrt(varp/avarx) # 
    ar2=ifelse(varp<1e-6, 1, varp/avarx) # leave out negligible contributions
    ans=ans+ar2/L
  }
  return(ans)
}
dip_expected_fr2_ind<-function(x,ind)
{
  vecG=sapply(x,function(y) dim(y)[3])
  sumG=sum(vecG)
  ans=0
  for (a in 1:L)
  {
    vecx=matrix(NaN,2,sumG)
    OFFSET=0
    for (ch in 1:nchrno)
    {
      hap<-c(ind*2-1,ind*2)
      vecx[1,(1+OFFSET):(vecG[ch]+OFFSET)]=x[[ch]][a,hap[1],] # all first haps
      vecx[2,(1+OFFSET):(vecG[ch]+OFFSET)]=x[[ch]][a,hap[2],] # all second haps
      OFFSET=OFFSET+vecG[ch]
    }
    p0<-(1-vecx[1,])*(1-vecx[2,]) # prob hom of not anc a
    p1<-vecx[1,]*(1-vecx[2,])+(1-vecx[1,])*vecx[2,] # prob het anc a
    p2<-vecx[1,]*vecx[2,] # prob hom of anc a
    px=p1+2*p2 # expected number of a alleles out of 2
    varp=sum((px)^2)/(sumG)-sum(px/sumG)^2
    varp=ifelse(varp<0,0,varp) # effectively ignore if no contribution due to no or all a ancestry here
    avarx = (sum(p1*(1-p1)+4*p0*p2)/(sumG)+varp)
    ar2=ifelse(varp<1e-6, 1, varp/avarx) # leave out negligible contributions
    ans=ans+ar2/L
  }
  return(ans)
}
dip_expected_fr2<-function(x) 
{
  vecG=sapply(x,function(y) dim(y)[3])
  NUMI=dim(x[[1]])[2]/2
  L=dim(x[[1]])[1]
  nchrno=length(x)
  sumG=sum(vecG)*NUMI
  ans=0
  for (a in 1:L)
  {
    vecx=matrix(NaN,2,sumG)
    OFFSET=0
    for (ch in 1:nchrno)
    {
      for (ind in 1:NUMI)
      {
        hap<-c(ind*2-1,ind*2)
        vecx[1,(vecG[ch]*(ind-1)+1+OFFSET):(vecG[ch]*ind+OFFSET)]=x[[ch]][a,hap[1],] # all first haps
        vecx[2,(vecG[ch]*(ind-1)+1+OFFSET):(vecG[ch]*ind+OFFSET)]=x[[ch]][a,hap[2],] # all second haps
      }
      OFFSET=OFFSET+vecG[ch]*NUMI
    }
    p0<-(1-vecx[1,])*(1-vecx[2,]) # prob hom of not anc a
    p1<-vecx[1,]*(1-vecx[2,])+(1-vecx[1,])*vecx[2,] # prob het anc a
    p2<-vecx[1,]*vecx[2,] # prob hom of anc a
    px=p1+2*p2 # expected number of a alleles out of 2
    varp=sum((px)^2)/(sumG)-sum(px/sumG)^2
    varp=ifelse(varp<0,0,varp) # effectively ignore if no contribution due to no or all a ancestry here
    avarx = (sum(p1*(1-p1)+4*p0*p2)/(sumG)+varp)
    ar2=ifelse(varp<1e-6, 1, varp/avarx) # leave out negligible contributions
    ans=ans+ar2/L
  }
  return(ans)
}
hap_expected_fr2_chr_k<-function(x,ch,k)
{
  tmpG=dim(x[[ch]])[3]
  ans=0
  for (a in 1:L) # average over choices of a and hap(ind)
  {
    px=x[[ch]][a,k,]
    varp=sum((px)^2)/(tmpG)-sum(px/tmpG)^2
    varp=ifelse(varp<0,0,varp) # effectively ignore if no contribution due to no or all a ancestry here
    avarx=(sum(px*(1-px))/(tmpG)+varp)
    #acor=sqrt(varp/avarx) 
    ar2=ifelse(varp<1e-6, 1, varp/avarx) # leave out negligible contributions
    ans=ans+ar2
  }
  return(ans/L)
}
hap_expected_fr2_hap<-function(x,k)
{
  vecG=sapply(x,function(y) dim(y)[3])
  sumG=sum(vecG)
  ans=0
  for (a in 1:L)
  {
    vecx=rep(NaN,sumG)
    OFFSET=0
    for (ch in 1:nchrno)
    {
      vecx[(1+OFFSET):(vecG[ch]+OFFSET)]=x[[ch]][a,k,] # all first haps
      OFFSET=OFFSET+vecG[ch]
    }
    p0<-(1-vecx)*(1-vecx) # prob hom of not anc a
    p1<-vecx*(1-vecx)+(1-vecx)*vecx # prob het anc a
    p2<-vecx*vecx # prob hom of anc a
    px=p1+2*p2 # expected number of a alleles out of 2
    varp=sum((px)^2)/(sumG)-sum(px/sumG)^2
    varp=ifelse(varp<0,0,varp) # effectively ignore if no contribution due to no or all a ancestry here
    avarx = (sum(p1*(1-p1)+4*p0*p2)/(sumG)+varp)
    ar2=ifelse(varp<1e-6, 1, varp/avarx) # leave out negligible contributions
    ans=ans+ar2/L
  }
  return(ans)
}
hap_expected_fr2<-function(x) 
{
  vecG=sapply(x,function(y) dim(y)[3])
  sumG=sum(vecG)*NUMA
  ans=0
  for (a in 1:L)
  {
    vecx=rep(NaN,sumG)
    OFFSET=0
    for (ch in 1:nchrno)
    {
      for (k in 1:NUMA)
        vecx[(vecG[ch]*(k-1)+1+OFFSET):(vecG[ch]*k+OFFSET)]=x[[ch]][a,k,] 
      OFFSET=OFFSET+vecG[ch]*NUMA
    }
    px=vecx
    varp=sum((px)^2)/(sumG)-sum(px/sumG)^2
    varp=ifelse(varp<0,0,varp) # effectively ignore if no contribution due to no or all a ancestry here
    avarx=(sum(px*(1-px))/(sumG)+varp)
    ar2=ifelse(varp<1e-6, 1, varp/avarx) # leave out negligible contributions
    ans=ans+ar2
  }
  return(ans/L)
}

dip_chr<-function(x) {ans=x[,seq(1,dim(x)[2],2),]+x[,seq(2,dim(x)[2],2),];dim(ans)=c(dim(x)[1],dim(x)[2]/2,dim(x)[3]);ans}
dip_chr_ind<-function(x,ind) {hap<-c(ind*2-1,ind*2);ans=x[,hap[1],]+x[,hap[2],];dim(ans)=c(dim(x)[1],1,dim(x)[3]);ans}

dip_fr2_chr_ind<-function(x,y,ch,ind) 
{
  return(cor(c(dip_chr_ind(x[[ch]],ind)), c(dip_chr_ind(y[[ch]],ind)))^2)
}
dip_fr2_chr<-function(x,y,ch) 
{
  ans=0
  for (ind in 1:NUMI)
    ans=ans+dip_fr2_chr_ind(x,y,ch,ind)
  return(ans/NUMI)
}
dip_fr2<-function(x,y) 
{
  return(cor(c(unlist(lapply(x,dip_chr))),c(unlist(lapply(y,dip_chr))))^2)
}

hap_fr2_chr_k<-function(x,y,ch,k) 
{
  return(cor(c(x[[ch]][,k,]), c(y[[ch]][,k,]))^2)
}
hap_fr2_chr<-function(x,y,ch) 
{
  ans=0
  for (k in 1:NUMA)
    ans=ans+hap_fr2_chr_k(x,y,ch,k)
  return(ans/NUMA)
}
hap_fr2<-function(x,y) 
{
  return(cor(c(unlist(x)),c(unlist(y)))^2)
}

