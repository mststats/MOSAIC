#ch=1
if (!exists("minS")) minS=10
if (!exists("ch")) ch="all"
if (!exists("pops")) pops="all"
#a=2
cols=c(1,4)
HLA=c(28510120,33480577)
require(boot)
#pathin="HGDP/";source("NAlabels.R")
whichhaps=1:NUMA
if (all(pops!="all"))
  whichhaps=which(hap_NAlabels%in%which(NAfs%in%pops))
N=length(localanc[[1]][a,,1])
par(mfrow=c(2,1),cex.lab=2,mar=c(6,5,2,2))
G=sapply(g.loc,length)
add.loc=0;if (nchrno>1) for (tch in 2:nchrno) add.loc[tch]=add.loc[tch-1]+g.loc[[tch-1]][G[tch-1]]-g.loc[[tch]][1]
m=nlp=list()
for (tch in 1:nchrno) 
  m[[tch]]=colMeans(localanc[[tch]][a,whichhaps,]) # mean over inds
mm=mean(unlist(m));Nm=sum(G)
sm=sd(unlist(m))
for (tch in 1:nchrno) 
{
  #nlp[[tch]]=tapply(1:G[tch],1:G[tch],function(g) -log(1-pt(abs(m[[tch]][g]-mm)/(sd(localanc[[tch]][a,,g])/sqrt(N)),N-1),10))
  nlp[[tch]]=tapply(1:G[tch],1:G[tch],function(g) -log(1-pnorm(abs(m[[tch]][g]-mm)/sm),10))
  nlp[[tch]][is.infinite(nlp[[tch]])]=20 # if too far out
  #tmpm=rowMeans(localanc[[tch]][a,,])
  #nlp[[tch]]=tapply(1:G[tch],1:G[tch],function(g) {tmp=localanc[[tch]][a,,g]-tmpm;-log(1-pt(mean(tmp)/(sd(tmp)/sqrt(N)),N-1),10)})
}
blockout=function(t.locs,t.minS,lower,upper,t.add.loc,t.col="white")
{
  for (s in 1:ceiling(locs[length(locs)])) # 1Mb regions
  {
    tmp.loc=which((locs>(s-1)) & (locs<s))
    if (length(tmp.loc) < minS)
      polygon(x=1e6*c(s-1,s,s,s-1)+t.add.loc,y=c(lower,lower,upper,upper),col=t.col,border=t.col)
  }
}
if (ch!="all")
{
  plot(g.loc[[ch]],m[[ch]],t='l',ylim=c(min(m[[ch]]),max(m[[ch]])),ylab="mean African ancestry",xlab=paste("Position on Chromosome",chrnos[ch]),main="")
  # block out 1Mb regions with fewer than minS markers
  snps=t(matrix(scan(paste0("HGDP/","snpfile.",ch),what="character"),nrow=6))
  locs=as.integer(snps[,4])*1e-6
  blockout(locs,minsS,min(m[[ch]]),max(m[[ch]]),0)
  abline(h=mm)
  abline(h=mm-2*sm,lty=2)
  abline(h=mm+2*sm,lty=2)
  if (chrnos[ch]==6) abline(v=HLA,col=2,lwd=2)
  plot(g.loc[[ch]],nlp[[ch]],t='l',ylab=bquote(paste(-log[.(10)],"p")),xlab=paste("Position on Chromosome",chrnos[ch]),main="")
  blockout(locs,minsS,min(nlp[[ch]]),max(nlp[[ch]]),0)
  if (chrnos[ch]==6) abline(v=HLA,col=2,lwd=2)
}

if (ch=="all")
{
  l=min(unlist(m));u=max(unlist(m))
  XLIM=c(g.loc[[1]][1],add.loc[nchrno]+g.loc[[nchrno]][G[nchrno]])
  par(xaxs="i",xaxt="n")
  plot(XLIM,c(l,u),t='n',ylab="mean African ancestry",xlab="chromosome")
  for (tch in 1:nchrno)
  {
    lines(g.loc[[tch]]+add.loc[tch], m[[tch]],col=cols[tch%%2+1])
    snps=t(matrix(scan(paste0("HGDP/","snpfile.",tch),what="character"),nrow=6))
    locs=as.integer(snps[,4])*1e-6
    blockout(locs,minsS,l,u,add.loc[[tch]])
    mtext(chrnos[tch],side=1,at=mean(g.loc[[tch]])+add.loc[tch],col=cols[tch%%2+1],cex=1.2,line=1.2)
    #if (chrnos[tch]==6) abline(v=add.loc[tch]+HLA,col=2,lwd=2)
  }
  abline(h=mm)
  abline(h=mm-2*sm,lty=2)
  abline(h=mm+2*sm,lty=2)
  
  pl=min(unlist(nlp));pu=max(unlist(nlp))
  plot(XLIM,c(pl,pu),t='n',ylab=bquote(paste(-log[.(10)],"p")),xlab="chromosome")
  for (tch in 1:nchrno)
  {
    lines(g.loc[[tch]]+add.loc[tch], nlp[[tch]],col=cols[tch%%2+1])
    snps=t(matrix(scan(paste0("HGDP/","snpfile.",tch),what="character"),nrow=6))
    locs=as.integer(snps[,4])*1e-6
    blockout(locs,minsS,pl,pu,add.loc[[tch]])
    mtext(chrnos[tch],side=1,at=mean(g.loc[[tch]])+add.loc[tch],col=cols[tch%%2+1],cex=1.2,line=1.2)
    #if (chrnos[tch]==6) abline(v=add.loc[tch]+HLA,col=2,lwd=2)
  }
}
