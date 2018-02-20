if (!exists("colvec")) colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")
happlot<-function(ch,k,x,probs,ylab,mlab=paste("Haplotype", k),xlab=paste("Position on Chromosome", chrnos[ch]),cexa=1) { # probs is L*K*length(x) in dimension
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,1)
  plot(x=range(x),y=c(0,1),xlim=xlim,ylim=ylim,t='n',axes=F,ylab=ylab,main=mlab,xlab=xlab)
  upper=lower=rep(0,G)
  for (i in 1:L) 
  {
    upper=lower+probs[i,]
    polygon(x=x,y=c(lower,rev(upper)),col=colvec[i],xlim=xlim,ylim=ylim,border=NA)
    lower=upper
  }
  axis(2)
}
dipplot<-function(ch,ind,x,probs,ylab,mlab=paste("Individual",xlab=paste("Position on Chromosome", chrnos[ch]), ind),cexa=1) { # probs is L*K*length(x) in dimension
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,2)
  plot(x=range(x),y=c(0,2),xlim=xlim,ylim=ylim,t='n',axes=F,ylab=ylab,main=mlab,xlab=xlab)
  upper=lower=rep(0,G)
  for (i in 1:L) 
  {
    upper=lower+probs[i,]
    polygon(x=x,y=c(lower,rev(upper)),col=colvec[i],xlim=xlim,ylim=ylim,border=NA)
    lower=upper
  }
  axis(2)
}

happlot_Mu<-function(ch,k,x,probs,ylab,mlab=paste("Haplotype", k),xlab=paste("Position on Chromosome", chrnos[ch]),t.Mu,pow=4,cexa=1) { # probs is L*K*length(x) in dimension
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,1)
  plot(x=range(x),y=c(0,1),xlim=xlim,ylim=ylim,t='n',axes=F,ylab=ylab,main=mlab,xlab=xlab)
  alpha.Mu<-t.Mu^pow
  alpha.Mu<-t(t(alpha.Mu)/apply(alpha.Mu,2,max))
  upper=lower=rep(0,G)
  for (i in 1:L) 
    for (ll in 1:kLL) 
    {
      tmprgb=rgb(t(col2rgb(colvec[i])/255),alpha=alpha.Mu[ll,i]) 
      upper=lower+probs[[k]][,(i-1)*kLL+ll]
      polygon(x=x,y=c(lower,rev(upper)),col=tmprgb,border=NA)
      lower=upper
    }
  axis(2)
}

dipplot_Mu<-function(ch,ind,x,probs,ylab,mlab=paste("Individual", ind),xlab=paste("Position on Chromosome", chrnos[ch]),t.Mu,pow=4,cexa=1) { # probs is L*K*length(x) in dimension
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  hap<-c(ind*2-1,ind*2)
  G=length(x)
  glocs=rep(x,100)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,2)
  plot(x=range(x),y=c(0,2),xlim=xlim,ylim=ylim,t='n',axes=F,ylab=ylab,main=mlab,xlab=xlab)
  alpha.Mu<-t.Mu^pow
  alpha.Mu<-t(t(alpha.Mu)/apply(alpha.Mu,2,max))
  upper=lower=rep(0,G)
  for (i in 1:L) 
    for (ll in 1:kLL) 
    {
      tmprgb=rgb(t(col2rgb(colvec[i])/255),alpha=alpha.Mu[ll,i]) 
      upper=lower+probs[[hap[1]]][,(i-1)*kLL+ll]+probs[[hap[2]]][,(i-1)*kLL+ll]
      polygon(x=x,y=c(lower,rev(upper)),col=tmprgb,border=NA)
      lower=upper
    }
  axis(2)
}

report_donors<-function(MODE) 
{
  w<-locator(n=2) # find 4 corners of region of interest; first is top left, second is bottom right
  # add lines to the plot to highlight region
  #lines(c(w$x[1],w$x[2]), c(w$y[2],w$y[2]), col=6, lwd=4) # bottom line
  #lines(c(w$x[1],w$x[2]), c(w$y[1],w$y[1]), col=6, lwd=4) # top line
  #lines(c(w$x[1],w$x[1]), c(w$y[1],w$y[2]), col=6, lwd=4) # left line
  #lines(c(w$x[2],w$x[2]), c(w$y[1],w$y[2]), col=6, lwd=4) # right line
  x<-which.min((g.loc[[ch]]*1e-6-w$x[1])^2):which.min((g.loc[[ch]]*1e-6-w$x[2])^2)
  # y is more complicated. End points not enough as vector will change across x
  if (MODE=="HAP") 
    probs<-gfbs[[ch]][[k]][x,]
  if (MODE=="DIP") 
    probs<-gfbs[[ch]][[hap[1]]][x,]+gfbs[[ch]][[hap[2]]][x,]
  cprobs<-apply(probs,1,cumsum)
  plower<-pupper<-ylower<-yupper<-NULL # these will hold the indices of the donors that fall inside the highlighted region
  for (g in 1:length(x))
  {
    ylower<-c(ylower,which.min((cprobs[,g]-w$y[2])^2))
    plower<-c(plower,cprobs[ylower[g],g])
    yupper<-c(yupper,which.min((cprobs[,g]-w$y[1])^2))
    pupper<-c(pupper,cprobs[yupper[g],g])
  }
  # more accurate boxing
  lines(g.loc[[ch]][x]*1e-6, plower, col=7, lwd=4) # bottom line
  lines(g.loc[[ch]][x]*1e-6, pupper, col=7, lwd=4) # top line
  lines(c(g.loc[[ch]][x[1]]*1e-6,g.loc[[ch]][x[1]]*1e-6), c(plower[1],pupper[1]), col=7, lwd=4) # left line
  g=length(x)
  lines(c(g.loc[[ch]][x[g]]*1e-6,g.loc[[ch]][x[g]]*1e-6), c(plower[g],pupper[g]), col=7, lwd=4) # right line
  ind.indices<-rep(seq(1:kLL),L) 
  f<-function(g) # donors and probabilities of copying the donors 
  {
    tmp<-c(ylower[g]-1):yupper[g]
    panels<-rownames(Mu)[ind.indices[tmp[-1]]]
    donprobs<-diff(cprobs[tmp,g])
    ord<-sort(donprobs,dec=T,index=T)$ix # reorder from high to low contributors
    ord<-ord[donprobs[ord]>0] # and drop the non-contributors
    panels<-panels[ord]
    donprobs<-donprobs[ord]
    cbind(panels=panels,probs=donprobs) 
  }
  donors<-tapply(1:length(x), 1:length(x), f)
  list(w=w,x=x,ylower=ylower,yupper=yupper,donors=donors)
}
all_donors<-function(t.gfbs,t.localanc) 
{
  donors=list()
  for (ch in 1:nchrno)
  {
    donors[[ch]]=array(0,c(length(g.loc[[ch]]),kLL,L)) # contains proportion copied from each panel across the gridpoints
    for (a in 1:L)
    {
      for (ll in 1:kLL)
	for (i in 1:NUMA)
	  donors[[ch]][,ll,a]=donors[[ch]][,ll,a]+t.gfbs[[ch]][[i]][,(a-1)*kLL+ll]
      donors[[ch]][,,a]=donors[[ch]][,,a]/colSums(t.localanc[[ch]][a,,]) # divide by prob of that ancestry 
    }
    donors[[ch]]=donors[[ch]]/NUMA
  }
  return(donors)
}
# fit=apply(donors[[ch]][8898:8909,,],c(2,3),mean);rownames(fit)=rownames(Mu); a=2;head(sort(fit[,a],decreasing=T),2)
plot_panel_dist=function(donors,ch,a)
{
  m=colMeans(donors[[ch]][,,a]) # note difference from Chromosomal mean rather than genomewide mean
  d2=sqrt(colMeans((t(donors[[ch]][,,a])-m)^2))
  plot(g.loc[[ch]], d2, t='l', ylab="abs change in group copying", xlab="position", main=paste("chromosome", chrnos[ch]))
}
plot_Mu<-function(t.Mu=Mu, t.alpha=alpha, t.NL=NL, MODE="scaled", showgradient=F, beside=F, ord=T, pow=1, cexa=1, shiftl=cexa, shiftt=cexa, cutoff=0,tol=1e-6)
{
  t.kLL=nrow(t.Mu)
  t.alpha=Reduce("+",t.alpha)/length(t.alpha)
  if (MODE=="copy")
    t.Mu=t.Mu # plot copying matrix as is i.e. P(panel | anc)
  if (MODE=="scaled" | MODE=="jointscaled")
  {
    t.Mu=t.Mu/t.NL[1:t.kLL] # copying matrix scaled by number in panel i.e. P(hap in panel | anc)=P(hap|panel)P(panel|anc)
    t.Mu=t(t(t.Mu)/colSums(t.Mu)) # resacale so columns sum to one
  }
  # jointscaled is copying matrix times alpha / number in panel i.e. P(hap in panel,anc)=P(hap in panel|anc)P(anc) 
  if (MODE=="joint" | MODE=="jointscaled") 
    t.Mu=t(t(t.Mu)*t.alpha) # plot copying matrix times alpha i.e. P(panel,anc)=P(panel|anc)P(anc) 
  L=ncol(t.Mu)
  par(mar=c(3,4+shiftl,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
  t.Mu[t.Mu<tol]<-tol
  if (ord)
  {
    vord<-NULL
    for (l in 1:L)
    {
      maxforthis<-which(apply(t.Mu,1,which.max)==l)
      vord<-c(vord,maxforthis[sort(t.Mu[maxforthis,l],index=T)$ix])
    }
    t.Mu<-as.matrix(t.Mu[vord,])
  }
  if (!beside)
  {
    #t.Mu=t.Mu[rowMeans(t.Mu)>=quantile(rowMeans(t.Mu),cutoff),] # only those above cutoff on the quantiles on average
    tmp=rep(F,nrow(t.Mu));for (l in 1:L) tmp=(tmp | (t.Mu[,l]>=quantile(t.Mu[,l],cutoff)));t.Mu=t.Mu[tmp,] # only those above cutoff on the quantiles in some anc
  }
  if (!showgradient & !beside)
  {
    barplot(t(t.Mu),space=F,col=colvec,horiz=T,las=T,cex.axis=cexa,cex.names=cexa)
    for (i in 1:L)
      mtext(side=3, at=(i)*0.2*max(rowSums(t.Mu)), round(t.alpha[i],3), cex=cexa, col=colvec[i])
  }
  if (!showgradient & beside)
  {
    if (ord) par(mfrow=c(1,L))
    if (!ord) 
    {
      par(mar=c(3,2,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
      nf <- layout(matrix(c(1:(L+1)),ncol=L+1), widths=c(1,rep(2,L)), TRUE);layout.show(nf)
      plot(c(-cexa,ncol(t.Mu)),c(0.5,nrow(t.Mu)+0.5),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
      text(x=2,pos=2,y=(1:nrow(t.Mu)),rownames(t.Mu),cex=1.75*cexa)
    }
    all.alpha=Reduce("+",alpha);all.alpha=all.alpha/sum(all.alpha);all.alpha=round(all.alpha,3)
    if (ord)
    {
      allord=apply(t.Mu,2,order)
      ordMu=NULL;for (a in 1:L) ordMu[[a]]=t.Mu[allord[,a],a] 
    } else 
    {
      ordMu=NULL;for (a in 1:L) ordMu[[a]]=t.Mu[,a] 
    }
    if (cutoff > 0)
    {
      for (a in 1:L)
	ordMu[[a]]=ordMu[[a]][ordMu[[a]]>=quantile(ordMu[[a]],cutoff)] # only show those that are above cutoff on the quantiles
    }
    xmax=max(unlist(ordMu))
    for (a in 1:L) 
    {
      if (ord) y=barplot(ordMu[[a]],horiz=T,las=1,col=colvec[a],xlim=c(0,xmax),cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa*2)
      if (!ord)
      {
	y=barplot(ordMu[[a]],horiz=T,las=1,col=colvec[a],xlim=c(0,xmax),cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa*2,names.arg=rep("",length(ordMu[[a]])))
      }
      text(xmax/4,max(y)+shiftt,all.alpha[a],col=colvec[a],cex=2*cexa)
    }
  }
  if (showgradient) # overrides beside
  {
    alpha.Mu<-t.Mu^pow
    alpha.Mu<-t(t(alpha.Mu)/apply(alpha.Mu,2,max))
    plot(c(-cexa,ncol(t.Mu)),c(0.5,nrow(t.Mu)+0.5),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
    mtext(side=3, at=1:L-0.5, round(t.alpha,3), cex=cexa)
    for (l in 1:ncol(t.Mu))
      for (ll in 1:nrow(t.Mu))
      {
	tmprgb=rgb(t(col2rgb(colvec[l])/255),alpha=alpha.Mu[ll,l])
	polygon(x=c(l,l+1,l+1,l)-1,y=c(ll,ll,ll+1,ll+1)-0.5,col=tmprgb,border=NA)
      }
    text(x=0,pos=2,y=(1:nrow(t.Mu)),rownames(t.Mu),cex=0.75*cexa)
  }
  if (ord) 
    if (!beside)
      return(t.Mu) # return as re-ordered version is useful for dipplot_Mu and happlot_Mu
    if (beside)
      return(ordMu) # return as re-ordered version is useful for dipplot_Mu and happlot_Mu
}
plot_Fst<-function(t.Fst, ord=T, cexa=1, shiftl=cexa, shiftt=cexa, cutoff=nrow(t.Fst))
{
  t.kLL=nrow(t.Fst)
  if ((!ord) & (cutoff!=t.kLL))
  {
    warning("showing all as re-ordering not allowed")
    cutoff=t.kLL
  }
  L=ncol(t.Fst)
  par(mar=c(3,4+shiftl,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
  if (ord) par(mfrow=c(1,L))
  if (!ord) 
  {
    par(mar=c(3,2,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
    nf <- layout(matrix(c(1:(L+1)),ncol=L+1), widths=c(1,rep(2,L)), TRUE);layout.show(nf)
    plot(c(-cexa,ncol(t.Fst)),c(0.5,nrow(t.Fst)+0.5),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
    text(x=2,pos=2,y=(1:nrow(t.Fst)),rownames(t.Fst),cex=1.75*cexa)
  }
  if (ord)
  {
    allord=apply(t.Fst,2,order,decreasing=T)
    ordFst=NULL;for (a in 1:L) ordFst[[a]]=t.Fst[allord[,a],a] 
  } else 
  {
    ordFst=NULL;for (a in 1:L) ordFst[[a]]=t.Fst[,a] 
  }
  for (a in 1:L)
    ordFst[[a]]=tail(ordFst[[a]],cutoff) # only show those that are below cutoff place
  for (a in 1:L) 
  {
    if (ord) y=barplot(ordFst[[a]],horiz=T,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa*2)
    if (!ord)
      y=barplot(ordFst[[a]],horiz=T,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa*2,names.arg=rep("",length(ordFst[[a]])))
  }
  if (ord) 
    return(ordFst)
}

plot_r2=function(r2file)
{
  tmp=read.table(r2file)
  r2=tmp[,2]
  names(r2)=paste0(tapply(1:length(r2),1:length(r2),function(i) strsplit(as.character(tmp[i,1]),"way")[[1]][1]),"way")
  par(mar=c(4,10,0,2))
  barplot(sort(r2),horiz=T,las=2,xlim=c(0,1))
  r2
}


# localanc is w.r.t. final phasing. The below undoes this to compare with input phasing
phase_localanc=function(t.localanc,t.flips) 
{
  for (ch in 1:nchrno)
    for (ind in 1:(NUMA/2))
    {
      haps=c(ind*2-1,ind*2)
      ans=t.localanc[[ch]][,haps,]
      for (k in haps)
      {
	otherhap=ifelse((k%%2)==0,1,2) # if even look at previous; if odd look at next
	kf=t.flips[[ind]][[ch]] # F and T values
	for (l in 1:L)
	  t.localanc[[ch]][l,k,kf]=ans[l,otherhap,kf]
      }
    }
  return(t.localanc)
}

# gfbs is w.r.t. final phasing. The below undoes this to compare with input phasing
phase_gfbs=function(t.gfbs,t.flips)
{
  for (ch in 1:nchrno)
    for (ind in 1:(NUMA/2))
    {
      haps=c(ind*2-1,ind*2)
      ans=list(t.gfbs[[ch]][[haps[1]]],t.gfbs[[ch]][[haps[2]]])
      for (k in haps)
      {
	otherhap=ifelse((k%%2)==0,1,2) # if even look at previous; if odd look at next
	kf=t.flips[[ind]][[ch]] # F and T values
	for (kl in 1:(L*kLL))
	  t.gfbs[[ch]][[k]][kf,kl]=ans[[otherhap]][kf,kl]
      }
    }
  return(t.gfbs)
}

plot_Mu_invFst=function(fst_panels,t.L=L,t.Mu=Mu,t.alpha=alpha,t.NL=NL,fst_tol=2e-3) 
{
  for (l in 1:t.L)
    if (min(fst_panels[,l])<0) # rescale this side's Fst values s.t. min is bigger than zero
    {
      fst_panels[,l]=fst_panels[,l]-min(fst_panels[,l])+fst_tol
    }
  t.kLL=nrow(t.Mu)
  t.alpha=Reduce("+",t.alpha)/length(t.alpha)
  t.Mu=t.Mu/t.NL[1:t.kLL] # copying matrix scaled by number in panel i.e. P(hap in panel | anc)=P(hap|panel)P(panel|anc)
  t.Mu=t(t(t.Mu)/colSums(t.Mu))
  #t.Mu=t(t(t.Mu)*t.alpha) # copying matrix times alpha i.e. P(panel,anc)=P(panel|anc)P(anc) 
  #t.Mu=t.Mu/sum(t.Mu) # normalise
  invFst=1/fst_panels;rownames(invFst)=rownames(t.Mu)
  #plot(c(t.Mu),c(invFst),pch=20,xlab="Mu",ylab=bquote(paste(F^{-1},""[.("st")])))  
  plot(c(t.Mu),c(invFst),pch=20,xlab="Mu",ylab=bquote(paste("(",F[.("st")],")",""^{-1})),col=c(rep(1,t.kLL),rep(2,t.kLL)))  
  return(list("tMu"=t.Mu, "invFst"=invFst))
}
#source("fst.R")
#fst_panels=matrix(NaN,t.kLL,t.L)
#for (l1 in 1:t.kLL)
#{
#  load(paste0("FREQS/", rownames(t.Mu)[l1],"_freqs.rdata"))
#  for (l2 in 1:t.L) 
#    fst_panels[l1,l2]=wc_fst(ancestral_freqs$freqs[[l2]],ancestral_freqs$counts[[l2]],pdata$freqs, pdata$counts)
#}

#load("all_Fst_2.rdata");i=49;filename=dir(pattern=glob2rx(paste0(paste0(strsplit(names(all_Fst)[i],"way")[[1]],collapse="*"),"*.RData")),path="RESULTS/");load(paste0("RESULTS/",filename))
#tmp=plot_Mu_invFst(all_Fst[[i]]$panels)
#cor.test(tmp$tMu,tmp$invFst)
