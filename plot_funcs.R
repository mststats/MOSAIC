# various functions that create useful plots of MOSAIC model fit
# such as local ancestry plots, plots of the inferred copying matrix, tables of top donors in selected regions, etc

happlot<-function(ch,k,x,L,probs,ylab,mlab=paste("Haplotype", k),xlab="",cexa=1,
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  # probs is L*K*length(x) in dimension
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,1)
  plot(x=range(x),y=c(0,1),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
  upper=lower=rep(0,G)
  for (i in 1:L) 
  {
    upper=lower+probs[i,]
    polygon(x=x,y=c(lower,rev(upper)),col=colvec[i],xlim=xlim,ylim=ylim,border=NA)
    lower=upper
  }
  axis(2)
}
dipplot<-function(ch,ind,x,L,probs,ylab,mlab=paste("Individual",ind),xlab="",cexa=1,
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,2)
  plot(x=range(x),y=c(0,2),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
  upper=lower=rep(0,G)
  for (i in 1:L) 
  {
    upper=lower+probs[i,]
    polygon(x=x,y=c(lower,rev(upper)),col=colvec[i],xlim=xlim,ylim=ylim,border=NA)
    lower=upper
  }
  axis(2)
}

happlot_Mu<-function(ch,k,x,L,probs,ylab,mlab=paste("Haplotype", k),xlab="",t.Mu,pow=4,cexa=1,
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  kLL=nrow(Mu)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,1)
  plot(x=range(x),y=c(0,1),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
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

dipplot_Mu<-function(ch,ind,x,L,probs,ylab,mlab=paste("Individual", ind),xlab="",t.Mu,pow=4,cexa=1, 
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  hap<-c(ind*2-1,ind*2)
  G=length(x)
  kLL=nrow(Mu)
  glocs=rep(x,100)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,2)
  plot(x=range(x),y=c(0,2),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
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
    ord<-sort(donprobs,dec=TRUE,index=TRUE)$ix # reorder from high to low contributors
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
# fit=apply(donors[[ch]][8898:8909,,],c(2,3),mean);rownames(fit)=rownames(Mu); a=2;head(sort(fit[,a],decreasing=TRUE),2)
plot_panel_dist=function(donors,ch,a)
{
  m=colMeans(donors[[ch]][,,a]) # note difference from Chromosomal mean rather than genomewide mean
  d2=sqrt(colMeans((t(donors[[ch]][,,a])-m)^2))
  plot(g.loc[[ch]], d2, t='l', ylab="abs change in group copying", xlab="position", main=paste("chromosome", t.chrnos[ch]))
}
plot_Mu<-function(t.Mu, t.alpha, t.NL, MODE="scaled", showgradient=FALSE, beside=TRUE, ord=TRUE, pow=1, cexa=1.5, shiftl=cexa*2, shiftt=cexa*2, 
		  cutoff=0,tol=1e-6, colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
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
  par(mar=c(3,4+shiftl,3+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
  t.Mu[t.Mu<tol]<-tol
  if (ord)
  {
    vord<-NULL
    for (l in 1:L)
    {
      maxforthis<-which(apply(t.Mu,1,which.max)==l)
      vord<-c(vord,maxforthis[sort(t.Mu[maxforthis,l],index=TRUE)$ix])
    }
    t.Mu<-as.matrix(t.Mu[vord,])
  }
  if (!beside)
  {
    #t.Mu=t.Mu[rowMeans(t.Mu)>=quantile(rowMeans(t.Mu),cutoff),] # only those above cutoff on the quantiles on average
    tmp=rep(FALSE,nrow(t.Mu));for (l in 1:L) tmp=(tmp | (t.Mu[,l]>=quantile(t.Mu[,l],cutoff)));t.Mu=t.Mu[tmp,] # only those above cutoff on the quantiles in some anc
  }
  if (!showgradient & !beside)
  {
    barplot(t(t.Mu),space=FALSE,col=colvec,horiz=TRUE,las=TRUE,cex.axis=cexa,cex.names=cexa)
    for (i in 1:L)
      mtext(side=3, at=(i)*0.2*max(rowSums(t.Mu)), round(t.alpha[i],3), cex=cexa, col=colvec[i])
  }
  if (!showgradient & beside)
  {
    if (ord) par(mfrow=c(1,L))
    par(mar=c(3,1.5*shiftl,1,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
    if (!ord) 
    {
      nf <- layout(matrix(c(1:(L+1)),ncol=L+1), widths=c(1,rep(2,L)), TRUE);#layout.show(nf)
      plot(c(-cexa*2,ncol(t.Mu)),c(1,nrow(t.Mu)+shiftt),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
      text(x=cexa*2,pos=2,y=0.5+1.2*(1:nrow(t.Mu)),rownames(t.Mu),cex=1.75*cexa)
    }
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
      if (ord) y=barplot(ordMu[[a]],horiz=TRUE,las=1,col=colvec[a],xlim=c(0,xmax),ylim=c(0,length(ordMu[[a]])+shiftt/2),cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa*2)
      if (!ord)
	y=barplot(ordMu[[a]],horiz=TRUE,las=1,col=colvec[a],xlim=c(0,xmax),cex.names=cexa,cex.axis=cexa,main="",ylim=c(0,length(ordMu[[a]])+shiftt/2),
		  cex.main=cexa*2,names.arg=rep("",length(ordMu[[a]])))
      text(xmax/4,(length(ordMu[[a]])+shiftt/2),round(t.alpha[a],3),col=colvec[a],cex=cexa)
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
plot_Fst<-function(t.Fst, ord=TRUE, cexa=1, shiftl=cexa, shiftt=cexa, cutoff=nrow(t.Fst),
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  t.kLL=nrow(t.Fst)
  if ((!ord) & (cutoff!=t.kLL))
  {
    warning("showing all as re-ordering not allowed")
    cutoff=t.kLL
  }
  if (cutoff>t.kLL) cutoff=t.kLL
  L=ncol(t.Fst)
  par(mar=c(6,4+shiftl,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
  if (ord) par(mfrow=c(1,L))
  if (!ord) 
  {
    par(mar=c(6,2,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
    nf <- layout(matrix(c(1:(L+1)),ncol=L+1), widths=c(1,rep(2,L)), TRUE);layout.show(nf)
    plot(c(-cexa,ncol(t.Fst)),c(0.5,nrow(t.Fst)+0.5),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
    text(x=2,pos=2,y=(1:nrow(t.Fst)),rownames(t.Fst),cex=cexa)
  }
  if (ord)
  {
    allord=apply(t.Fst,2,order,decreasing=TRUE)
    ordFst=NULL;for (a in 1:L) ordFst[[a]]=t.Fst[allord[,a],a] 
  } else 
  {
    ordFst=NULL;for (a in 1:L) ordFst[[a]]=t.Fst[,a] 
  }
  rangeFst=NULL
  for (a in 1:L) {
    ordFst[[a]]=tail(ordFst[[a]],cutoff) # only show those that are below cutoff place
    rangeFst[[a]]=range(1-ordFst[[a]])
  }
  for (a in 1:L) 
  {
    if (ord) {
      plot(c(0,1),c(0,cutoff+1),t='n',yaxt='n',ylab="",xlab="1-Fst",cex.main=cexa,main="",cex.axis=cexa,xlim=c(rangeFst[[a]][1],rangeFst[[a]][2]))
      tmp=1-ordFst[[a]]-rangeFst[[a]][1]
      y=barplot(tmp,horiz=TRUE,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa,add=T,offset=rangeFst[[a]][1])
    }
    if (!ord)
      y=barplot(1-ordFst[[a]]-rangeFst[[a]][1],horiz=TRUE,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",xlab="1-Fst",
		cex.main=cexa,names.arg=rep("",length(ordFst[[a]])),offset=rangeFst[[a]][1])
  }
  if (ord) 
    return(ordFst)
}

# localanc is w.r.t. final phasing. The below undoes this to compare with input phasing (necessary for Fst calculations)
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
	kf=t.flips[[ind]][[ch]] # FALSE and TRUE values
	for (l in 1:L)
	  t.localanc[[ch]][l,k,kf]=ans[l,otherhap,kf]
      }
    }
  return(t.localanc)
}

# function to plot the local ancestry of each target admixed genome along each chromosome
plot_localanc=function(t.chrnos, t.g.loc, t.localanc, t.g.true_anc=NULL,cexa=2,pow=1,y.lab="expected",MODE="BAR",NCHR=2,PAUSE=T,t.Mu=NULL,t.gfbs=NULL,
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  if (is.null(t.Mu) & MODE=="GRAD") stop("Please supply a copying matrix Mu for use with this plot")
  if (is.null(t.gfbs) & MODE=="GRAD") stop("Please supply a full array of posterior probabilities for use with this plot")
  G=sapply(t.localanc, function(x) dim(x)[3])
  NUMA=dim(t.localanc[[1]])[2]
  L=dim(t.localanc[[1]])[1]
  nchrno=length(t.chrnos)
  NUMI=NUMA/2
  if (NCHR==2 & NUMA==1)
  {
    warning("changing NCHR to HAP as only on hap in target",immediate.=T)
    NCHR=1
  }
  if (NCHR==1)
  {
    expected_r2<-matrix(NaN,nchrno,NUMA)
    if (!is.null(t.g.true_anc))
    {
      perc<-matrix(NaN,nchrno,NUMA)
      mse<-matrix(NaN,nchrno,NUMA)
      r2<-matrix(NaN,nchrno,NUMA)
    }
    for (ch in 1:nchrno)
      for (k in 1:NUMA)
      {
	ind=as.integer((k+1)/2)
	expected_r2[ch,k]=hap_expected_fr2_chr_k(t.localanc, ch, k)
	cat("E[r2]: ", expected_r2[ch,k], "  ", sep="")
	if (is.null(t.g.true_anc)) par(mfrow=c(1,1),mar=c(4, 5.2, 1, 0)+0.1)
	if (!is.null(t.g.true_anc)) par(mfrow=c(2,1),mar=c(4, 5.2, 1, 0)+0.1)
	if (!is.null(t.g.true_anc))
	{
	  # how many maximal ancestry predictions are the same as the truth
	  perc[ch,k]<-round(100*mean(apply(matrix(t.localanc[[ch]][,k,],L),2,which.max)==
				     apply(matrix(t.g.true_anc[[ch]][,k,],L),2,which.max)),2)
	  mse[ch,k]<-mean((t.localanc[[ch]][,k,]-t.g.true_anc[[ch]][,k,])^2)
	  r2[ch,k]=hap_fr2_chr_k(t.localanc, t.g.true_anc, ch, k)
	  cat("r2: ",r2[ch,k],"  ", sep="")
	  cat("correct: ", perc[ch,k], "%   ", sep="")
	  cat("MSE:", mse[ch,k], "  ", sep="")
	  cat("\n")
	  par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	  if (MODE=="LINE")
	  {
	    plot(range(t.g.loc[[ch]])*1e-6,c(0,1),axes=F,t='n',ylab="truth",main=paste("Haplotype", k), xlab=paste("Position on Chromosome",t.chrnos[ch]))
	    for (i in 1:L) lines(t.g.loc[[ch]]*1e-6, t.g.true_anc[[ch]][i,k,], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa);	
	    axis(2)
	  }
	  if (MODE=="BAR" | MODE=="GRAD")
	    happlot(ch,k,t.g.loc[[ch]]*1e-6,L,t.g.true_anc[[ch]][,k,],ylab="truth",cexa=cexa)
	}
	par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	if (MODE=="LINE")
	{
	  plot(range(t.g.loc[[ch]])*1e-6,c(0,1),axes=F,t='n',ylab=y.lab,main=paste("Haplotype", k), xlab=paste("Position on Chromosome",t.chrnos[ch]))
	  for (i in 1:L) lines(t.g.loc[[ch]]*1e-6, t.localanc[[ch]][i,k,], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa);	
	  axis(2)
	}
	if (MODE=="BAR")
	  happlot(ch,k,t.g.loc[[ch]]*1e-6,L,t.localanc[[ch]][,k,],xlab=paste("Position on Chromosome",t.chrnos[ch]),ylab=y.lab,cexa=cexa)
	if (MODE=="GRAD")
	  happlot_Mu(ch,k,t.g.loc[[ch]]*1e-6,L,t.gfbs[[ch]],xlab=paste("Position on Chromosome",t.chrnos[ch]),ylab=y.lab,t.Mu=t.Mu,pow=pow,cexa=cexa)
	mp<-axTicks(1,round(axp=c(min(t.g.loc[[ch]])*1e-6,max(t.g.loc[[ch]])*1e-6,5)))
	axis(1,at=mp,labels=signif(mp,3))
	if (PAUSE) readline()
      }
  }
  if (NCHR==2)
  {
    expected_r2<-matrix(NaN,nchrno,NUMI)
    if (!is.null(t.g.true_anc))
    {
      perc<-matrix(NaN,nchrno,NUMI)
      mse<-matrix(NaN,nchrno,NUMI)
      r2<-matrix(NaN,nchrno,NUMI)
    }
    for (ch in 1:nchrno)
      for (ind in 1:(NUMI))
      {
	hap<-c(ind*2-1,ind*2)
	expected_r2[ch,ind]=dip_expected_fr2_chr_ind(t.localanc, ch, ind)
	cat("E[r2]: ", expected_r2[ch,ind], "  ", sep="")
	if (is.null(t.g.true_anc)) par(mfrow=c(1,1),mar=c(4, 5.2, 1, 0)+0.1)
	if (!is.null(t.g.true_anc)) par(mfrow=c(2,1),mar=c(4, 5.2, 1, 0)+0.1)
	if (!is.null(t.g.true_anc))
	{
	  perc[ch,ind]<-0 # how many maximal ancestry counts 0,1,2 are the same as the truth
	  for (i in 1:L)
	    perc[ch,ind]<-perc[ch,ind]+mean(as.integer(t.localanc[[ch]][i,hap[1],]+t.localanc[[ch]][i,hap[2],]+0.5)==
					    as.integer(t.g.true_anc[[ch]][i,hap[1],]+t.g.true_anc[[ch]][i,hap[2],]+0.5)) # push to int required here due to shift to the grid
	  perc[ch,ind]<-round(100*perc[ch,ind]/L,2)
	  mse[ch,ind]<-mean((t.localanc[[ch]][,hap[1],]+t.localanc[[ch]][,hap[2],]-t.g.true_anc[[ch]][,hap[1],]-t.g.true_anc[[ch]][,hap[2],])^2)
	  r2[ch,ind]=dip_fr2_chr_ind(t.localanc, t.g.true_anc, ch, ind)
	  cat("r2: ",r2[ch,ind],"  ", sep="")
	  cat("% correct: ", perc[ch,ind], "  ", sep="")
	  cat("MSE: ", mse[ch,ind], "  ", sep="")
	  cat("\n")
	  par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	  if (MODE=="LINE")
	  {
	    plot(range(t.g.loc[[ch]])*1e-6,c(0,2),axes=F,t='n',ylab="truth",main=paste("Individual", ind),xlab=paste("Position on Chromosome",t.chrnos[ch]))
	    for (i in 1:L) lines(t.g.loc[[ch]]*1e-6, t.g.true_anc[[ch]][i,hap[1],]+t.g.true_anc[[ch]][i,hap[2],], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa)
	    axis(2)
	  }
	  if (MODE=="BAR" | MODE=="GRAD")
	    dipplot(ch,ind,t.g.loc[[ch]]*1e-6,L,t.g.true_anc[[ch]][,hap[1],]+t.g.true_anc[[ch]][,hap[2],],xlab=paste("Position on Chromosome",t.chrnos[ch]),
		    ylab="truth",cexa=cexa)
	}
	if (MODE=="LINE")
	{
	  plot(range(t.g.loc[[ch]])*1e-6,c(0,2),axes=F,t='n',ylab=y.lab,main=paste("Individual", ind),xlab=paste("Position on Chromosome",t.chrnos[ch]))
	  for (i in 1:L) lines(t.g.loc[[ch]]*1e-6, t.localanc[[ch]][i,hap[1],]+t.localanc[[ch]][i,hap[2],], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa)
	  axis(2)
	}
	if (MODE=="BAR")
	  dipplot(ch,ind,t.g.loc[[ch]]*1e-6,L,t.localanc[[ch]][,hap[1],]+t.localanc[[ch]][,hap[2],],xlab=paste("Position on Chromosome",t.chrnos[ch]),ylab=y.lab,cexa=cexa)
	if (MODE=="GRAD")
	  xy<-dipplot_Mu(ch,ind,t.g.loc[[ch]]*1e-6,L,t.gfbs[[ch]],xlab=paste("Position on Chromosome",t.chrnos[ch]),ylab=y.lab,t.Mu=t.Mu,pow=pow,cexa=cexa)
	mp<-axTicks(1,round(axp=c(min(t.g.loc[[ch]])*1e-6,max(t.g.loc[[ch]])*1e-6,5)))
	axis(1,at=mp,labels=signif(mp,3))
	if (PAUSE) readline()
      }
  }
}

plot_mean_localanc=function(ch, chrnos, g.loc, localanc, whichhaps=1:dim(localanc[[ch]])[2], whicha=1:dim(localanc[[ch]])[1], cexa=2, lega=FALSE, 
			    ret=FALSE, colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  m=list()
  A=dim(localanc[[ch]])[1]
  for (a in 1:A)
    m[[a]]=colMeans(localanc[[ch]][a,whichhaps,]) # mean over target haplotypes
  plot(c(g.loc[[ch]][1],g.loc[[ch]][length(g.loc[[ch]])]),c(0,1),t='n',
       ylab="mean ancestry",xlab=paste("Position on Chromosome",chrnos[ch]),main="")  
  for (a in whicha)
  {
    lines(g.loc[[ch]],m[[a]],lwd=cexa,col=colvec[a])
    abline(h=c(mean(m[[a]])-1.96*sd(m[[a]]),mean(m[[a]])+1.96*sd(m[[a]])),lwd=cexa,lty=2,col=colvec[a])
  }
  if (lega)
    legend("topright", legend=whicha, col=colvec[whicha], lwd=cexa)
  if (ret)
    return(m)
}

# function to plot most useful figures
plot_all_mosaic=function(pathout="MOSAIC_PLOTS/",pathin="HGDP/",GpcM=60,doFst=TRUE) {
  targetdetails=paste0(target, "_", L, "way_", NUMA, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
		       "_",sum(NL),"_",GpcM)
  pdf(file=paste0(pathout,targetdetails,"_Mu.pdf"), width=12, height=7)
  ord.Mu=plot_Mu(Mu,alpha,NL,cexa=1.5,beside=T,shiftl=5,shiftt=2,cutoff=0,ord=T)
  dev.off()
  
  # note that it takes a while to calculate frequencies, etc
  if (doFst) {
    flocalanc=phase_localanc(localanc,final.flips) 
    if (target=="simulated")
      write_admixed_summary(targetdatasource="MOSAIC_RESULTS/",datasource=pathin,t.localanc=flocalanc,chrnos=chrnos)
    if (target!="simulated")
      write_admixed_summary(targetdatasource=pathin,datasource=pathin,t.localanc=flocalanc,chrnos=chrnos)
    write_panel_summaries(panels=rownames(Mu),chrnos=chrnos)
    all_Fst=Fst_combos(target, L, sum(NL), rownames(Mu))
    pdf(file=paste0(pathout,targetdetails,"_Fst.pdf"), width=21, height=28)
    ord.Fst=plot_Fst(all_Fst$panels,cexa=2,ord=T, shiftl=6, cutoff=10)
    dev.off()
  }
  
  # dimensions of plots
  d1=switch(L,NaN,1,2,2,3,3) 
  d2=switch(L,NaN,3,3,5,5,7) 
  pdf(file=paste0(pathout,targetdetails,"_acoanc.pdf"), width=5*d2,height=5*d1)
  this_acoplots=plot_coanccurves(acoancs,dr,lwd=4,cexa=2,verbose=F,axisall=F,samedates=F,asym=F,min.cM=0.5)
  dev.off()
  EMlog=extract_log(logfile)
  pdf(file=paste0(pathout,targetdetails,"_EMlog.pdf"), width=14,height=14)
  plot_loglike(EMlog)
  dev.off()
}

