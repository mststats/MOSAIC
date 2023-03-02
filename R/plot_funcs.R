# various functions that create useful plots of MOSAIC model fit
# such as local ancestry plots, plots of the inferred copying matrix, tables of top donors in selected regions, etc

happlot<-function(ch,k,x,A,probs,ylab,mlab=paste("Haplotype", k),xlab="",cexa=1,
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  # probs is A*K*length(x) in dimension
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,1)
  plot(x=range(x),y=c(0,1),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
  upper=lower=rep(0,G)
  for (i in 1:A) 
  {
    upper=lower+probs[[ch]][i,k,]
    polygon(x=x,y=c(lower,rev(upper)),col=colvec[i],xlim=xlim,ylim=ylim,border=NA)
    lower=upper
  }
  axis(2)
}
dipplot<-function(ch,ind,x,A,probs,ylab,mlab=paste("Individual",ind),xlab="",cexa=1,
		  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,2)
  hap=c(ind*2-1,ind*2)
  plot(x=range(x),y=c(0,2),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
  upper=lower=rep(0,G)
  for (i in 1:A) 
  {
    upper=lower+probs[[ch]][i,hap[1],]+probs[[ch]][i,hap[2],]
    polygon(x=x,y=c(lower,rev(upper)),col=colvec[i],xlim=xlim,ylim=ylim,border=NA)
    lower=upper
  }
  axis(2)
}

happlot_Mu<-function(ch,k,x,A,probs,ylab,mlab=paste("Haplotype", k),xlab="",t.Mu,pow=4,cexa=1,
		     colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  G=length(x)
  kLL=nrow(t.Mu)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,1)
  plot(x=range(x),y=c(0,1),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
  alpha.Mu<-t.Mu^pow
  alpha.Mu<-t(t(alpha.Mu)/apply(alpha.Mu,2,max))
  upper=lower=rep(0,G)
  for (i in 1:A) 
    for (ll in 1:kLL) 
    {
      tmprgb=rgb(t(col2rgb(colvec[i])/255),alpha=alpha.Mu[ll,i]) 
      upper=lower+probs[[ch]][[k]][,(i-1)*kLL+ll]
      polygon(x=x,y=c(lower,rev(upper)),col=tmprgb,border=NA)
      lower=upper
    }
  axis(2)
}

dipplot_Mu<-function(ch,ind,x,A,probs,ylab,mlab=paste("Individual", ind),xlab="",t.Mu,pow=4,cexa=1, 
		     colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  par(mar=c(4, 1.5*cexa+2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
  hap<-c(ind*2-1,ind*2)
  G=length(x)
  kLL=nrow(t.Mu)
  glocs=rep(x,100)
  xlim=range(x)
  x=c(x,rev(x))
  ylim=c(0,2)
  plot(x=range(x),y=c(0,2),xlim=xlim,ylim=ylim,t='n',axes=FALSE,ylab=ylab,main=mlab,xlab=xlab)
  alpha.Mu<-t.Mu^pow
  alpha.Mu<-t(t(alpha.Mu)/apply(alpha.Mu,2,max))
  upper=lower=rep(0,G)
  for (i in 1:A) 
    for (ll in 1:kLL) 
    {
      tmprgb=rgb(t(col2rgb(colvec[i])/255),alpha=alpha.Mu[ll,i]) 
      upper=lower+probs[[ch]][[hap[1]]][,(i-1)*kLL+ll]+probs[[ch]][[hap[2]]][,(i-1)*kLL+ll]
      polygon(x=x,y=c(lower,rev(upper)),col=tmprgb,border=NA)
      lower=upper
    }
  axis(2)
}

plot_Mu<-function(t.Mu, t.alpha, t.NL, MODE="scaled", showgradient=FALSE, beside=TRUE, ord=TRUE, pow=1, cexa=1.5, 
		  shiftl=ifelse(showgradient,0,max(sapply(rownames(t.Mu),nchar))/2*cexa), shiftt=ifelse(!showgradient,cexa,cexa*2),
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
  A=ncol(t.Mu)
  par(mar=c(3,1+shiftl,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
  t.Mu[t.Mu<tol]<-tol
  if (ord)
  {
    vord<-NULL
    for (l in 1:A)
    {
      maxforthis<-which(apply(t.Mu,1,which.max)==l)
      vord<-c(vord,maxforthis[sort(t.Mu[maxforthis,l],index=TRUE)$ix])
    }
    t.Mu<-as.matrix(t.Mu[vord,])
  }
  if (!beside)
  {
    #t.Mu=t.Mu[rowMeans(t.Mu)>=quantile(rowMeans(t.Mu),cutoff),] # only those above cutoff on the quantiles on average
    tmp=rep(FALSE,nrow(t.Mu));for (l in 1:A) tmp=(tmp | (t.Mu[,l]>=quantile(t.Mu[,l],cutoff)));t.Mu=t.Mu[tmp,] # only those above cutoff on the quantiles in some anc
  }
  if (!showgradient & !beside)
  {
    barplot(t(t.Mu),space=FALSE,col=colvec,horiz=TRUE,las=TRUE,cex.axis=cexa,cex.names=cexa)
    for (i in 1:A)
      text(x=(i)*0.2*max(rowSums(t.Mu)), y=0.5+1.2*(nrow(t.Mu)+1), round(t.alpha[i],3), cex=cexa, col=colvec[i])
  }
  if (!showgradient & beside)
  {
    if (ord) par(mfrow=c(1,A))
    par(mar=c(3,1+shiftl,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
    if (!ord) 
    {
      nf <- layout(matrix(c(1:(A+1)),ncol=A+1), widths=c(1,rep(2,A)), TRUE);#layout.show(nf)
      plot(c(-cexa*2,ncol(t.Mu)),c(1,nrow(t.Mu)+shiftt),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
      text(x=cexa*2,pos=2,y=0.5+1.2*(1:nrow(t.Mu)),rownames(t.Mu),cex=1.75*cexa)
    }
    if (ord)
    {
      allord=apply(t.Mu,2,order)
      ordMu=NULL;for (a in 1:A) ordMu[[a]]=t.Mu[allord[,a],a] 
    } else 
    {
      ordMu=NULL;for (a in 1:A) ordMu[[a]]=t.Mu[,a] 
    }
    if (cutoff > 0)
    {
      for (a in 1:A)
	ordMu[[a]]=ordMu[[a]][ordMu[[a]]>=quantile(ordMu[[a]],cutoff)] # only show those that are above cutoff on the quantiles
    }
    xmax=max(unlist(ordMu))
    for (a in 1:A) 
    {
      tmpMu=c(ordMu[[a]],NaN);names(tmpMu)[length(ordMu[[a]])+1]=round(t.alpha[a],3)
      if (ord) y=barplot(tmpMu,horiz=TRUE,las=1,col=colvec[a],xlim=c(0,xmax),ylim=c(0,1+length(tmpMu)),cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa*2)
      if (!ord)
	y=barplot(tmpMu,horiz=TRUE,las=1,col=colvec[a],xlim=c(0,xmax),cex.names=cexa,cex.axis=cexa,main="",ylim=c(0,1+length(tmpMu)),
		  cex.main=cexa*2,names.arg=rep("",length(tmpMu)))
    }
  }
  if (showgradient) # overrides beside
  {
    alpha.Mu<-t.Mu^pow
    alpha.Mu<-t(t(alpha.Mu)/apply(alpha.Mu,2,max))
    plot(c(-shiftl,ncol(t.Mu)),c(0.5,nrow(t.Mu)+0.5),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
    mtext(side=3, at=1:A-0.5, round(t.alpha,3), cex=cexa)
    for (l in 1:ncol(t.Mu))
      for (ll in 1:nrow(t.Mu))
      {
	tmprgb=rgb(t(col2rgb(colvec[l])/255),alpha=alpha.Mu[ll,l])
	polygon(x=c(l,l+1,l+1,l)-1,y=c(ll,ll,ll+1,ll+1)-0.5,col=tmprgb,border=NA)
      }
    text(x=0,pos=2,y=(1:nrow(t.Mu)),rownames(t.Mu),cex=0.75*cexa)
  }
  if (showgradient) 
    return(NULL)
  if (ord) 
    if (!beside)
      return(t.Mu) # return as re-ordered version is useful for dipplot_Mu and happlot_Mu
  if (beside)
    return(ordMu) # return as re-ordered version is useful for dipplot_Mu and happlot_Mu
}
plot_Fst<-function(t.Fst, ord=TRUE, cexa=1, shiftl=cexa, shiftt=cexa, cutoff=nrow(t.Fst), reverse=TRUE,
		   colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  t.kLL=nrow(t.Fst)
  if ((!ord) & (cutoff!=t.kLL))
  {
    warning("########## showing all as re-ordering not allowed ##########", immediate.=TRUE)
    cutoff=t.kLL
  }
  if (cutoff>t.kLL) cutoff=t.kLL
  A=ncol(t.Fst)
  par(mar=c(6,4+shiftl,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
  if (ord) par(mfrow=c(1,A))
  if (!ord) 
  {
    par(mar=c(6,2,1+shiftt,1),bty='n', cex.axis=cexa, cex.lab=cexa, cex.main=cexa)
    nf <- layout(matrix(c(1:(A+1)),ncol=A+1), widths=c(1,rep(2,A)), TRUE);#layout.show(nf)
    plot(c(-cexa,ncol(t.Fst)),c(0.5,nrow(t.Fst)+0.5),t='n',yaxt='n',xaxt='n',main="",xlab="",ylab="")
    text(x=2,pos=2,y=(1:nrow(t.Fst)),rownames(t.Fst),cex=cexa)
  }
  if (ord)
  {
    allord=apply(t.Fst,2,order,decreasing=TRUE)
    ordFst=NULL;for (a in 1:A) ordFst[[a]]=t.Fst[allord[,a],a] 
  } else 
  {
    ordFst=NULL;for (a in 1:A) ordFst[[a]]=t.Fst[,a] 
  }
  rangeFst=NULL
  for (a in 1:A) {
    ordFst[[a]]=tail(ordFst[[a]],cutoff) # only show those that are below cutoff place
    if (reverse) 
      rangeFst[[a]]=range(1-ordFst[[a]])
    if (!reverse) 
      rangeFst[[a]]=range(ordFst[[a]])
  }
  for (a in 1:A) 
  {
    if (ord) {
      if (reverse)
	plot(c(0,1),c(0,cutoff+1),t='n',yaxt='n',ylab="",xlab="1-Fst",cex.main=cexa,main="",cex.axis=cexa,xlim=c(rangeFst[[a]][1],rangeFst[[a]][2]))
      if (!reverse)
	plot(c(0,1),c(0,cutoff+1),t='n',yaxt='n',ylab="",xlab="Fst",cex.main=cexa,main="",cex.axis=cexa,xlim=c(0,rangeFst[[a]][2]))
      if (reverse) {
	tmp=1-ordFst[[a]]-rangeFst[[a]][1]
	y=barplot(tmp,horiz=TRUE,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa,add=TRUE,offset=rangeFst[[a]][1])
      }
      if (!reverse) 
	y=barplot(ordFst[[a]],horiz=TRUE,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",cex.main=cexa,add=TRUE)
    }
    if (!ord & reverse)
      y=barplot(1-ordFst[[a]]-rangeFst[[a]][1],horiz=TRUE,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",xlab="1-Fst",
		cex.main=cexa,names.arg=rep("",length(ordFst[[a]])),offset=rangeFst[[a]][1])
    if (!ord & !reverse)
      y=barplot(ordFst[[a]],horiz=TRUE,las=1,col=colvec[a],cex.names=cexa,cex.axis=cexa,main="",xlab="Fst",
		cex.main=cexa,names.arg=rep("",length(ordFst[[a]])))
  }
  if (ord) 
    return(ordFst)
}

# localanc is w.r.t. final phasing. The below undoes this to compare with input phasing (necessary for Fst calculations)
phase_localanc=function(t.localanc,t.flips) 
{
  nchrno=length(t.localanc)
  A=dim(t.localanc[[1]])[1]
  NUMA=dim(t.localanc[[1]])[2]
  for (ch in 1:nchrno)
    for (ind in 1:(NUMA/2))
    {
      haps=c(ind*2-1,ind*2)
      ans=t.localanc[[ch]][,haps,]
      for (k in haps)
      {
	otherhap=ifelse((k%%2)==0,1,2) # if even look at previous; if odd look at next
	kf=t.flips[[ind]][[ch]] # FALSE and TRUE values
	for (l in 1:A)
	  t.localanc[[ch]][l,k,kf]=ans[l,otherhap,kf]
      }
    }
  return(t.localanc)
}

# function to plot the local ancestry of each target admixed genome along each chromosome
plot_localanc=function(t.chrnos, t.g.loc, t.localanc, t.g.true_anc=NULL,cexa=2,pow=1,y.lab="expected",MODE="BAR",NCHR=2,PAUSE=TRUE,t.Mu=NULL,
		       colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")) { 
  G=sapply(t.localanc, function(x) dim(x)[3])
  NUMA=dim(t.localanc[[1]])[2]
  A=dim(t.localanc[[1]])[1]
  nchrno=length(t.chrnos)
  NUMI=NUMA/2
  if (NCHR==2 & NUMA==1)
  {
    warning("########## changing NCHR to HAP as only on hap in target ##########",immediate.=TRUE)
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
	  perc[ch,k]<-round(100*mean(apply(matrix(t.localanc[[ch]][,k,],A),2,which.max)==
				     apply(matrix(t.g.true_anc[[ch]][,k,],A),2,which.max)),2)
	  mse[ch,k]<-mean((t.localanc[[ch]][,k,]-t.g.true_anc[[ch]][,k,])^2)
	  r2[ch,k]=hap_fr2_chr_k(t.localanc, t.g.true_anc, ch, k)
	  cat("r2: ",r2[ch,k],"  ", sep="")
	  cat("correct: ", perc[ch,k], "%   ", sep="")
	  cat("MSE:", mse[ch,k], "  ", sep="")
	  cat("\n")
	  par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	  if (MODE=="LINE")
	  {
	    plot(range(t.g.loc[[ch]])*1e-6,c(0,1),axes=FALSE,t='n',ylab="truth",main=paste("Haplotype", k), xlab=paste("Mb Position on Chromosome",t.chrnos[ch]))
	    for (i in 1:A) lines(t.g.loc[[ch]]*1e-6, t.g.true_anc[[ch]][i,k,], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa);	
	    axis(2)
	  }
	  if (MODE=="BAR")
	    happlot(ch,k,t.g.loc[[ch]]*1e-6,A,t.g.true_anc,ylab="truth",cexa=cexa)
	}
	par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	if (MODE=="LINE")
	{
	  plot(range(t.g.loc[[ch]])*1e-6,c(0,1),axes=FALSE,t='n',ylab=y.lab,main=paste("Haplotype", k), xlab=paste("Mb Position on Chromosome",t.chrnos[ch]))
	  for (i in 1:A) lines(t.g.loc[[ch]]*1e-6, t.localanc[[ch]][i,k,], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa);	
	  axis(2)
	}
	if (MODE=="BAR")
	  happlot(ch,k,t.g.loc[[ch]]*1e-6,A,t.localanc,xlab=paste("Mb Position on Chromosome",t.chrnos[ch]),ylab=y.lab,cexa=cexa)
	mp<-axTicks(1,axp=round(c(min(t.g.loc[[ch]])*1e-6,max(t.g.loc[[ch]])*1e-6,5)))
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
	  for (i in 1:A)
	    perc[ch,ind]<-perc[ch,ind]+mean(as.integer(t.localanc[[ch]][i,hap[1],]+t.localanc[[ch]][i,hap[2],]+0.5)==
					    as.integer(t.g.true_anc[[ch]][i,hap[1],]+t.g.true_anc[[ch]][i,hap[2],]+0.5)) # push to int required here due to shift to the grid
	  perc[ch,ind]<-round(100*perc[ch,ind]/A,2)
	  mse[ch,ind]<-mean((t.localanc[[ch]][,hap[1],]+t.localanc[[ch]][,hap[2],]-t.g.true_anc[[ch]][,hap[1],]-t.g.true_anc[[ch]][,hap[2],])^2)
	  r2[ch,ind]=dip_fr2_chr_ind(t.localanc, t.g.true_anc, ch, ind)
	  cat("r2: ",r2[ch,ind],"  ", sep="")
	  cat("% correct: ", perc[ch,ind], "  ", sep="")
	  cat("MSE: ", mse[ch,ind], "  ", sep="")
	  cat("\n")
	  par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	  if (MODE=="LINE")
	  {
	    plot(range(t.g.loc[[ch]])*1e-6,c(0,2),axes=FALSE,t='n',ylab="truth",main=paste("Individual", ind),xlab=paste("Mb Position on Chromosome",t.chrnos[ch]))
	    for (i in 1:A) lines(t.g.loc[[ch]]*1e-6, t.g.true_anc[[ch]][i,hap[1],]+t.g.true_anc[[ch]][i,hap[2],], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa)
	    axis(2)
	  }
	  if (MODE=="BAR")
	    dipplot(ch,ind,t.g.loc[[ch]]*1e-6,A,t.g.true_anc,xlab=paste("Mb Position on Chromosome",t.chrnos[ch]),
		    ylab="truth",cexa=cexa)
	}
	if (MODE=="LINE")
	{
	  plot(range(t.g.loc[[ch]])*1e-6,c(0,2),axes=FALSE,t='n',ylab=y.lab,main=paste("Individual", ind),xlab=paste("Mb Position on Chromosome",t.chrnos[ch]))
	  for (i in 1:A) lines(t.g.loc[[ch]]*1e-6, t.localanc[[ch]][i,hap[1],]+t.localanc[[ch]][i,hap[2],], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa)
	  axis(2)
	}
	if (MODE=="BAR")
	  dipplot(ch,ind,t.g.loc[[ch]]*1e-6,A,t.localanc,xlab=paste("Mb Position on Chromosome",t.chrnos[ch]),ylab=y.lab,cexa=cexa)
	mp<-axTicks(1,axp=round(c(min(t.g.loc[[ch]])*1e-6,max(t.g.loc[[ch]])*1e-6,5)))
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
       ylab="mean ancestry",xlab=paste("Mb Position on Chromosome",chrnos[ch]),main="")  
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

# function to plot most useful figures to PDF. Note defaults are 
plot_all_mosaic=function(pathout,target,EM,PHASE,t.GpcM=GpcM,t.all_Fst=all_Fst,t.A=A,t.NUMA=NUMA,
			 t.Mu=Mu,t.chrnos=chrnos,t.alpha=alpha,t.NL=NL,t.acoancs=acoancs,t.dr=dr,t.logfile=logfile, 
			 t.localanc, t.g.loc=g.loc, g.true_anc=NULL){
  if (!dir.exists(pathout))
    dir.create(file.path(pathout))
  nchrno=length(t.chrnos)
  targetdetails=paste0(target, "_", t.A, "way_", t.NUMA, "_", paste(t.chrnos[c(1,nchrno)],collapse="-"),
		       "_",sum(t.NL),"_",t.GpcM)
  pdf(file=file.path(pathout,paste0(targetdetails,"_Mu.pdf")), width=12, height=7)
  ord.Mu=plot_Mu(t.Mu,t.alpha,t.NL,cexa=1.5,beside=TRUE,shiftl=5,shiftt=2,cutoff=0,ord=TRUE)
  dev.off()

  pdf(file=file.path(pathout,paste0(targetdetails,"_karyograms.pdf")), width=12, height=12)
  for (ind in 1:(t.NUMA/2)) 
    karyogram(t.A, t.chrnos, t.localanc, t.g.loc, t.GpcM, ind, dist="genetic", g.true_anc=g.true_anc) # one per PDF page
  dev.off()

  # note that it takes a while to calculate frequencies, etc
  if (!is.null(t.all_Fst))
    if (all(!is.nan(t.all_Fst$panels))) {
      pdf(file=file.path(pathout,paste0(targetdetails,"_Fst.pdf")), width=21, height=28)
      ord.Fst=plot_Fst(t.all_Fst$panels,cexa=2,ord=TRUE, shiftl=6, cutoff=10)
      dev.off()
    }

  # dimensions of plots
  d1=switch(t.A,NaN,1,2,2,3,3) 
  d2=switch(t.A,NaN,3,3,5,5,7) 
  pdf(file=file.path(pathout,paste0(targetdetails,"_acoanc.pdf")), width=5*d2,height=5*d1)
  this_acoplots=plot_coanccurves(t.acoancs,t.dr,lwd=4,cexa=2,verbose=FALSE,axisall=FALSE,samedates=FALSE,asym=FALSE)
  dev.off()
  if (EM | PHASE) {
    EMlog=extract_log(t.logfile)
    pdf(file=file.path(pathout,paste0(targetdetails,"_EMlog.pdf")), width=14,height=14)
    plot_loglike(EMlog)
    dev.off()
  }
}

# create a map showing closest populations to mixing groups based on MOSAIC output
# sources should be a matrix NxA where N is #sources and A is #ancestries
# geolocs is matrix with each row a name/Longitude/Latitude trio. Must contain same population names as sources matrix.
# cexa scaling factor for pie chart sizes
# byFst specifies if we're using Fst or copying matrix values to plot
plot_admix_map=function(sources, geolocs, cexa=1, byFst=TRUE) { 
  sourcelocs=matrix(NaN,nrow(sources),2);
  for (i in 1:nrow(sources)) {
    k=which(geolocs[,1]==(rownames(sources)[i]))
    sourcelocs[i,1]=as.numeric(geolocs[k,2])
    sourcelocs[i,2]=as.numeric(geolocs[k,3])
  }
  if (byFst) {
    sources[sources<0]=0 # negative estimates are really zeros
    sources=1/sources # get 1-Fst
    sources=t(t(sources)/apply(sources,2,max)) # rescale as a ratio to closest population to each ancestry
  }
  sources=cexa*sources*ncol(sources)/sum(sources) # rescale to sum to A
  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")
  drawpie=function(center,radius,probs,n=50,colours=colvec[1:length(probs)],bord=FALSE,...)
  {
    x <- c(0,cumsum(probs)/sum(probs))
    dx <- diff(x)
    np <- length(probs)
    for (i in 1:np)
    {
      t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
      xc <- center[1] + c(cos(t2p), 0) * radius
      yc <- center[2] + c(sin(t2p), 0) * radius
      polygon(xc, yc, border = bord, col = colours[i],...)
    }
  }  

  ii=which(geolocs[,1]==target)
  par(bg="lightblue",mar=c(0,0,0,0))
  plot(sourcelocs,t='n',xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(-125,155),ylim=c(-46,63))
  map(col="lightgray",add=TRUE,fill=TRUE,border="darkgray")
  showsource=rep(FALSE,nrow(sources))
  tmax=0 # will be as 1.5 times big as biggest source pie chart
  for (i in 1:nrow(sources)) 
    if (any(tapply(1:A,1:A,function(l) (rank(-sources[,l])[i])<=5))) {  # take top 5 in each side
      tmax=max(tmax,sum(sources[i,]))
      drawpie(sourcelocs[i,],sum(sources[i,])*10,sources[i,],col=colvec[1:A])
      showsource[i]=TRUE
    }
  # add in target pie chart with border
  rad=1.5*tmax*10
  drawpie(as.numeric(geolocs[ii,2:3]),rad,Reduce("+",alpha)/NUMA*2,col=colvec[1:A],bord=TRUE,lwd=2)
  for (i in 1:nrow(sources)) {
    if (showsource[i]) {
      xy=c(geolocs[ii,2]-sourcelocs[i,1],geolocs[ii,3]-sourcelocs[i,2])
      theta=atan2(xy[2],xy[1]) # tan(theta)=y/x
      endp=c(geolocs[ii,2]-rad*cos(theta),geolocs[ii,3]-rad*sin(theta))
      arrows(sourcelocs[i,1],sourcelocs[i,2],endp[1],endp[2],lwd=2,length=0.1)
    }
  }
}

# function to plot local ancestry karyogram for one individual based on estimated local ancestry 
karyogram=function(t.A, chrnos, localanc, g.loc, GpcM, ind, dist="genetic", cexa=1.5, m.lab=paste("individual", ind), g.true_anc=NULL) {
  colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999")
  hap=c(ind*2-1,ind*2)
  if (is.null(g.true_anc))
    par(mar=c(4, 5.2, cexa*2, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa, yaxt='n')
  if (!is.null(g.true_anc)) {
    par(mar=c(4, 5.2, cexa*2, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa, yaxt='n',mfrow=c(1,2))
    m.lab=c(paste0("E[",m.lab,"]"),paste0("truth[",m.lab,"]"))
  }
  if (dist=="physical")
    plot(range(sapply(g.loc,range))*1e-6, c(min(chrnos)-0.5,max(chrnos)+0.5), ,t='n',xlab="position (Mb)",ylab="chromosome", main=m.lab[1])
  if (dist=="genetic")
    plot(c(0,max(sapply(g.loc,length))/GpcM), c(min(chrnos)-0.5,max(chrnos)+0.5), ,t='n',xlab="genetic position (cM)",ylab="chromosome", main=m.lab[1])
  mtext(chrnos, 2, at=chrnos,cex=cexa,las=1)
  for (ch in 1:length(chrnos)) {
    G=length(g.loc[[ch]])
    if (dist=="physical") {
      x=c(g.loc[[ch]],rev(g.loc[[ch]]))*1e-6
      upper=lower=rep(0,G)
    }
    if (dist=="genetic") {
      x=c(0:(G-1),(G-1):0)/GpcM
      upper=lower=rep(0,G)
    }
    for (i in 1:t.A) 
    {
      upper=lower+localanc[[ch]][i,hap[1],]+localanc[[ch]][i,hap[2],]
      polygon(x=x,y=chrnos[ch]-0.5+0.4*c(lower,rev(upper)),col=colvec[i],border=NA)
      lower=upper
    }
  }
  if (!is.null(g.true_anc)) {
    if (dist=="physical")
      plot(range(sapply(g.loc,range))*1e-6, c(min(chrnos)-0.5,max(chrnos)+0.5), ,t='n',xlab="position (Mb)",ylab="chromosome", main=m.lab[2])
    if (dist=="genetic")
      plot(c(0,max(sapply(g.loc,length))/GpcM), c(min(chrnos)-0.5,max(chrnos)+0.5), ,t='n',xlab="genetic position (cM)",ylab="chromosome", main=m.lab[2])
    mtext(chrnos, 2, at=chrnos,cex=cexa,las=1)
    for (ch in 1:length(chrnos)) {
      G=length(g.loc[[ch]])
      if (dist=="physical") {
        x=c(g.loc[[ch]],rev(g.loc[[ch]]))*1e-6
        upper=lower=rep(0,G)
      }
      if (dist=="genetic") {
        x=c(0:(G-1),(G-1):0)/GpcM
        upper=lower=rep(0,G)
      }
      for (i in 1:t.A) 
      {
        upper=lower+g.true_anc[[ch]][i,hap[1],]+g.true_anc[[ch]][i,hap[2],]
        polygon(x=x,y=chrnos[ch]-0.5+0.4*c(lower,rev(upper)),col=colvec[i],border=NA)
        lower=upper
      }
    }
  }
}

