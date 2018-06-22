# script to plot the local ancestry of each target admixed genome along each chromosome
plot_localanc=function(t.chrnos, t.g.loc, t.localanc, t.g.true_anc=NULL,cexa=2,pow=1,y.lab="expected",MODE="BAR",NCHR=2,PAUSE=T,t.Mu=NULL,t.gfbs=NULL) {
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
