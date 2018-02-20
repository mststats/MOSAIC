source("plot_funcs.R")
source("calc_r2.R")
if (!exists("G")) G=sapply(localanc, function(x) dim(x)[3])
NN=sum(NL)
NUMI=NUMA/2
if (!exists("cexa")) cexa=2
if (!exists("pow")) pow=1
if (!exists("y.lab")) y.lab="expected"
if (!exists("MODE")) MODE="BAR" # or GRAD or LINE
if (!exists("NCHR")) NCHR=2
if (!exists("PAUSE")) PAUSE=T
if (!exists("PNG")) PNG=F
#source("localanc.R")
if (NCHR==2 & NUMA==1)
{
  warning("changing NCHR to HAP as only on hap in target",immediate.=T)
  NCHR=1
}
#if (target=="simulated") # calculate the true number of gridpoint anc switches => rate of anc switching => #gens since admixture
#{
  #obs.rate=sapply(g.true_anc, function(x) {m<-(apply(abs(apply(x[,k,],1,diff)),1,sum)>tol);mean(m)}) # switch rate / chr
  #obs.lambda=-log(1-obs.rate)/dr # lambda / chr
  #obs.lambda=-log(1-sum(obs.sum)/sum(G-nchrno))/dr
  # SM way
  #obs.sum=sapply(g.true_anc, function(x) {m<-(apply(abs(apply(x[,k,],1,diff)),1,sum)>tol);sum(m)})
  #obs.lambda=sum(obs.sum)*L/(L-1)/2/(sum(G)*dr)
  #obs.lambda=mean(sapply(g.true_anc, function(x) {m<-(apply(abs(apply(x[,k,],1,diff)),1,sum)>tol);mean(m)}))*L/2/2/dr
  #cat("true #gens since admixture =", mean(obs.lambda), "\n")
#}
if (NCHR==1)
{
  expected_r2<-matrix(NaN,nchrno,NUMA)
  if (target=="simulated") 
  {
    perc<-matrix(NaN,nchrno,NUMA)
    mse<-matrix(NaN,nchrno,NUMA)
    r2<-matrix(NaN,nchrno,NUMA)
  }
  for (ch in 1:nchrno)
    for (k in 1:NUMA)
    {
      ind=as.integer((k+1)/2)
      if (PNG)
	png(file=paste0("PLOTS/Localanc/", target, "_HAP_", NN, "_", k,".",chrnos[ch],".png"),width=1920,height=960)
      expected_r2[ch,k]=hap_expected_fr2_chr_k(localanc, ch, k)
      cat("E[r2]: ", expected_r2[ch,k], "  ", sep="")
      if (target!="simulated") par(mfrow=c(1,1),mar=c(4, 5.2, 1, 0)+0.1)
      if (target=="simulated") par(mfrow=c(2,1),mar=c(4, 5.2, 1, 0)+0.1)
      if (target=="simulated")
      {
	# how many maximal ancestry predictions are the same as the truth
	perc[ch,k]<-round(100*mean(apply(matrix(localanc[[ch]][,k,],L),2,which.max)==
				   apply(matrix(g.true_anc[[ch]][,k,],L),2,which.max)),2)
	mse[ch,k]<-mean((localanc[[ch]][,k,]-g.true_anc[[ch]][,k,])^2)
	r2[ch,k]=hap_fr2_chr_k(localanc, g.true_anc, ch, k)
	cat("r2: ",r2[ch,k],"  ", sep="")
	cat("correct: ", perc[ch,k], "%   ", sep="")
	cat("MSE:", mse[ch,k], "  ", sep="")
	cat("\n")
        par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
        if (MODE=="LINE")
          {
	    plot(range(g.loc[[ch]])*1e-6,c(0,1),axes=F,t='n',ylab="truth",main=paste("Haplotype", k), xlab=paste("Position on Chromosome",chrnos[ch]))
            for (i in 1:L) lines(g.loc[[ch]]*1e-6, g.true_anc[[ch]][i,k,], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa);	
	    axis(2)
	  }
        if (MODE=="BAR" | MODE=="GRAD")
	  happlot(ch,k,g.loc[[ch]]*1e-6,g.true_anc[[ch]][,k,],"truth",cexa=cexa)
      }
      par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
      if (MODE=="LINE")
        {
	  plot(range(g.loc[[ch]])*1e-6,c(0,1),axes=F,t='n',ylab=y.lab,main=paste("Haplotype", k), xlab=paste("Position on Chromosome",chrnos[ch]))
          for (i in 1:L) lines(g.loc[[ch]]*1e-6, localanc[[ch]][i,k,], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa);	
	  axis(2)
	}
      if (MODE=="BAR")
	happlot(ch,k,g.loc[[ch]]*1e-6,localanc[[ch]][,k,],ylab=y.lab,cexa=cexa)
      if (MODE=="GRAD")
	happlot_Mu(ch,k,g.loc[[ch]]*1e-6,gfbs[[ch]],ylab=y.lab,t.Mu=Mu,pow=pow,cexa=cexa)
      mp<-axTicks(1,round(axp=c(min(g.loc[[ch]])*1e-6,max(g.loc[[ch]])*1e-6,5)))
      axis(1,at=mp,labels=signif(mp,3))
      if (RPE>0) 
	if (length(phase.error.locs[[as.integer((k+1)/2)]][[ch]])>0)
	  mtext("|",side=3,at=g.loc[[ch]][phase.error.locs[[as.integer((k+1)/2)]][[ch]]]*1e-6,cex=0.5,col=8,line=-0.5)
      #if (chrnos[ch]==6) abline(v=c(28510120,33480577),col=4,lwd=cexa)
      if (PNG) dev.off()
      if (PAUSE) readline()
    }
}
if (NCHR==2)
{
  expected_r2<-matrix(NaN,nchrno,NUMI)
  if (target=="simulated") 
  {
    perc<-matrix(NaN,nchrno,NUMI)
    mse<-matrix(NaN,nchrno,NUMI)
    r2<-matrix(NaN,nchrno,NUMI)
  }
  for (ch in 1:nchrno)
    for (ind in 1:(NUMI))
    {
      hap<-c(ind*2-1,ind*2)
      if (PNG)
	png(file=paste0("PLOTS/Localanc/", target, "_DIP_", NN, "_", ind,".",chrnos[ch],".png"),width=1920,height=960)
      expected_r2[ch,ind]=dip_expected_fr2_chr_ind(localanc, ch, ind)
      cat("E[r2]: ", expected_r2[ch,ind], "  ", sep="")
      if (target!="simulated") par(mfrow=c(1,1),mar=c(4, 5.2, 1, 0)+0.1)
      if (target=="simulated") par(mfrow=c(2,1),mar=c(4, 5.2, 1, 0)+0.1)
      if (target=="simulated")
      {
	perc[ch,ind]<-0 # how many maximal ancestry counts 0,1,2 are the same as the truth
	for (i in 1:L)
	  perc[ch,ind]<-perc[ch,ind]+mean(as.integer(localanc[[ch]][i,hap[1],]+localanc[[ch]][i,hap[2],]+0.5)==
			                  as.integer(g.true_anc[[ch]][i,hap[1],]+g.true_anc[[ch]][i,hap[2],]+0.5)) # push to int required here due to shift to the grid
	perc[ch,ind]<-round(100*perc[ch,ind]/L,2)
	mse[ch,ind]<-mean((localanc[[ch]][,hap[1],]+localanc[[ch]][,hap[2],]-g.true_anc[[ch]][,hap[1],]-g.true_anc[[ch]][,hap[2],])^2)
	r2[ch,ind]=dip_fr2_chr_ind(localanc, g.true_anc, ch, ind)
	cat("r2: ",r2[ch,ind],"  ", sep="")
	cat("% correct: ", perc[ch,ind], "  ", sep="")
	cat("MSE: ", mse[ch,ind], "  ", sep="")
	cat("\n")
        par(mar=c(4, 5.2, cexa, 0), cex.main=cexa, cex.axis=cexa, cex.lab=cexa)
	if (MODE=="LINE")
	{
	  plot(range(g.loc[[ch]])*1e-6,c(0,2),axes=F,t='n',ylab="truth",main=paste("Individual", ind),xlab=paste("Position on Chromosome",chrnos[ch]))
          for (i in 1:L) lines(g.loc[[ch]]*1e-6, g.true_anc[[ch]][i,hap[1],]+g.true_anc[[ch]][i,hap[2],], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa)
          axis(2)
	}
        if (MODE=="BAR" | MODE=="GRAD")
	  dipplot(ch,ind,g.loc[[ch]]*1e-6,g.true_anc[[ch]][,hap[1],]+g.true_anc[[ch]][,hap[2],],ylab="truth",cexa=cexa)
      }
      if (MODE=="LINE")
      {
	plot(range(g.loc[[ch]])*1e-6,c(0,2),axes=F,t='n',ylab=y.lab,main=paste("Individual", ind),xlab=paste("Position on Chromosome",chrnos[ch]))
	for (i in 1:L) lines(g.loc[[ch]]*1e-6, localanc[[ch]][i,hap[1],]+localanc[[ch]][i,hap[2],], t='l', col=rgb(t(col2rgb(colvec[i])/255),alpha=0.5), lwd=cexa)
        axis(2)
      }
      if (MODE=="BAR")
	dipplot(ch,ind,g.loc[[ch]]*1e-6,localanc[[ch]][,hap[1],]+localanc[[ch]][,hap[2],],,ylab=y.lab,cexa=cexa)
      if (MODE=="GRAD")
	xy<-dipplot_Mu(ch,ind,g.loc[[ch]]*1e-6,gfbs[[ch]],ylab=y.lab,t.Mu=Mu,pow=pow,cexa=cexa)
      mp<-axTicks(1,round(axp=c(min(g.loc[[ch]])*1e-6,max(g.loc[[ch]])*1e-6,5)))
      axis(1,at=mp,labels=signif(mp,3))
      if (RPE>0) 
	if (length(phase.error.locs[[ind]][[ch]])>0)
	  mtext("|",side=3,at=g.loc[[ch]][phase.error.locs[[ind]][[ch]]]*1e-6,cex=cexa/2,col=8,line=-0.5)
      #if (chrnos[ch]==6) abline(v=c(28510120,33480577),col=4,lwd=cexa)
      if (PNG) dev.off()
      if (PAUSE) readline()
    }
}
