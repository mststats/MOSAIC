a=1
if (HapMixMODE=="HAP")
{
  if (showunphased)
  cat("HapMix r^2 = ", hap_fr2(hap_hmixlocalanc,true_anc), 
      " unphased MOSAIC r^2 = ", hap_fr2(noanc_unphased_localanc_gs,g.true_anc_gs),
      " phased MOSAIC r^2 = ", hap_fr2(localanc_gs,g.true_anc_gs), sep="", "\n")
  if (!showunphased)
  cat("HapMix r^2 = ", hap_fr2(hap_hmixlocalanc,true_anc), 
      " MOSAIC r^2 = ", hap_fr2(localanc_gs,g.true_anc_gs), sep="", "\n")
}
if (HapMixMODE=="DIP" | HapMixMODE=="BOTH") 
{
  if (showunphased)
  cat("HapMix r^2 = ", dip_fr2(hap_hmixlocalanc,true_anc), 
      " unphased MOSAIC r^2 = ", dip_fr2(noanc_unphased_localanc_gs,g.true_anc_gs), 
      " phased MOSAIC r^2 = ", dip_fr2(localanc_gs,g.true_anc_gs), sep="", "\n")
  if (!showunphased)
  cat("HapMix r^2 = ", dip_fr2(hap_hmixlocalanc,true_anc), 
      " MOSAIC r^2 = ", dip_fr2(localanc_gs,g.true_anc_gs), sep="", "\n")
}
readline("Hit return to see plots")
for (ch in 1:nchrno)
{
  if (HapMixMODE=="HAP")
  for (k in 1:NUMA) 
  {
    if (PNG) png(file=paste0("PLOTS/", NAPOP, k, "_G:", G, "_EM:",EM,"_",k,"_",".png"),width=480*2,height=480)
    if (PLOTMODE=="LINE")
    {
      par(mar=c(5,5,4,2),mfrow=c(1,1))
      plot(g.loc[[ch]][gs[[ch]]],g.true_anc_gs[[ch]][a,k,],col=3,lwd=cexa*2,t='l',
	   ylab=paste("E[anc",a,"alleles]"), xlab=paste("Position on Chromosome",chrnos[ch]),
	   main=paste("haplotype", k),cex.axis=cexa,cex.lab=cexa,cex.main=cexa)
      #lines(hmixloc[[ch]],true_anc[[ch]][a,k,],col=4,lwd=cexa*2,lty=2) # sanity check; looks ok
      lines(hmixloc[[ch]],hap_hmixlocalanc[[ch]][a,k,],col=2,lwd=cexa*2)
      if (showunphased)
        lines(g.loc[[ch]][gs[[ch]]],noanc_unphased_localanc_gs[[ch]][a,k,],lwd=cexa*2,col=4)
      lines(g.loc[[ch]][gs[[ch]]],localanc_gs[[ch]][a,k,],lwd=cexa*2)
      if (!showunphased)
        legend("topleft",c("Truth","HapMix","MOSAIC"),col=3:1,lty=c(1,1),lwd=cexa*2,cex=cexa)
      if (showunphased)
        legend("topleft",c("Truth","HapMix","unphased MOSAIC","MOSAIC"),col=c(3,2,4,1),lty=c(1,1),lwd=cexa*2,cex=cexa)
    }
    if (PLOTMODE=="BAR")
    {
      if (!showunphased)
        par(mfrow=c(3,1))
      if (showunphased)
        par(mfrow=c(4,1))
      happlot(ch,k,g.loc[[ch]][gs[[ch]]]*1e-6,g.true_anc_gs[[ch]][,k,],ylab="truth",cexa=cexa)
      happlot(ch,k,hmixloc[[ch]]*1e-6,hap_hmixlocalanc[[ch]][,k,],ylab="HapMix",cexa=cexa,mlab="")
      if (showunphased)
        happlot(ch,k,g.loc[[ch]][gs[[ch]]]*1e-6,noanc_unphased_localanc_gs[[ch]][,k,],ylab="unphased MOSAIC",cexa=cexa,mlab="")
      happlot(ch,k,g.loc[[ch]][gs[[ch]]]*1e-6,localanc_gs[[ch]][,k,],ylab="MOSAIC",cexa=cexa,mlab="")
      mp<-axTicks(1,round(axp=c(min(g.loc[[ch]])*1e-6,max(g.loc[[ch]])*1e-6,5)))
      axis(1,at=mp,labels=signif(mp,3))
    }
    if (!showunphased)
      cat("HapMix r^2 = ", hap_fr2_chr_k(hap_hmixlocalanc,true_anc,ch,k), " MOSAIC r^2 = ", hap_fr2_chr_k(localanc_gs,g.true_anc_gs,ch,k), sep="", "\n")
    if (showunphased)
      cat("HapMix r^2 = ", hap_fr2_chr_k(hap_hmixlocalanc,true_anc,ch,k), 
	  " unphased MOSAIC r^2 = ", hap_fr2_chr_k(noanc_unphased_localanc_gs,g.true_anc_gs,ch,k),
	  " phased MOSAIC r^2 = ", hap_fr2_chr_k(localanc_gs,g.true_anc_gs,ch,k), sep="", "\n")
    if (!PNG) readline()
    if(PNG) dev.off()
  }
  if (HapMixMODE=="DIP" | HapMixMODE=="BOTH") 
  for (k in 1:NUMI) 
  {
    if (PNG) png(file=paste0("PLOTS/", NAPOP, k, "_G:", G, "_EM:",EM,"_",k,"_",".png"),width=480*2,height=480)
    haps=c(k*2-1,k*2)
    if (PLOTMODE=="LINE")
    {
      par(mar=c(5,5,4,2),mfrow=c(1,1))
      plot(g.loc[[ch]][gs[[ch]]],g.true_anc_gs[[ch]][a,haps[1],]+g.true_anc_gs[[ch]][a,haps[2],],col=3,lwd=cexa*2,t='l',
	   ylab=paste("E[anc",a,"alleles]"), xlab=paste("Position on Chromosome",chrnos[ch]),main=paste("Individual", k),cex.axis=cexa,cex.lab=cexa,cex.main=cexa)
      #lines(hmixloc[[ch]],true_anc[[ch]][a,haps[1],]+true_anc[[ch]][a,haps[2],],col=3,lwd=cexa*2,lty=2) # sanity check; looks ok
      lines(hmixloc[[ch]],hap_hmixlocalanc[[ch]][a,haps[1],]+hap_hmixlocalanc[[ch]][a,haps[2],],col=2,lwd=cexa*2)
      lines(g.loc[[ch]][gs[[ch]]],localanc_gs[[ch]][a,haps[1],]+localanc_gs[[ch]][a,haps[2],],lwd=cexa*2)
      if (showunphased)
        lines(g.loc[[ch]][gs[[ch]]],noanc_unphased_localanc_gs[[ch]][a,haps[1],]+noanc_unphased_localanc_gs[[ch]][a,haps[2],],lwd=cexa*2,col=4)
      #lines(hmixloc[[ch]],hmixlocalanc[[ch]][a,k,],col=2,lwd=cexa*2) # same as above line only for DIP HapMixMODE
      if (!showunphased)
        legend("topleft",c("Truth","HapMix","MOSAIC"),col=3:1,lty=c(1,1),lwd=cexa*2,cex=cexa)
      if (showunphased)
        legend("topleft",c("Truth","HapMix","unphased MOSAIC","MOSAIC"),col=c(3,2,4,1),lty=c(1,1),lwd=cexa*2,cex=cexa)
    }
    if (PLOTMODE=="BAR")
    {
      if (!showunphased)
        par(mfrow=c(3,1))
      if (showunphased)
        par(mfrow=c(4,1))
      dipplot(ch,k,g.loc[[ch]][gs[[ch]]]*1e-6,g.true_anc_gs[[ch]][,haps[1],]+g.true_anc_gs[[ch]][,haps[2],],ylab="truth",cexa=cexa)
      dipplot(ch,k,hmixloc[[ch]]*1e-6,hap_hmixlocalanc[[ch]][,haps[1],]+hap_hmixlocalanc[[ch]][,haps[2],],ylab="HapMix",cexa=cexa,mlab="")
      #dipplot(ch,k,hmixloc[[ch]]*1e-6,hmixlocalanc[[ch]][,k,],ylab="HapMix",cexa=cexa,mlab="") # same as above line only for DIP HapMixMODE
      if (showunphased)
        dipplot(ch,k,g.loc[[ch]][gs[[ch]]]*1e-6,noanc_unphased_localanc_gs[[ch]][,haps[1],]+noanc_unphased_localanc_gs[[ch]][,haps[2],],ylab="unphased MOSAIC",cexa=cexa,mlab="")
      dipplot(ch,k,g.loc[[ch]][gs[[ch]]]*1e-6,localanc_gs[[ch]][,haps[1],]+localanc_gs[[ch]][,haps[2],],ylab="MOSAIC",cexa=cexa,mlab="")
      mp<-axTicks(1,round(axp=c(min(g.loc[[ch]])*1e-6,max(g.loc[[ch]])*1e-6,5)))
      axis(1,at=mp,labels=signif(mp,3))
    }
    if (!showunphased)
      cat("HapMix r^2 = ", dip_fr2_chr_ind(hap_hmixlocalanc,true_anc,ch,k), " MOSAIC r^2 = ", dip_fr2_chr_ind(localanc_gs,g.true_anc_gs,ch,k), sep="", "\n")
    if (showunphased)
      cat("HapMix r^2 = ", dip_fr2_chr_ind(hap_hmixlocalanc,true_anc,ch,k), 
	  " unphased MOSAIC r^2 = ", dip_fr2_chr_ind(noanc_unphased_localanc_gs,g.true_anc_gs,ch,k),
	  " phased MOSAIC r^2 = ", dip_fr2_chr_ind(localanc_gs,g.true_anc_gs,ch,k), sep="", "\n")
    if (!PNG) readline()
    if(PNG) dev.off()
  }
}
