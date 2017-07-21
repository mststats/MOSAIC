samelambda=F
asym=F
min.cM=0.5
targetdetails=gsub(".RData","",filename)
source("plot_funcs.R") # needs to be here as some data files include old plot code
source("coancestry.R")# needs to be here as some data files include old plot code
#source("fst.R");
# order via continents, etc
source("regions.R");ord.pops=rev(unlist(regional.pops));
#for (MODE in c("copy", "joint", "scaled"))
#{
  png(filename=paste0(pathout,targetdetails,"_",MODE,"_Mu.png"), width=1200, height=1920)
  tmp=match(ord.pops,rownames(Mu));tmp=tmp[!is.na(tmp)];Mu=Mu[tmp,]
  ord.Mu<-plot_Mu(Mu,alpha,MODE=MODE,cexa=2,shiftl=11,cutoff=0,ord=F)
  dev.off()
  png(filename=paste0(pathout,targetdetails,"_",MODE,"_Mu_beside.png"), width=800*L, height=1920)
  tmp=match(ord.pops,rownames(Mu));tmp=tmp[!is.na(tmp)];Mu=Mu[tmp,]
  ord.Mu<-plot_Mu(Mu,alpha,MODE=MODE,cexa=2,beside=T,shiftl=11,cutoff=0,ord=F)
  dev.off()
#}
#if ((L==2) | (L==3)) 
#{
  #this_Fst=Fst_combos(target,L,sum(NL),rownames(Mu))
  this_Fst=all_Fst[[which(names(all_Fst)==paste0(target,"_",L,"way_",sum(NL)))]]
  tmp=match(ord.pops,rownames(this_Fst$panels));tmp=tmp[!is.na(tmp)];this_Fst$panels=this_Fst$panels[tmp,]
  tmp_Fst=this_Fst # don't change this_Fst as used below for Rst calcs
  tmp_Fst$panels[is.nan(tmp_Fst$panels)]=1 # farthest possible if unknown
  png(filename=paste0(pathout,targetdetails,"_Fst.png"), width=1200, height=1920)
  ord.Fst<-plot_Fst(tmp_Fst$panels,cexa=3,ord=T, shiftl=14, cutoff=10)
  dev.off()
#}
if (dolocal | dostats)
{
  load(paste0(pathin, "localanc_", filename))
  localanc[[22]]=NULL # FLAG due to Myanmar bug
}
if (dolocal)
{
  y.lab="expected";cexa=2
  for (ch in 1:4) for (ind in 1:NUMI)
  {
    png(filename=paste0(pathout,"Localanc/",targetdetails,"_DIP_localanc.png"), width=1960, height=920)
    dipplot(ch,ind,g.loc[[ch]]*1e-6,localanc[[ch]],y.lab,cexa=cexa)
    mp<-axTicks(1,axp=c(min(g.loc[[ch]])*1e-6,max(g.loc[[ch]])*1e-6,5))
    axis(1,at=mp,labels=signif(mp,3))
    dev.off()
  }
}
# co-ancestry plots
if (L==2) {d1=1;d2=3} # dimensions of plots
if (L==3) {d1=2;d2=3} # dimensions of plots
if (L==4) {d1=2;d2=5} # dimensions of plots
if (L==5) {d1=3;d2=5} # dimensions of plots
if (L==6) {d1=3;d2=7} # dimensions of plots
acoancs=create_coancs(localanc,dr,"DIP",max.cM=50) # FLAG
png(filename=paste0(pathout,targetdetails,"_acoanc.png"), width=320*d2,height=320*d1)
this_acoplots=plot_coanccurves(acoancs,dr,lwd=4,cexa=2,verbose=F,axisall=F,samelambda=samelambda,asym=asym,min.cM=min.cM)
dev.off()
png(filename=paste0(pathout,targetdetails,"_coanc.png"), width=320*d2,height=320*d1)
this_coplots=plot_coanccurves(coancs,dr,lwd=4,cexa=2,verbose=F,axisall=F,samelambda=samelambda,asym=asym,min.cM=min.cM)
dev.off()
# EM log plot
doEM=T
if (doEM)
{
  EMlogfile=dir(pathin,glob2rx(paste0(strsplit(filename,sum(NL))[[1]][1],sum(NL),"*merged_EMlog.out")))
  source("log_funcs.R")
  EMlog<-extract_log(paste0(pathin,EMlogfile))
  png(filename=paste0(pathout,targetdetails,"_EMlog.png"), width=960,height=960)
  plot_loglike(EMlog)
  dev.off()
}
if (dostats)
{
  # report expected r^2 values
  source("calc_r2.R")
  Rst=Q_Fst(this_Fst)
  tmp=c(targetdetails,min(Reduce("+",alpha)/NUMI),mean(this_acoplots$params[,,3]),mean(this_coplots$params[,,3]),signif(dip_expected_fr2(localanc),3),this_Fst$anc,Rst)
  write(tmp, file=paste0(pathout, "stats",L,"way.txt"), append=T,ncol=length(tmp))
}
