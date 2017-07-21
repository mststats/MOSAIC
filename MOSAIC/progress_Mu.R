PDF=T
PNG=F
cexa=1
#filename="simulated_3way_1-4_15-22_770_60_0.99_100.RData"
pathin="RESULTS/"
pathout="Muplots/"
load(paste0(pathin,filename))
source("plot_funcs.R")
source("log_funcs.R")
panels=rownames(Mu)
NUMI=NUMA/2
ord=apply(Mu/NL[1:kLL],2,order) # ordering of panels in each ancestry by final copying matrix
EMlogfile=dir(pathin,glob2rx(paste0(strsplit(filename,sum(NL))[[1]][1],sum(NL),"*_EMlog.out")))
EMlog=extract_log(paste0(pathin,EMlogfile))
iter=6 # skip Q only updates
iterEM=1
for (iter in iter:nrow(EMlog))
{
  if (as.character(EMlog[iter,1])=="EM")
  {
    paras=extract_paras(EMlog,iter,rownames(Mu))
    all.alpha=Reduce("+",paras$alpha);all.alpha=all.alpha/sum(all.alpha);all.alpha=round(all.alpha,3)
    # note we are looking at Mu/NL*alpha
    tmpMu=paras$Mu/NL[1:kLL];tmpMu=t(t(tmpMu)/colSums(tmpMu))
    ordMu=NULL;for (a in 1:L) ordMu[[a]]=tmpMu[ord[,a],a]#*all.alpha[a] # order this iteration by final Mu
    xmax=max(unlist(ordMu))
    if (PDF) pdf(file=paste0(pathout,paste0(target,collapse="-"),"_iter_",iterEM,".pdf"),width=10)
    if (PNG) png(file=paste0(pathout,paste0(target,collapse="-"),"_iter_",iterEM,".png"))
    par(mfrow=c(1,L),mar=c(2,10,2,1))
    for (a in 1:L) 
    {
      mlab=ifelse(a==1,iterEM, "")
      y=barplot(ordMu[[a]],horiz=T,las=1,col=colvec[a],xlim=c(0,xmax),cex.names=cexa,cex.axis=cexa,main=mlab,cex.main=cexa*2)
      text(xmax/4,max(y)+cexa*3,all.alpha[a],col=colvec[a],cex=cexa)
    }
    if (PDF | PNG) dev.off()
    if (!PDF & !PNG) readline()
    iterEM=iterEM+1
  }
}
