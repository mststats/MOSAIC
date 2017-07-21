#useIC="AIC"
#source("sourceCppcalls.R")
verbose=T
source("init_Mu.R")
#windowed_copying<-window_chunks(nswitches=noanc_gswitches,ww=ww,verbose=verbose)
maxic=-1e10
allIC<-NULL
all.Mu<-list()
all.alpha<-list()
for (L in 1:maxL) 
{
  tmp<-cluster_windows(windowed_copying,PLOT=F,t.L=L,method="EMmult",verbose=verbose)
  nparas=((kLL-1)*L+L-1) # Mu paras + alpha paras
  tmpBIC<-tmp$ll-log(sum(sapply(windowed_copying$wmat,ncol)))*nparas
  tmpAIC<-tmp$ll-2*nparas
  tmpIC=tmp$ll-2*L*(nparas) # a load of shite
  allIC<-rbind(allIC,c(L,tmp$ll,tmpBIC,tmpAIC,tmpIC))
  newIC=ifelse(useIC=="AIC", tmpAIC, tmpIC)
  all.Mu[[L]]<-tmp$Mu
  all.alpha[[L]]<-tmp$alpha
  if (newIC>maxic)
  {
    Mu=tmp$Mu
    alpha=tmp$alpha
    maxic=newIC
    bestL=L
  }
}
par(mfrow=c(1,2))
source("plot_funcs.R");
ord.Mu=plot_Mu(Mu,list(alpha),NL,cexa=2)
par(mar=c(4,0,2,4))
plot(allIC[,c(1,3)],t='b',xlab="#ancs", ylab="BIC",main=paste0(target,"_",bestL,"_",NUMA),lwd=2)
par(new=T,mar=c(4,0,2,4))
plot(allIC[,c(1,4)],t='b',xaxt='n',yaxt='n',xlab='',ylab='',col=2)
axis(4,col=2)
par(new=T,mar=c(4,0,2,4))
plot(allIC[,c(1,5)],t='b',xaxt='n',yaxt='n',xlab='',ylab='',col=3)
axis(4,col=3,line=2)
legend("bottomright", c("BIC","AIC","IC"),col=1:3,lty=c(1,1))
