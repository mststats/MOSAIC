# extract parameters from a specified log file
extract_log=function(logfile)
{
  EMlog=read.delim(logfile, sep=" ")
  colnames(EMlog)<-gsub("[.]","",colnames(EMlog))
  return(EMlog)
}
extract_paras=function(EMlog, iter, panelnames=NULL)
{
  paras=as.numeric(c(EMlog[iter,-c(1,2,ncol(EMlog))])) # remove mode, time, and log-likelihood
  t.Mu=t(matrix(paras[1:(L*kLL)],L));paras=paras[-(1:(L*kLL))] 
  rownames(t.Mu)=panelnames
  t.rho=paras[1:L];paras=paras[-(1:L)]
  t.Q=list();for (ind in 1:NUMI) {t.Q[[ind]]=matrix(paras[1:(L*L)],L);paras=paras[-(1:(L*L))]}
  t.alpha=list();for (ind in 1:NUMI) {t.alpha[[ind]]=paras[1:L];paras=paras[-(1:L)]}
  t.lambda=list();for (ind in 1:NUMI) {t.lambda[[ind]]=paras[1];paras=paras[-1]}
  t.theta=paras[1:L];paras=paras[-(1:L)]
  return(list(Mu=t.Mu, rho=t.rho, Q=t.Q, alpha=t.alpha, lambda=t.lambda, theta=t.theta))
}
# EMlog<-extract_log(logfile)
# e.g. paras=extract_paras(EMlog, nrow(EMlog))
plot_loglike=function(EMlog,cexa=2.5,
		      colvec=c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#0072B2", "#999999"),...)
{

  par(cex.lab=cexa, mar=c(5,5,0,0))
  plot(cumsum(EMlog$time),EMlog$loglikelihood,col=8,t='l',ylab="log-likelihood", xlab="time in seconds",lwd=cexa,...)
  points(cumsum(EMlog$time),EMlog$loglikelihood,col=colvec[EMlog$mode],t='p',pch=20,cex=cexa)
  legend("bottomright",levels(EMlog$mode),col=colvec[1:length(levels(EMlog$mode))],pch=20,cex=cexa)
}

