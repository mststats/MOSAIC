# function to create a logfile for a MOSAIC model run
create_logfile=function(resultsdir,target,kLL,A,NUMI,firstind,chrnos,nchrno,NN,GpcM) {
  lognames<-c("mode", "time")
  for (j in 1:kLL) for (i in 1:A) lognames<-c(lognames, paste0("Mu[",j,",",i,"]"))
  for (i in 1:A) lognames<-c(lognames, paste0("rho[",i,"]"))
  for (ind in 1:NUMI) for (i in 1:A) for (j in 1:A) lognames<-c(lognames, paste0("PI[",ind,",",i,",",j,"]"))
  for (ind in 1:NUMI) for (i in 1:A) lognames<-c(lognames, paste0("alpha[",ind,",",i,"]"))
  for (ind in 1:NUMI) lognames<-c(lognames, paste0("lambda[",ind,"]"))
  for (i in 1:A) lognames<-c(lognames, paste0("theta[",i,"]"))
  lognames<-c(lognames, "loglikelihood")
  logfile=paste0(resultsdir, target, "_", A, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",format(Sys.time(), "%Y_%m_%d_%H:%M:%S"),"_EMlog.out")
  len=length(lognames) # total number of items to log
  write(file=logfile, lognames, ncol=length(lognames)) # start the EM log file
  rtime<-as.numeric(Sys.time())
  return(list(runtime=rtime,logfile=logfile,len=len))
}

# function to extract parameters from a specified log file
extract_log=function(logfile)
{
  EMlog=read.delim(logfile, sep=" ")
  EMlog$mode=as.factor(EMlog$mode) # required since R 4.0 as stringsAsFactors now defaults to FALSE
  colnames(EMlog)<-gsub("[.]","",colnames(EMlog))
  return(EMlog)
}
extract_paras=function(EMlog, iter, panelnames=NULL)
{
  paras=as.numeric(c(EMlog[iter,-c(1,2,ncol(EMlog))])) # remove mode, time, and log-likelihood
  t.Mu=t(matrix(paras[1:(A*kLL)],A));paras=paras[-(1:(A*kLL))] 
  rownames(t.Mu)=panelnames
  t.rho=paras[1:A];paras=paras[-(1:A)]
  t.PI=list();for (ind in 1:NUMI) {t.PI[[ind]]=matrix(paras[1:(A*A)],A);paras=paras[-(1:(A*A))]}
  t.alpha=list();for (ind in 1:NUMI) {t.alpha[[ind]]=paras[1:A];paras=paras[-(1:A)]}
  t.lambda=list();for (ind in 1:NUMI) {t.lambda[[ind]]=paras[1];paras=paras[-1]}
  t.theta=paras[1:A];paras=paras[-(1:A)]
  return(list(Mu=t.Mu, rho=t.rho, PI=t.PI, alpha=t.alpha, lambda=t.lambda, theta=t.theta))
}
# EMlog<-extract_log(logfile)
# e.g. paras=extract_paras(EMlog, nrow(EMlog))
plot_loglike=function(EMlog,cexa=2.5,
		      colvec=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#0072B2", "#999999"),...)

{

  par(cex.lab=cexa, mar=c(5,5,0,0))
  plot(cumsum(EMlog$time),EMlog$loglikelihood,col=8,t='l',ylab="log-likelihood", xlab="time in seconds",lwd=cexa,...)
  points(cumsum(EMlog$time),EMlog$loglikelihood,col=colvec[EMlog$mode],t='p',pch=20,cex=cexa)
  legend("bottomright",levels(EMlog$mode),col=colvec[1:length(levels(EMlog$mode))],pch=20,cex=cexa)
}

# function to write a row in the log file
  writelog<-function(t.logfile,t.alg,t.diff.time,t.len,t.Mu,t.rho,t.PI,t.alpha,t.lambda,t.theta,t.cloglike) # single consistent function to write to logfile
    write(file=t.logfile,c(t.alg,signif(t.diff.time,4),signif(t(t.Mu),4),signif(t.rho,4),c(sapply(t.PI, function(x) signif(t(x),4))),
			   sapply(t.alpha, function(x) signif(x,4)),sapply(t.lambda,function(x) round(x,4)),signif(t.theta,4),round(t.cloglike,4)),ncol=t.len,append=T)
