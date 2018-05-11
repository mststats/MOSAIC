# function to create a logfile for a MOSAIC model run
create_logfile=function(resultsdir,target,kLL,L,NUMI,firstind,chrnos,nchrno,NN,GpcM) {
  lognames<-c("mode", "time")
  for (j in 1:kLL) for (i in 1:L) lognames<-c(lognames, paste0("Mu[",j,",",i,"]"))
  for (i in 1:L) lognames<-c(lognames, paste0("rho[",i,"]"))
  for (ind in 1:NUMI) for (i in 1:L) for (j in 1:L) lognames<-c(lognames, paste0("PI[",ind,",",i,",",j,"]"))
  for (ind in 1:NUMI) for (i in 1:L) lognames<-c(lognames, paste0("alpha[",ind,",",i,"]"))
  for (ind in 1:NUMI) lognames<-c(lognames, paste0("lambda[",ind,"]"))
  for (i in 1:L) lognames<-c(lognames, paste0("theta[",i,"]"))
  lognames<-c(lognames, "loglikelihood")
  EMlogfile=paste0(resultsdir, target, "_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",format(Sys.time(), "%Y_%m_%d_%H:%M:%S"),"_EMlog.out")
  len=length(lognames) # total number of items to log
  write(file=EMlogfile, lognames, ncol=length(lognames)) # start the EM log file
  rtime<-as.numeric(Sys.time())
  return(list(runtime=rtime,logfile=EMlogfile,len=len))
}
