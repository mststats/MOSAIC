require(doParallel)
registerDoParallel(cores=8)
shargs<-commandArgs(trailingOnly=TRUE)
filename=shargs[1]
samelambda=as.logical(shargs[2])
recalc=as.logical(shargs[3])
pathin="RESULTS/"
load(paste0(pathin,filename))
source("coancestry.R")
if (recalc)
{
  load(paste0(pathin,"localanc_",filename))
  load(paste0(pathin,"noanc_unphased_localanc_",filename))
  new_acoancs=create_coancs(localanc,dr)
  new_coancs=create_coancs(noanc_unphased_localanc,dr)
}
cat("acoancs =     ", mean(plot_coanccurves(acoancs,dr,samelambda=samelambda,PLOT=F)$params[,,3]), "\n")
if (recalc)
  cat("new acoancs = ", mean(plot_coanccurves(new_acoancs,dr,samelambda=samelambda,PLOT=F)$params[,,3]), "\n")
cat("coancs =      ", mean(plot_coanccurves(coancs,dr,samelambda=samelambda,PLOT=F)$params[,,3]), "\n")
if (recalc)
  cat("new coancs =  ", mean(plot_coanccurves(new_coancs,dr,samelambda=samelambda,PLOT=F)$params[,,3]), "\n")

