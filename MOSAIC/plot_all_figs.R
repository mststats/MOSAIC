#rm(list=ls());chrnos=22:1;target="Bedouin";L=2;
require(doParallel);MC=4;registerDoParallel(cores=MC)
MODE="scaled"
dolocal=F
dostats=T
source("fst.R")
source("plot_funcs.R")
#pathin="/data/smew2/salter/HGDP_RESULTS/"
pathin="RESULTS/"
#pathin="/media/mst/9C66B0D766B0B37E/sanderling/smew_MOSAIC_RESULTS/RESULTS/"
#pathout="/data/smew2/salter/HGDP_PLOTS/"
pathout="PLOTS/"
GpcM=60;chrnos=1:22
filenames=dir(pathin,glob2rx(paste0("*_", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData"))) # all targets
if (length(grep("localanc_", filenames))>0) 
  filenames=filenames[-grep("localanc_", filenames)] # remove localanc filenames
if (length(grep("gfbs_", filenames))>0) 
  filenames=filenames[-grep("gfbs_", filenames)] # remove gfbs filenames
#target="Kalash"
#filenames=dir(pathin,glob2rx(paste0(target, "*_", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData"))) # particular targets

# only take smaller of each as this is the post group dropping result
#tmp<-NULL
#for (firstind in unique(as.integer(sapply(strsplit(sapply(strsplit(filenames,"way_"),function(x) x[[2]]),"_"),function(x) x[1]))))
#  tmp<-c(tmp,sort(filenames[grep(paste0(target,"_",L,"way_",firstind),filenames)])[2])
#filenames<-filenames[-match(tmp,filenames)]
load("all_Fst_2.rdata");all_Fst2=all_Fst;load("all_Fst_3.rdata");all_Fst3=all_Fst; load("all_Fst_4.rdata");all_Fst4=all_Fst;
if (dostats)
{
  write(paste("TargetDetails", "min.alpha", "mean.a.lambda", "mean.lambda", "r2","Fst", "Rst"), file=paste0(pathout, "stats2way.txt"))
  write(paste("TargetDetails", "min.alpha", "mean.a.lambda", "mean.lambda", "r2","Fst_1x2", "Fst_1x3", "Fst_2x3", "Rst_1x2", "Rst_1x3", "Rst_2x3"), file=paste0(pathout, "stats3way.txt"))
  write(paste("TargetDetails", "min.alpha", "mean.a.lambda", "mean.lambda", "r2","Fst_1x2", "Fst_1x3", "Fst_1x4", "Fst_2x3", "Fst_2x4", "Fst_3x4",
                                                                "Rst_1x2", "Rst_1x3", "Rst_1x4", "Rst_2x3", "Rst_2x4", "Rst_3x4"), file=paste0(pathout, "stats4way.txt"))
}
for (filename in filenames)
{
  cat("Looking at ", filename, "\n")
  if (exists("colvec")) rm(colvec)
  load(paste0(pathin,filename))
  if (L==2) all_Fst=all_Fst2
  if (L==3) all_Fst=all_Fst3
  if (L==4) all_Fst=all_Fst4
  NUMI=length(alpha)
  #neworder<-sort(Reduce("+",alpha)/NUMI,index=T)$ix;source("reorder_ancs.R") # order from low to high
  source("plot_figs.R")
}

