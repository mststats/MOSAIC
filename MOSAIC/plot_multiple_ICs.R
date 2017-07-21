useIC="AIC"
maxL=5
path="RESULTS/"
#filenames=dir(path,glob2rx(paste0("init*.RData")))
#filenames=dir(path,glob2rx(paste0("init_Yoruba*.RData")))
filenames=dir(path,glob2rx(paste0("init_SpainPopn_*.RData")))
#filenames="init_simulated_2_2.RData"
for (filename in filenames)
{
  load(paste0(path,"/",filename))
  PNG=T # important that this comes after load()
  #if (PNG) png(file=paste0("PLOTS/",strsplit(filename,".RData")[[1]],".png"),width=1840,height=1840) # use for HGDp
  if (PNG) png(file=paste0("PLOTS/",strsplit(filename,".RData")[[1]],".png"),width=1840,height=920) # use for Spain
  source("plot_ICs.R")
  if (PNG) dev.off()
}
