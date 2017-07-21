# recalculate and save localanc that is ancestry unaware
pathin=resultsdir="RESULTS/"
GpcM=60;chrnos=1:22
filenames=dir(pathin,glob2rx(paste0("*_", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData"))) # all targets
if (length(grep("localanc_", filenames))>0) 
  filenames=filenames[-grep("localanc_", filenames)] # remove localanc filenames
if (length(grep("gfbs_", filenames))>0) 
  filenames=filenames[-grep("gfbs_", filenames)] # remove gfbs filenames
require(doParallel)
registerDoParallel(cores=16)

#filenames="simulated_2way_1-1_22-22_396_60_0.99_100.RData";chrnos=22 # FLAG
tmp=c("BantuSouthAfrica","Brahui","Bulgarian","Cambodian","Daur","Druze","Georgian","Greek","HanNchina","Hazara","Hezhen","Hungarian","Indian","Makrani","Mandenka",
      "Melanesian","Mozabite","NorthItalian","Oroqen","Pima","Polish","Romanian","SanNamibia","Tu","Turkish","Tuscan","Uygur","Uzbekistani","WestSicilian","Yemeni")
tmp2=NULL
for (i in 1:length(tmp))
  tmp2=c(tmp2,grep(glob2rx(paste0(tmp[i],"_2way*", "1-22_2980*")),filenames))
#filenames=c(filenames[tmp2])
filenames=c(filenames[-tmp2])

for (filename in filenames[8:length(filenames)])
{
  mask=NULL
  cat("Looking at ", filename, "\n")
  target=strsplit(filename,"_")[[1]][1]
  L=as.integer(strsplit(strsplit(filename,"_")[[1]][2],"way")[[1]][1])
  firstind=as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]][1])
  NUMI=diff(as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]]))+1
  tmp=as.integer(strsplit(strsplit(filename,"_")[[1]][4],"-")[[1]]);chrnos=tmp[1]:tmp[2]
  NN=as.integer(strsplit(filename,"_")[[1]][5])
  GpcM=as.integer(strsplit(filename,"_")[[1]][6])
  prop.don=as.numeric(strsplit(filename,"_")[[1]][7])
  max.donors=as.integer(strsplit(strsplit(filename,"_")[[1]][8],"R")[[1]][1])
  source("reload.R") # recompute umatch matrices
  get_switches=F;LOG=T;PLOT=F;eps.lower=log(2);min.bg=0.1;max.bg=1.0;require(boot);o.LOG=F
  for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # undo phase flips
  EM=F # turn off EM updates
  getnoancgfbs=T
  samp_chrnos=chrnos;subNUMA=NUMA;subNL=max(NL) # use them all
  source("noanc.R")
  source("cleanup.R")
  Mu=a.Mu;rho=a.rho;theta=a.theta;Q=a.Q;alpha=a.alpha;lambda=a.lambda
  source("ancunaware.R")
  noanc_unphased_localanc=get_ancunaware_localanc() # works off noanc_gfbs
  cat("saving results to file\n")
  # resave over original file
  save(file=paste0(pathin,"noanc_unphased_localanc_",target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),
   	     "_",NN,"_",GpcM,"_",prop.don,"_",max.donors,".RData"), noanc_unphased_localanc)
}

