L=3
#num=3 # how many panels to use to summarise each ancestry
source("fst.R")
pathin="RESULTS/"
datasource="HGDP/"
GpcM=60;chrnos=1:22
filenames=dir(pathin,glob2rx(paste0("*_", L,"way*", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData"))) # all targets
if (length(grep("localanc_", filenames))>0) 
  filenames=filenames[-grep("localanc_", filenames)] # remove localanc filenames
if (length(grep("gfbs_", filenames))>0) 
  filenames=filenames[-grep("gfbs_", filenames)] # remove gfbs filenames

#remfiles=c("NorthAfrican") # special cases to remove (if there)
#tmp2=NULL
#for (i in 1:length(remfiles))
#  tmp2=c(tmp2,grep(remfiles[i],filenames))
#if (length(tmp2)>0)
#  filenames=filenames[-tmp2]

#tmp=c("BantuSouthAfrica", "Brahui","Bulgarian","Cambodian","Daur","Druze","Georgian","Greek","HanNchina","Hazara","Hezhen","Hungarian","Indian","Makrani","Mandenka",
#      "Melanesian","Mozabite","NorthItalian","Oroqen","Pima","Polish","Romanian","SanNamibia","Tu","Turkish","Tuscan","Uygur","Uzbekistani","WestSicilian","Yemeni")
#tmp2=NULL
#for (i in 1:length(tmp))
#  tmp2=c(tmp2,grep(paste0(tmp[i],"_2way"),filenames))
#filenames=c(filenames[tmp2])#,filenames[-tmp2])

all_Fst=list()
source("plot_funcs.R")
for (i in 1:length(filenames))
{
  filename=filenames[i]
  cat("Looking at ", filename, "\n")
  load(paste0(pathin,filename))
  target=strsplit(filename,"_")[[1]][1]
  L=as.integer(strsplit(strsplit(filename,"_")[[1]][2],"way")[[1]][1])
  firstind=as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]][1])
  NUMI=diff(as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]]))+1
  tmp=as.integer(strsplit(strsplit(filename,"_")[[1]][4],"-")[[1]]);chrnos=tmp[1]:tmp[2]
  NN=as.integer(strsplit(filename,"_")[[1]][5])
  GpcM=as.integer(strsplit(filename,"_")[[1]][6])
  prop.don=as.numeric(strsplit(filename,"_")[[1]][7])
  max.donors=as.integer(strsplit(strsplit(filename,"_")[[1]][8],"R")[[1]][1])
  ##########################
  all_Fst[[i]]=Fst_combos(target,L,NN,rownames(Mu))
  names(all_Fst)[i]=paste0(target,"_",L,"way_",NN)
}
save(all_Fst, file=paste0("all_Fst_", L, ".rdata"))
