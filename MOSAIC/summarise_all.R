#L=2
# script to summarise all model outputs and all panels as frequency/counts data along chromosomes
source("fst.R")
source("plot_funcs.R")
pathin="RESULTS/"
#pathout="/data/smew2/salter/HGDP_PLOTS/"
pathout="FREQS/"
datasource="HGDP/"
GpcM=60;chrnos=1:22
filenames=dir(pathin,glob2rx(paste0("*_", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData"))) # all targets
if (length(grep("localanc_", filenames))>0) 
  filenames=filenames[-grep("localanc_", filenames)] # remove localanc filenames
if (length(grep("gfbs_", filenames))>0) 
  filenames=filenames[-grep("gfbs_", filenames)] # remove gfbs filenames

#tmp=c("BantuSouthAfrica", "Brahui","Bulgarian","Cambodian","Daur","Druze","Georgian","Greek","HanNchina","Hazara","Hezhen","Hungarian","Indian","Makrani","Mandenka",
#      "Melanesian","Mozabite","NorthItalian","Oroqen","Pima","Polish","Romanian","SanNamibia","Tu","Turkish","Tuscan","Uygur","Uzbekistani","WestSicilian","Yemeni")

#tmp2=NULL
#for (i in 1:length(tmp))
#  tmp2=c(tmp2,grep(paste0(tmp[i],"_", L, "way"),filenames))

#tmp=scan("tmp_done.out",what="character",quiet=T)
#for (i in 1:length(tmp))
#  tmp2=c(tmp2,grep(tmp[i],filenames))

#tmp=c("Sardinian","SanKhomani", "Han", "Spanish", "Druze", "Bedouin", "Palestinian")
#tmp2=NULL
#for (i in 1:length(tmp))
#  tmp2=c(tmp2,grep(paste0(tmp[i],"_3way"),filenames))


#filenames=c(filenames[tmp2])#,filenames[-tmp2])
#filenames=c(filenames[-tmp2])

for (filename in filenames)
{
  cat("Looking at ", filename, "\n")
  if (exists("colvec")) rm(colvec)
  if (exists("localanc")) rm(localanc) # remove as otherwise can re-use from different target if missing. Prefer error in that case
  if (exists("finalflips")) rm(finalflips) # remove as otherwise can re-use from different target if missing. Prefer error in that case
  load(paste0(pathin,filename))
  load(paste0(pathin,"localanc_",filename))
  NUMI=length(alpha)
  flocalanc=phase_localanc(localanc,final.flips)
  ancestral_freqs=maximal_alleles(target,chrnos,flocalanc,datasource,datasource) 
  save(ancestral_freqs, file=paste0(pathout, target, "_", L, "way_", sum(NL), "_freqs.rdata"))
}
