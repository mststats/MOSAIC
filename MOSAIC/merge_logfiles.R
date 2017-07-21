source("log_funcs.R")
pathin="RESULTS/"
GpcM=60
filenames=dir(pattern=glob2rx("*out"),path=pathin)
# minus already merged ones
tmp1=grep("merged",filenames)
if (length(tmp1)>0)
  filenames=filenames[-tmp1]

for (filename in filenames)
{
  target=strsplit(filename,"_")[[1]][1]
  L=as.integer(strsplit(strsplit(filename,"_")[[1]][2],"way")[[1]][1])
  firstind=as.integer(strsplit(strsplit(filename,"way_")[[1]][2],"-")[[1]][1])
  NUMI=as.integer(strsplit(strsplit(filename,"-")[[1]][2],"_")[[1]][1])-firstind+1
  chrnos=as.integer(c(strsplit(strsplit(filename,"-")[[1]][2],"_")[[1]][2],strsplit(strsplit(filename,"-")[[1]][3],"_")[[1]][1]))
  nchrno=diff(chrnos)+1
  NN=as.integer(tail(strsplit(strsplit(filename,paste0("_",GpcM))[[1]][1],"_")[[1]],1))
  tmp=paste0(target,"_",L,"way_",firstind,"-",NUMI,"_",paste(chrnos,collapse="-"),"_",NN)
  #filenames=filenames[-grep(tmp,filenames)] # remove the matching ones from the list
  cat("Looking at ", tmp, "\n")
  mfilenames=sort(dir(pattern=glob2rx(paste0(tmp,"*out")),path=pathin))
  tmp1=grep("merged",mfilenames)
  if (length(tmp1)>0)
    mfilenames=mfilenames[-tmp1]
  EMlog=NULL
  for (mfilename in mfilenames)
  {
    tmp=extract_log(paste0(pathin,mfilename))
    # now remove overlaps
    if (!is.null(EMlog))
      tmp=tmp[tmp[,ncol(tmp)]>EMlog[nrow(EMlog),ncol(EMlog)],]
    tmp$time[is.na(tmp$time)]=0
    EMlog=rbind(EMlog,tmp)
  }
  EMlogfile=paste0(pathin, target, "_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos,collapse="-"),"_",NN,"_",GpcM,"_merged_EMlog.out")
  write.table(EMlog,file=EMlogfile,row.names=F,quote=F)
}
