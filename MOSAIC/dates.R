L=2
samelambda=T
asym=F
min.cM=0.5
thresh=1e-4 # minimum amount of minor ancestry to consider using an individual for the calculation
source("coancestry.R")
require(doParallel);MC=as.integer(detectCores()/2);registerDoParallel(cores=MC)
source("fst.R")
pathin="RESULTS/"
pathout="PLOTS/"
datasource="HGDP/"
GpcM=60;chrnos=1:22
filenames=dir(pathin,glob2rx(paste0("*_", L,"way*", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData"))) # all targets
tmp=grep("localanc_", filenames);if (length(tmp)>0) filenames=filenames[-tmp] # remove localanc filenames
tmp=grep("gfbs_", filenames);if (length(tmp)>0) filenames=filenames[-tmp] # remove gfbs filenames
remfiles=c("Spanish","NorthAfrican") # special cases to remove (if there)

tmp=c("BantuSouthAfrica","Brahui","Bulgarian","Cambodian","Daur","Druze","Georgian","Greek","HanNchina","Hazara","Hezhen","Hungarian","Indian","Makrani","Mandenka",
      "Melanesian","Mozabite","NorthItalian","Oroqen","Pima","Polish","Romanian","SanNamibia","Tu","Turkish","Tuscan","Uygur","Uzbekistani","WestSicilian","Yemeni")
tmp2=NULL
for (i in 1:length(tmp))
  tmp2=c(tmp2,grep(glob2rx(paste0(tmp[i],"_2way*", "1-22_2980*")),filenames))
filenames=c(filenames[tmp2])#,filenames[-tmp2])

stats=read.table(paste0(pathout,"stats2way.txt"),stringsAsFactors=F,header=T)
stats=stats[grep("2980_60",stats[,1]),]
stats[,1]=sapply(strsplit(stats[,1],"_2way"),function(x) x[[1]])

load("all_Fst_2.rdata")
for (remfile in remfiles)
{
  tmp=grep(remfile,filenames)
  if (length(tmp)!=0)
    filenames=filenames[-tmp]
}
targets=NULL
source("plot_funcs.R")
tmp=c("target", "min alpha", "Fst", "Rst", "r2", "anc mean","anc SE","anc boot mean", "anc boot sd", "noanc mean","noanc SE","noanc boot mean", "noanc boot sd")
write(tmp,file="mosaic_dates.csv",sep=",",ncol=length(tmp))
for (i in 1:length(filenames))
{
  filename=filenames[i]
  cat("Looking at ", filename, "\n")
  load(paste0(pathin,filename))
  target=strsplit(filename,"_")[[1]][1]
  L=as.integer(strsplit(strsplit(filename,"_")[[1]][2],"way")[[1]][1])
  firstind=as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]][1])
  NUMI=diff(as.integer(strsplit(strsplit(filename,"_")[[1]][3],"-")[[1]]))+1
  NUMA=NUMI*2
  tmp=as.integer(strsplit(strsplit(filename,"_")[[1]][4],"-")[[1]]);chrnos=tmp[1]:tmp[2]
  NN=as.integer(strsplit(filename,"_")[[1]][5])
  GpcM=as.integer(strsplit(filename,"_")[[1]][6])
  prop.don=as.numeric(strsplit(filename,"_")[[1]][7])
  max.donors=as.integer(strsplit(strsplit(filename,"_")[[1]][8],"R")[[1]][1])
  ##########################
  targets[i]=target
  load(paste0(pathin,"localanc_",filename))
  load(paste0(pathin,"noanc_unphased_localanc_",filename))
  new_acoancs=create_coancs(localanc,dr)
  new_coancs=create_coancs(noanc_unphased_localanc,dr)
  tmp_fst=all_Fst[[which(names(all_Fst)==paste0(target,"_",L,"way_",NN))]]
  tmpdates=c(target, min(Reduce("+",alpha)/NUMI),tmp_fst$anc,stats$Rst[match(target,stats[,1])],stats$r2[match(target,stats[,1])])
  kgens=rep(NaN,NUMI);for (k in 1:NUMI) if (min(alpha[[k]])>thresh) kgens[k]=mean(plot_coanccurves(new_coancs,dr,k=k,PLOT=F,samelambda=samelambda,asym=asym,min.cM=min.cM)$params[,,3],na.rm=T)
  akgens=rep(NaN,NUMI);for (k in 1:NUMI) if (min(alpha[[k]])>thresh) akgens[k]=mean(plot_coanccurves(new_acoancs,dr,k=k,PLOT=F,samelambda=samelambda,asym=asym,min.cM=min.cM)$params[,,3],na.rm=T)
  nsamps=100
  G=sapply(localanc,function(x) dim(x)[3])
  aboot.gens=boot.gens=rep(NaN,nsamps)
  aboot.localanc=boot.localanc=list()
  schrno=21 # FLAG use 1:schrno; set to 3=21 to avoid bug in data
  for (t.ch in 1:schrno) 
  {
    aboot.localanc[[t.ch]]=boot.localanc[[t.ch]]=array(NaN,c(L,NUMA,G[t.ch]))
  }
  pb<-txtProgressBar(min=0,max=nsamps,style=3)
  for (r in 1:nsamps) ## nsamps bootstrap samples
  {
    setTxtProgressBar(pb, r)
    #bootstrap chromosomes; generate NUMI pseudo-individuals using random chromosomes drawn from all inds w/ replacement; see GT S4.4 for details
    boot.haps=matrix(sample(1:NUMA,NUMA*schrno,replace=T),NUMA)
    boot.inds=matrix(sample(1:NUMI,NUMI*schrno,replace=T),NUMI)
    aboot.haps=matrix(NaN,NUMA,schrno);for (t.ch in 1:schrno) for (ind in 1:NUMI) aboot.haps[c(ind*2-1,ind*2),t.ch]=c(boot.inds[ind,t.ch]*2-1,boot.inds[ind,t.ch]*2)
    boot.haps=aboot.haps
    for (t.ch in 1:schrno) 
    {
      for (l in 1:L) for (hap in 1:NUMA)
      {
        aboot.localanc[[t.ch]][l,hap,]=localanc[[t.ch]][l,aboot.haps[hap,t.ch],]
        boot.localanc[[t.ch]][l,hap,]=noanc_unphased_localanc[[t.ch]][l,boot.haps[hap,t.ch],]
      }
    }
    aboot.coancs=create_coancs(aboot.localanc,dr,"DIP",max.cM=50)#*mean(unlist(lambda))/100);
    boot.coancs=create_coancs(boot.localanc,dr,"DIP",max.cM=50)#*mean(unlist(lambda))/100);
    aboot.gens[r]=mean(plot_coanccurves(aboot.coancs,dr,PLOT=F,samelambda=samelambda,asym=asym,min.cM=min.cM)$params[,,3],na.rm=T)
    boot.gens[r]=mean(plot_coanccurves(boot.coancs,dr,PLOT=F,samelambda=samelambda,asym=asym,min.cM=min.cM)$params[,,3],na.rm=T)
  }  
  close(pb)
  tmpdates=c(tmpdates,mean(akgens,na.rm=T),sd(akgens,na.rm=T)/sqrt(sum(!is.na(akgens))),mean(aboot.gens,na.rm=T),sd(aboot.gens,na.rm=T))
  tmpdates=c(tmpdates,mean(kgens,na.rm=T),sd(kgens,na.rm=T)/sqrt(sum(!is.na(kgens))),mean(boot.gens,na.rm=T),sd(boot.gens,na.rm=T))
  write(tmpdates,file="mosaic_dates.csv",sep=",",append=T,ncol=length(tmpdates))
}
gens_to_date=function(x) 1950-28*(x+1) # note that GlobeTrotter www states this wrongly; see page 28 of supplement
date_to_gens=function(x) (1950-x)/28-1 # note that GlobeTrotter www states this wrongly; see page 28 of supplement



