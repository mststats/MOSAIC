# code to load multiple individuals from a target population
path="RESULTS/"
filenames=dir(path,glob2rx(paste0("reduced_", target, "_", L, "way_*", paste(chrnos[c(1,length(chrnos))],collapse="-"),"_*_",GpcM,"*.RData")))
NUMI=length(filenames)
load(paste0(path,filenames[1])) # load the first one just for dimensions
G<-sapply(localanc,function(x) dim(x)[[3]])
allinds.localanc<-list()
allinds.Mu<-NULL
allinds.NL<-NULL
allinds.alpha<-NULL
for (ch in 1:nchrno)
  allinds.localanc[[ch]]<-array(NaN, c(L,NUMI*2,G[ch]))
for (ind in 1:NUMI) 
{
  cat("loading ind", ind, "\n")
  hap<-c(ind*2-1,ind*2)
  load(paste0(path,filenames[ind]))
  for (ch in 1:nchrno) 
    allinds.localanc[[ch]][,hap,]<-localanc[[ch]]
  allinds.Mu<-rbind(allinds.Mu,Mu)
  allinds.NL<-c(allinds.NL,NL[1:kLL])
  allinds.alpha<-rbind(allinds.alpha,alpha[[1]])
}
NUMA=NUMI*2
# permute to best label matching
#use reorder_ancs.R

# now average over panels for Mu and drop any that don't always appear
summed.Mu<-matrix(NaN, length(unique(rownames(allinds.Mu))), L)
for (ll in 1:length(unique(rownames(allinds.Mu))))
  for (i in 1:L)
    summed.Mu[ll,i]<-sum(allinds.Mu[rownames(allinds.Mu)==unique(rownames(allinds.Mu))[ll],i]) # sum rather than mean as some appear more often than others

summed.Mu<-t(t(summed.Mu)/colSums(summed.Mu))
averaged.NL<-c(tapply(allinds.NL, rownames(allinds.Mu), mean))
rownames(summed.Mu)<-unique(rownames(allinds.Mu))
averaged.alpha<-apply(allinds.alpha,2,mean)

#source("plot_funcs.R");MODE="copy";cexa=1.3;for (ind in 1:30) {pdf(file=paste0("PLOTS/Mu_",target,"_",L,"way_",ind,".pdf"),width=14,height=21);ind.Mu=allinds.Mu[((ind-1)*kLL+1):(ind*kLL),];ind.alpha=allinds.alpha[ind,];ord.Mu<-plot_Mu(ind.Mu,list(ind.alpha),MODE=MODE,cexa=cexa);dev.off();}

