path<-"/data/sanderling/salter/Data/SpanishDataFromClare/data/"
tmp<-scan(paste0(path,"SP_PR_HM_NA.chr",chrno,"_v1.chromo.haps"),what="character")
allNN<-as.integer(tmp[2])*2 # #inds X 2 for #haps
S<-as.integer(tmp[3])
locs<-as.integer(tmp[5:(4+S)])
# now read in population information
allpops<-read.table(paste0(path,"LOO-SP_PR_HM_NA.gtv14.input.file.ids"))
# now reduce to parts we need
rem<-rep(F,nrow(allpops)); # remove none to start
rem[agrep("exclude",allpops[,2])]=T
rem[allpops[,3]==0]=T # remove individuals with 0 in third column
keep<-(!rem)
pops=unique(allpops[,2])
NN<-sum(keep)*2 # #haps to keep
keep<-which(keep)
allpops=allpops[keep,]
write.table(allpops[,c(2,1,3)],file="spanish/sample.names",row.names=F,col.names=F,quote=F)
Y<-matrix(0L,NN,S)
f<-function(n) as.integer(strsplit(tmp[5+n+S],fixed=T,split="")[[1]])
 # this has been checked by comparison with finding all Y and subsetting
for (i in 1:(NN/2))
{
  for (h in 1:2) # get both haps
  {
    n=(keep[i]-1)*2+h # find index
    Y[(i-1)*2+h,]<-f(n)
  }
}
rm(tmp)
rownames(Y)<-rep(allpops[,2],each=2)
pops<-unique(rownames(Y))

# create list to store populations
for (i in 1:length(pops))
{
 y=Y[rownames(Y)==pops[i],]
 write.table(t(y),file=paste0("spanish/",pops[i],"genofile.",chrno),sep="",col.names=F,row.names=F)
}
# create matrix of snps for which we have haplotypes
snps<-matrix(NaN, S, 6) # same size as hapmix files, most will be left blank here
snps[,2]<-chrno
snps[,4]<-locs
write.table(snps, file=paste0("spanish/snpfile.", chrno)) 


rates<-read.table(paste0(path,"genetic_map_chr", chrno, "_combined_b37.txt"),header=T)
write(paste0(":sites:",nrow(rates)), paste0("spanish/rates.",chrno))
write(rates[,1],paste0("spanish/rates.",chrno),append=T,ncol=nrow(rates))
write(rates[,3],paste0("spanish/rates.",chrno),append=T,ncol=nrow(rates)) # in centimorgans 
