# script to convert shapeit2 .haps output to MOSAIC input type 1
# reads in rates file and .haps files (one per chromosome) and outputs rates file and haps files (one per population per chromosome)
pathin="~/Data/HGDP_shapeit/"
chrno=22
hapsfile=paste0("Hellenthal2013IIChrom",chrno,"phased.haps")
#ratesfile=""
inds.data="sample.names"
pathout="TMP/"

shapeithaps=read.table(paste0(pathin,hapsfile))
S=nrow(shapeithaps)
NN=ncol(shapeithaps)
locs=shapeithaps[,3]

# now read in population information
allpops=read.table(paste0(pathin,inds.data),header=TRUE)[-1,]
# now reduce to parts we need
pops=as.character(unique(allpops[,1]))
keep=rep(T,nrow(allpops)); # remove none to start
keep=which(keep)
allpops=allpops[keep,]
write.table(allpops[,c(1,2,3)],file=paste0(pathout,inds.data),row.names=F,col.names=F,quote=F)
hap.pops=rep(allpops[,1],each=2)

# create list to store populations
for (i in 1:length(pops))
{
 y=shapeithaps[,5+which(hap.pops==pops[i])]
 write.table(y,file=paste0(pathout,pops[i],"genofile.",chrno),sep="",col.names=F,row.names=F)
}
# create matrix of snps for which we have haplotypes
snps=matrix(NaN, S, 6) # same size as hapmix files, most will be left blank here
snps[,2]=chrno
snps[,4]=locs
write.table(snps, file=paste0(pathout,"snpfile.", chrno)) 


#rates=read.table(paste0(pathin,ratesfile),header=T)
#write(paste0(":sites:",nrow(rates)), paste0(pathout,"rates.",chrno))
#write(rates[,1],paste0(pathout,"rates.",chrno),append=T,ncol=nrow(rates))
#write(rates[,3],paste0(pathout,"rates.",chrno),append=T,ncol=nrow(rates)) # in centimorgans 
