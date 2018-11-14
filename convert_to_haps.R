#!/usr/bin/env Rscript
# reads in .haps files (one per chromosome) and outputs MOSAIC haps files (one per population per chromosome)
require(argparser, quiet=TRUE) 
m.args=arg_parser("script to convert shapeit2 .haps output to MOSAIC input type 1")
######################## required arguments ###############################
m.args=add_argument(m.args, "pathin", help="full path to MOSAIC inputs", type="character")
m.args=add_argument(m.args, "chrno", help="which chromosome to convert", type="integer")
m.args=add_argument(m.args, "haps_stem", help="part of .haps filename before chrno", type="character")
m.args=add_argument(m.args, "haps_end", help="part of .haps filename after chrno", type="character")
m.args=add_argument(m.args, "samples", help="name of file with details on samples", type="character")
m.args=add_argument(m.args, "pathout", help="full path to folder to store output", type="character")
m.args=add_argument(m.args, "--pops", help="donor populations to convert", default="NULL", type="character",short="-p")
argv=parse_args(m.args)
pathin=argv$pathin 
chrno=argv$chrno 
hapsfile=paste0(argv$haps_stem,chrno,argv$haps_end)
inds.data=argv$samples 
pathout=argv$pathout 
pops=strsplit(argv$pops," ")[[1]] # split the space separated group names
snps=read.table(paste0(pathin,"snpfile.",chrno),as.is=TRUE)

# read in population information; assumption is that this is correct in terms of order and size of inds in each population
allpops=read.table(paste0(pathin,inds.data),header=FALSE)
allpops=allpops[which(!is.na(match(allpops[,1],pops))),] # reduce to specified populations
hap.pops=rep(allpops[,1],each=2)
NN=length(hap.pops) # total number of haplotypes
S=nrow(snps)

# format of .haps files is RSID RSID position ref alt haps
Y=matrix(NaN,S,5+NN)
Y[,1]=snps[,1];Y[,2]=snps[,1];Y[,3]=snps[,4];Y[,4]=snps[,5];Y[,5]=snps[,6]
for (i in 1:length(pops)) {
  tmp<-scan(paste0(pathin,pops[i],"genofile.",chrno),what="character",quiet=T)
  tmp<-strsplit(tmp,"")
  allS<-length(tmp)
  N2<-length(tmp[[1]])
  Y[,5+which(hap.pops==pops[i])]=matrix(sapply(tmp, as.double), N2, allS)
}
write.table(Y,file=paste0(pathout,hapsfile),sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
