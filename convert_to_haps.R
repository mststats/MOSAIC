#!/usr/bin/env Rscript
# reads in MOSAIC haps files (one per population per chromosome) and outputs .haps files (one per chromosome) 
require(argparser, quiet=TRUE) 
require(LaF, quiet=TRUE) 
m.args=arg_parser("script to convert MOSAIC input files to shapeit2 .haps files")
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
sampsfile=paste0(argv$haps_stem,chrno,inds.data)
pathout=argv$pathout 
pops=argv$pops
if (pops=="NULL") pops=NULL
if (!is.null(pops)) pops=strsplit(pops," ")[[1]] # split the space separated group names
snps=read.table(file.path(pathin,paste0("snpfile.",chrno)),as.is=TRUE)

# read in population information; assumption is that this is correct in terms of order and size of inds in each population
allpops=read.table(file.path(pathin,sampsfile),header=FALSE)
if (!is.null(pops))
  allpops=allpops[which(!is.na(match(allpops[,1],pops))),] # reduce to specified populations
pops=unique(allpops[,1])
hap.pops=rep(allpops[,1],each=2)
NN=length(hap.pops) # total number of haplotypes
S=nrow(snps)

# format of .haps files is RSID RSID position ref alt haps
Y=matrix(NaN,S,5+NN)
Y[,1]=chrno;Y[,2]=snps[,1];Y[,3]=snps[,4];Y[,4]=snps[,5];Y[,5]=snps[,6]
multigenos=list()
for (i in 1:length(pops)) {
  tmp<-scan(file.path(pathin,paste0(pops[i],"genofile.",chrno)),what="character",quiet=TRUE)
  tmp<-strsplit(tmp,"")
  allS<-length(tmp)
  N2<-length(tmp[[1]])
  tmpfilename<-file.path(pathin,paste0(pops[i],"genofile.",chrno))
  tmp<-scan(tmpfilename,what="character",quiet=TRUE,nlines=1)
  N2<-nchar(tmp)
  multigenos[[i]]=laf_open_fwf(tmpfilename, column_widths=rep(1,N2),column_types=rep("character",N2))
  Y[,5+which(hap.pops==pops[i])]=matrix(sapply(multigenos[[i]][,], as.double), allS, N2)
}
write.table(Y,file=file.path(pathout,hapsfile),sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)

# format of .samples file is ID_1, ID_2, missing
allpops=allpops[,1:3] # only first 3 columns
allpops=rbind(c("ID_1", "ID_2", "missing"), rep(0,3), allpops)
write.table(allpops ,file=file.path(pathout,sampsfile),sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
