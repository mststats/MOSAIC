#!/usr/bin/env Rscript
# script to run all code required to fit MOSAIC. Reads in data, initialises, performs thin->phase->EM cycles and outputs results.
# real example: Rscript mosaic.R -t Moroccan -d example_data/ -a 2 -n 4 -c 18:22
# simulated example: Rscript mosaic.R -t simulated -d example_data/ -a 2 -n 4 -c 18:22 -p "English Mandenka"

require(MOSAIC)
require(argparser, quiet=TRUE) 
m.args=arg_parser("run MOSAIC to model admixture and infer local ancestry without knowledge of mixing groups")
######################## required arguments ###############################
m.args=add_argument(m.args, "target", help="name of target population", type="character")
m.args=add_argument(m.args, "data", help="folder containing data", type="character")
######################## optional arguments ###############################
m.args=add_argument(m.args, "--ancestries", help="number of mixing ancestries", default=2, type="integer",short="-a")
m.args=add_argument(m.args, "--number", help="number of target individuals", default=1000, type="integer",short="-n")
m.args=add_argument(m.args, "--maxcores", help="maximum number of cores to use (will grab half of all available if set to 0)", default=0, type="integer",short="-m")
m.args=add_argument(m.args, "--noEM", help="whether to perform EM inference of model parameters",flag=TRUE,short="-noEM")
m.args=add_argument(m.args, "--nophase", help="whether to re-phase",flag=TRUE,short="-nophase")
m.args=add_argument(m.args, "--gens", help="generations since mixing", default=0, type="integer",short="-gens")
m.args=add_argument(m.args, "--ratios", help="ratios of ancestral mixing groups", default="NULL", type="character",short="-ratios")
m.args=add_argument(m.args, "--chromosomes", help="chromosomes as c_start:c_end", default="1:22", type="character",short="-c")
m.args=add_argument(m.args, "--rounds", help="number of inference rounds", default=5, type="integer",short="-r")
m.args=add_argument(m.args, "--GpcM", help="number of gridpoints per centiMorgan", default=60, type="integer",short="-g")
m.args=add_argument(m.args, "--donors_per_group", help="maximum number of donors used per panel", default=1000, type="integer",short="-dpg")
m.args=add_argument(m.args, "--maxdonors", help="maximum number of donors used per gridpoint", default=100, type="integer", short="-k")
m.args=add_argument(m.args, "--prop", help="proportion of probability required per gridpoint", default=0.99, type="numeric",short="-prop")
m.args=add_argument(m.args, "--index", help="index of first individual in the target data", default=1, type="integer",short="-i")
m.args=add_argument(m.args, "--fastfiles", help="location of fast-files", default="/dev/shm/", type="character",short="-f")
m.args=add_argument(m.args, "--panels", help="listed donor groups to use; if running a simulation first #ancestries are simulated from", default="NULL", type="character",short="-p")
m.args=add_argument(m.args, "--mask", help="listed groups to remove", default="NULL", type="character",short="-mask")
m.args=add_argument(m.args, "--Ne", help="Effective Population size", default=90000, type="integer",short="-Ne")
m.args=add_argument(m.args, "--Fst", help="calculate Fst summaries",flag=TRUE,short="-Fst")
m.args=add_argument(m.args, "--plots", help="create summary plots",flag=TRUE,short="-plots")
# simulated target and panels=c(a,b,c,d,.,.) vector will admix the first A populations and fit using the rest when running MOSAIC on simulated data

# note that flags are flipped from the default (which is FALSE) if included in the command

argv=parse_args(m.args)
target=argv$target
datasource=argv$data
A=argv$ancestries
firstind=argv$index
NUMI=argv$number
chrnos=argv$chromosomes
chrnos=strsplit(chrnos,":")[[1]];chrnos=as.integer(chrnos[1]):as.integer(chrnos[2])
MC=argv$maxcores
pops=argv$panels
mask=argv$mask
REPS=argv$rounds
GpcM=argv$GpcM
EM=!argv$noEM # run EM algorithm?
PHASE=!argv$nophase # rephase using local ancestry model?
doFst=!argv$Fst # calculate Fst summaries
doplots=!argv$plots # create summary plots?
gens=argv$gens
Ne=argv$Ne
ratios=argv$ratios
dpg=argv$donors_per_group
max.donors=argv$maxdonors
prop.don=argv$prop
ffpath=argv$fastfiles
if (pops=="NULL") pops=NULL
if (mask=="NULL") mask=NULL
if (ratios=="NULL") ratios=NULL
if (!is.null(pops)) pops=strsplit(pops," ")[[1]] # split the space separated group names
if (!is.null(mask)) mask=strsplit(mask," ")[[1]] # split the space separated group names
if (!is.null(ratios)) {ratios=strsplit(ratios," ")[[1]];ratios=as.numeric(ratios)}
if (is.null(pops) & target=="simulated") stop("Please provide at least ", A, " named populations to admix for this simulation using the -p flag", "\n")

# the rest are mostly used in debugging, etc
verbose=T # print certain statements of progress as algorithm runs?
doMu=T # update copying matrix parameters?
doPI=T # update ancestry switching parameters parameters?
dorho=T # update recombination w/in same ancestry parameters? 
dotheta=T # update error / mutation parameters?
return.res=TRUE #interactive() # whether to return results in a list; for use within an interactive R session

# this function includes saving results to disk
mosaic.result=run_mosaic(target,datasource,chrnos,A,NUMI,pops,mask=mask,PLOT=doplots,doFst=doFst,PHASE=PHASE,gens=gens,ratios=ratios,EM=EM,
			 ffpath=ffpath,MC=MC,return.res=return.res,REPS=REPS,GpcM=GpcM,nl=dpg,max.donors=max.donors,prop.don=prop.don,
			 doMu=doMu,doPI=doPI,dorho=dorho,dotheta=dotheta,firstind=firstind,verbose=verbose,Ne=Ne)
