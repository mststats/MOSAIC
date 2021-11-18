#!/usr/bin/env Rscript
# script to run all code required to fit MOSAIC. Reads in data, initialises, performs thin->phase->EM cycles and outputs results.
# real example: Rscript mosaic.R Moroccan example_data/ -a 2 -n 4 -c 18:22
# simulated example: Rscript mosaic.R simulated example_data/ -a 2 -n 4 -c 18:22 -p "English Mandenka"

require(MOSAIC)
require(argparser, quiet=TRUE) 
m.args=arg_parser("run MOSAIC to model admixture and infer local ancestry without knowledge of mixing groups")
######################## required arguments ###############################
m.args=add_argument(m.args, "target", help="name of target population", type="character")
m.args=add_argument(m.args, "data", help="folder containing data", type="character")
######################## optional arguments ###############################
m.args=add_argument(m.args, "--ancestries", help="number of mixing ancestries", default=2, type="integer",short="-a")
m.args=add_argument(m.args, "--number", help="number of target individuals", default=1000, type="integer",short="-n")
m.args=add_argument(m.args, "--panels", help="listed donor groups to use; if running a simulation first #ancestries are simulated from", default="NULL", type="character",short="-p")
m.args=add_argument(m.args, "--maxcores", help="maximum number of cores to use (will grab half of all available if set to 0)", default=0, type="integer",short="-m")
m.args=add_argument(m.args, "--chromosomes", help="chromosomes as c_start:c_end", default="1:22", type="character",short="-c")
m.args=add_argument(m.args, "--rounds", help="number of inference rounds", default=5, type="integer",short="-r")
m.args=add_argument(m.args, "--GpcM", help="number of gridpoints per centiMorgan", default=60, type="integer",short="-g")
m.args=add_argument(m.args, "--maxdonors", help="maximum number of donors used per gridpoint", default=100, type="integer", short="-k")
m.args=add_argument(m.args, "--index", help="index of first individual in the target data", default=1, type="integer",short="-i")
m.args=add_argument(m.args, "--fastfiles", help="location of fast-files", default="/dev/shm/", type="character",short="-f")
m.args=add_argument(m.args, "--Ne", help="Effective Population size", default=90000, type="integer",short="-N")
m.args=add_argument(m.args, "--noEM", help="whether to perform EM inference of model parameters",flag=TRUE)
m.args=add_argument(m.args, "--nophase", help="whether to re-phase",flag=TRUE)
m.args=add_argument(m.args, "--gens", help="generations since mixing", default="NULL", type="integer")
m.args=add_argument(m.args, "--ratios", help="ratios of ancestral mixing groups", default="NULL", type="character")
m.args=add_argument(m.args, "--donors_per_group", help="maximum number of donors used per panel", default=1000, type="integer")
m.args=add_argument(m.args, "--prop", help="proportion of probability required per gridpoint", default=0.99, type="numeric")
m.args=add_argument(m.args, "--mask", help="listed groups to remove", default="NULL", type="character")
m.args=add_argument(m.args, "--noFst", help="don't calculate Fst summaries",flag=TRUE)
m.args=add_argument(m.args, "--noreturn", help="don't return result in R",flag=TRUE)
m.args=add_argument(m.args, "--noMu", help="don't update copying matrix Mu via EM",flag=TRUE)
m.args=add_argument(m.args, "--noPI", help="don't update ancestry transitions matrix PI via EM",flag=TRUE)
m.args=add_argument(m.args, "--norho", help="don't update within-ancestry recombination rates rho via EM",flag=TRUE)
m.args=add_argument(m.args, "--notheta", help="don't update emission probabilities / mutation and error rates theta via EM",flag=TRUE)
m.args=add_argument(m.args, "--noplots", help="don't create summary plots",flag=TRUE)
m.args=add_argument(m.args, "--MODE", help="diploid or haploid mode", default="DIP", type="character")
m.args=add_argument(m.args, "--singlePI", help="use a common ancestry switching matrix for all targets i.e. assume a shared history",flag=TRUE)
m.args=add_argument(m.args, "--init.rho", help="initial values for rho as space separated values", default="NULL", type="character")
m.args=add_argument(m.args, "--init.theta", help="initial values for theta as space separated values", default="NULL", type="character")
m.args=add_argument(m.args, "--init.Mu", help="initial values for Mu: Mu should be #panels rows and #A columns presented as a single space separated vector comprising the rows", default="NULL", type="character")
m.args=add_argument(m.args, "--init.PI", help="initial values for PI: one AxA matrix per admixed target, presented as a single space separated vector", default="NULL", type="character")
# simulated target and panels=c(a,b,c,d,.,.) vector will admix the first A populations and fit using the rest when running MOSAIC on simulated data

# note that flags are flipped from the default (which is FALSE) if included in the command

argv=parse_args(m.args)
target=argv$target
datasource=argv$data
chrnos=argv$chromosomes
chrnos=strsplit(chrnos,":")[[1]];chrnos=as.integer(chrnos[1]):as.integer(chrnos[2])
A=argv$ancestries
NUMI=argv$number
pops=argv$panels
mask=argv$mask
doplots=!argv$noplots # create summary plots?
doFst=!argv$noFst # calculate Fst summaries
PHASE=!argv$nophase # rephase using local ancestry model?
gens=argv$gens
ratios=argv$ratios
EM=!argv$noEM # run EM algorithm?
ffpath=argv$fastfiles
MC=argv$maxcores
return.res=!argv$noreturn #interactive() # whether to return results in a list; for use within an interactive R session 
REPS=argv$rounds
GpcM=argv$GpcM
dpg=argv$donors_per_group
max.donors=argv$maxdonors
prop.don=argv$prop
doMu=!argv$noMu # update copying matrix parameters? 
doPI=!argv$noPI # update ancestry switching parameters parameters? 
dorho=!argv$norho # update recombination w/in same ancestry parameters?  
dotheta=!argv$notheta # update error / mutation parameters? 
firstind=argv$index
Ne=argv$Ne
MODE=argv$MODE 
singlePI=argv$singlePI
init.rho=argv$init.rho 
init.theta=argv$init.theta 
init.Mu=argv$init.Mu 
init.PI=argv$init.PI 
if (!file.exists(ffpath))
  stop(paste0("requested location '", ffpath, "' for storage of fast files doesn't exist. Please use any temporary folder on your system. \n"))
if (pops=="NULL") pops=NULL
if (mask=="NULL") mask=NULL
if (gens=="NULL") gens=NULL else gens=as.numeric(gens)
if (ratios=="NULL") ratios=NULL
if (init.rho=="NULL") init.rho=NULL
if (init.theta=="NULL") init.theta=NULL
if (init.Mu=="NULL") init.Mu=NULL
if (init.PI=="NULL") init.PI=NULL
if (!is.null(pops)) pops=strsplit(pops," ")[[1]] # split the space separated group names
if (!is.null(mask)) mask=strsplit(mask," ")[[1]] # split the space separated group names
if (!is.null(ratios)) {ratios=strsplit(ratios," ")[[1]];ratios=as.numeric(ratios)}
if (!is.null(init.rho)) {init.rho=strsplit(init.rho," ")[[1]];init.rho=as.numeric(init.rho)}; 
if (!is.null(init.theta)) {init.theta=strsplit(init.theta," ")[[1]];init.theta=as.numeric(init.theta)}
if (!is.null(init.Mu)) {
  init.Mu=matrix(as.numeric(strsplit(init.Mu," ")[[1]]),ncol=A)
}
if (!is.null(init.PI)) {
  tmp=as.numeric(strsplit(init.PI, " ")[[1]]);init.PI=list();for (i in 1:NUMI) init.PI[[i]]=matrix(tmp[((i-1)*A*A+1):(i*A*A)],A)
}
if (is.null(pops) & target=="simulated") stop("Please provide at least ", A, " named populations to admix for this simulation using the -p flag", "\n")

# this function includes saving results to disk
mosaic.result=run_mosaic(target,datasource,chrnos,A,NUMI,pops,mask=mask,PLOT=doplots,doFst=doFst,PHASE=PHASE,gens=gens,ratios=ratios,EM=EM,
			 ffpath=ffpath,MC=MC,return.res=return.res,REPS=REPS,GpcM=GpcM,nl=dpg,max.donors=max.donors,prop.don=prop.don,
			 doMu=doMu,doPI=doPI,dorho=dorho,dotheta=dotheta,firstind=firstind,verbose=TRUE,Ne=Ne,
			 MODE=MODE,singlePI=singlePI,init.rho=init.rho,init.theta=init.theta,init.Mu=init.Mu,init.PI=init.PI)
