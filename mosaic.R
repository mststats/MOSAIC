#!/usr/bin/env Rscript
# script to run all code required to fit MOSAIC. Reads in data, initialises, performs thin->phase->EM cycles and outputs results.
# real example: Rscript mosaic.R -t Moroccan -d example_data/ -a 2 -n 4 -c 18:22
# simulated example: Rscript mosaic.R -t simulated -d example_data/ -a 2 -n 4 -c 18:22 -k TRUE

require(MOSAIC)
require(argparser, quiet=TRUE) 
m.args=arg_parser("run MOSAIC to model admixture and infer local ancestry without knowledge of mixing groups")
######################## required arguments ###############################
m.args=add_argument(m.args, "target", help="name of target population", type="character")
m.args=add_argument(m.args, "data", help="folder containing data", type="character")
######################## optional arguments ###############################
m.args=add_argument(m.args, "--ancestries", help="number of mixing ancestries", default=2, type="integer",short="-a")
m.args=add_argument(m.args, "--number", help="number of target haplotypes (double the number of individuals assuming diploid)", default=1000, type="integer",short="-n")
m.args=add_argument(m.args, "--maxcores", help="maximum number of cores to use (will grab half of all available if set to 0)", default=0, type="integer",short="-m")
m.args=add_argument(m.args, "--noEM", help="whether to perform EM inference of model parameters",flag=TRUE,short="-noEM")
m.args=add_argument(m.args, "--nophase", help="whether to re-phase",flag=TRUE,short="-nophase")
m.args=add_argument(m.args, "--gens", help="generations since mixing", default=0, type="integer",short="-gens")
m.args=add_argument(m.args, "--ratios", help="ratios of ancestral mixing groups", default="NULL", type="character",short="-ratios")
m.args=add_argument(m.args, "--chromosomes", help="chromosomes as c_start:c_end", default="1:22", type="character",short="-c")
m.args=add_argument(m.args, "--rounds", help="number of inference rounds", default=5, type="integer",short="-r")
m.args=add_argument(m.args, "--GpcM", help="number of gridpoints per centiMorgan", default=60, type="integer",short="-g")
m.args=add_argument(m.args, "--donors_per_group", help="maximum number of donors used per panel", default=1000, type="integer",short="-dpg")
m.args=add_argument(m.args, "--maxdonors", help="maximum number of donors used per gridpoint", default=100, type="integer")
m.args=add_argument(m.args, "--prop", help="proportion of probability required per gridpoint", default=0.99, type="numeric",short="-p")
m.args=add_argument(m.args, "--index", help="index of first individual in the target data", default=1, type="integer",short="-i")
m.args=add_argument(m.args, "--fastfiles", help="location of fast-files", default="/dev/shm/", type="character",short="-f")
m.args=add_argument(m.args, "--known", help="a-priori mixing groups", default="NULL", type="character",short="-k")
# known=NULL is no a-priori knowledge of which panels to use for which ancestries and all panels are used
# known=TRUE will run a pre-programmed simulation 
# known=c(a,b,c,d,.,.) vector will admix the first A populations and fit using the rest when running MOSAIC on simulated data

# note that flags are flipped from the default (which is FALSE) if included in the command

argv=parse_args(m.args)
target=argv$target
datasource=argv$data
A=argv$ancestries
firstind=argv$index
NUMA=argv$number
chrnos=argv$chromosomes
chrnos=strsplit(chrnos,":")[[1]];chrnos=as.integer(chrnos[1]):as.integer(chrnos[2])
MC=argv$maxcores
ANC=argv$known
REPS=argv$rounds
GpcM=argv$GpcM
EM=!argv$noEM # run EM algorithm?
PHASE=!argv$nophase # rephase using local ancestry model?
gens=argv$gens
ratios=argv$ratios
dpg=argv$donors_per_group
max.donors=argv$maxdonors
prop.don=argv$prop
ffpath=argv$fastfiles
if (ANC=="TRUE") ANC=TRUE
if (ANC=="NULL") ANC=NULL
if (ratios=="NULL") ratios=NULL
if (!is.null(ratios)) {ratios=strsplit(ratios," ")[[1]];ratios=as.numeric(ratios)}
if (is.null(ANC) & target=="simulated") {
  warning('use --known TRUE or --known "vector of populations" when running a simulation\n
	  setting --known TRUE', immediate.=T)
}

# the rest are mostly used in debugging, etc
verbose=T # print certain statements of progress as algorithm runs?
doMu=T # update copying matrix parameters?
doPI=T # update ancestry switching parameters parameters?
dorho=T # update recombination w/in same ancestry parameters? 
dotheta=T # update error / mutation parameters?
return.res=TRUE #interactive() # whether to return results in a list; for use within an interactive R session

# this function includes saving results to disk
mosaic.result=run_mosaic(target,datasource,chrnos,A,NUMA,ANC,REPS=REPS,GpcM=GpcM,PHASE=PHASE,nl=dpg,max.donors=max.donors,prop.don=prop.don,
			 return.res=return.res,ffpath=ffpath,doMu=doMu,doPI=doPI,dorho=dorho,dotheta=dotheta,EM=EM,gens=gens,ratios=ratios,
			 firstind=firstind,MC=MC,verbose=verbose) 

plot_all_mosaic(mosaic.result,pathout="MOSAIC_PLOTS/", EM=EM)
