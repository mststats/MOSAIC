# script to run all code required to fit MOSAIC. Reads in data, initialises, performs thin->phase->EM cycles and outputs results.
# example usage: Rscript mosaic.R Moroccan example_data/ 2 1 2 2 20:22
# second example: Rscript mosaic.R simulated example_data/ 
require(MOSAIC)
######################## first set some options ###############################
# FLAG use - flags for arguments
# important things to set
shargs=commandArgs(trailingOnly=TRUE) # read in arguments from the command line; 
target=shargs[1]  # e.g. Moroccan
datasource=shargs[2]; # e.g. example_data/ ; note the / at the end is crucial
L=as.integer(shargs[3]) # number of admixing groups
firstind=as.integer(shargs[4]); # which target individual to start from. If NUMA=2 then only this ind is run
NUMA=as.integer(shargs[5]) # total number of target admixed haplotypes 
# optional list of chromosomes on which to run mosaic; defaults to 1:22 if not supplied
chrnos=shargs[6]
if (!is.na(chrnos)) {chrnos=strsplit(chrnos,":")[[1]];chrnos=as.integer(chrnos[1]):as.integer(chrnos[2])}
if (is.na(chrnos)) chrnos=1:22 # which chromosomes to run on
MC=as.integer(shargs[7]) # optional number of cores to use for parallelized code; will grab half of all available if not supplied
ANC=NULL; # no a-priori knowledge of which panels to use for which ancestries

#comment this out for demo
#chrnos=21:22;firstind=1;NUMA=2;L=2;datasource="example_data/";target="Moroccan";ANC=NULL
#comment this out for simulated data demo
chrnos=21:22;firstind=1;NUMA=2;L=2;datasource="example_data/";target="simulated";ANC=TRUE

nchrno=length(chrnos) # number of chromosomes for these target haplotypes
ffpath="/dev/shm/" # location of fast-files

# the rest are mostly used in debugging, etc
verbose=T # print certain statements of progress as algorithm runs?
EM=T # run EM algorithm?
doMu=T # update copying matrix parameters?
doPI=T # update ancestry switching parameters parameters?
dorho=T # update recombination w/in same ancestry parameters? 
dotheta=T # update error / mutation parameters?
PLOT=F
return.res=interactive() # whether to return results in a list; for use within an interactive R session

# this function includes saving results to disk
mosaic.result=run_mosaic(ANC,chrnos,datasource,doMu,doPI,dorho,dotheta,EM,ffpath,firstind,
			 L,MC,nchrno,NUMA,PLOT,target,verbose,return.res) 

cat("Expected r-squared (genomewide):", dip_expected_fr2(mosaic.result$localanc),"\n")
if (target=="simulated") 
  cat("Actual r-squared (genomewide):", dip_fr2(mosaic.result$localanc,mosaic.result$g.true_anc),"\n")

plot_all_mosaic(mosaic.result,pathout="PLOTS/")
