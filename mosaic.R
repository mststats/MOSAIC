# script to run all code required to fit MOSAIC. Reads in data, initialises, performs thin->phase->EM cycles and outputs results.
# example usage: Rscript mosaic.R SanKhomani HGDP/ 4 1 60 16
require(MOSAIC)
######################## first set some options ###############################
# important things to set
shargs<-commandArgs(trailingOnly=TRUE) # read in arguments from the command line; 
target=shargs[1]  # e.g. Moroccan
datasource=shargs[2]; # e.g. example_data/ ; note the / at the end is crucial
L=as.integer(shargs[3]) # number of admixing groups
firstind<-as.integer(shargs[4]); # which target individual to start from. If NUMA=2 then only this ind is run
NUMA=as.integer(shargs[5]) # total number of target admixed haplotypes 

MC=as.integer(shargs[6]) # optional number of cores to use for parallelized code; will grab half of all available if not supplied
chrnos=1:22 # which chromosomes to run on

# FLAG
#chrnos=21:22;firstind=1;NUMA=2;L=2;datasource="example_data/";target="Moroccan";ANC=NULL
chrnos=20:22;firstind=1;NUMA=2;L=2;datasource="HGDP/";target="simulated";RPE=0.0;ANC=T;

nchrno=length(chrnos) # number of chromosomes for these target haplotypes
ffpath="/dev/shm/" # location of fast-files

# the rest are mostly used in debugging, etc
#ANC=NULL; # no a-priori knowledge of which panels to use for which ancestries
verbose=T # print certain statements of progress as algorithm runs?
EM=T # run EM algorithm?
doMu=T # update copying matrix parameters?
doPI=T # update ancestry switching parameters parameters?
dorho=T # update recombination w/in same ancestry parameters? 
dotheta=T # update error / mutation parameters?
PLOT=F
return.res=TRUE # whether to return results in a list; for use within an interactive R session

# this function includes saving results to disk
mosaic.result=run_mosaic(ANC,chrnos,datasource,doMu,doPI,dorho,dotheta,EM,ffpath,firstind,
			 L,MC,nchrno,NUMA,PLOT,target,verbose,return.res) 

cat("Expected r-squared (genomewide):", dip_expected_fr2(mosaic.result$localanc),"\n")
if (target=="simulated") 
  cat("Actual r-squared (genomewide):", dip_fr2(mosaic.result$localanc,mosaic.result$g.true_anc),"\n")

# get de-phased local ancestry
#flocalanc=phase_localanc(mosaic.result$localanc,mosaic.result$final.flips)
#for (panel in rownames(mosaic.result$Mu))
#  summarise_panels(panel,datasource,mosaic.result$chrnos)
plot_all(mosaic.result,pathout="PLOTS/","FREQS/")
