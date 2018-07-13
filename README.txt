######  OVERVIEW   ###############################################################################################
MOSAIC is designed to run on the linux command line. 
You will need to install "argparser" and "MOSAIC" first. Try 
> Rscript mosaic.R --help 
to see how to run it. The only two arguments that must be provided is the name of the target population and the folder
in which the data are stored (see below for details).
MOSAIC may also be used in an interactive R session. First load the package
  > require(MOSAIC)
then use 
  > run_mosaic(target,folder,chromosomes,A,n) 
# where A is the number of unseen mixing groups and n is the number of target haplotypes. Specifying n larger than
the number of haplotypes in the target file results in running MOSAIC on all of them. A defaults to 2 and n to 1000.
Run 
  > plot_all_mosaic(mosaic.result,pathout="PLOTS/")
to output model plots to the folder "PLOTS/"

#######  INPUTS   #################################################################################################
There should be a folder with 4 types of input file:
1. phased haplotypes: "pop.genofile.chr" in the format #snps rows and #haps columns.
2. pop names: "sample.names" format unimportant apart from first column should have all the population names.
3. snp files: "snpfile.chr" #snps rowns and 6 columns of rsID, chr, position, distance, allele ?, allele ?. 
4. recombination map: "rates.chr" 3 rows of #sites, position, recombination rate. 

#######  EXAMPLES  ################################################################################################

example_data contains example data for chromosomes 18 to 22 and a real-data example run of mosaic can be done using this data via:
> Rscript mosaic.R Moroccan example_data/ -a 2 -n 4 -c 18:22
or equivalently in an interactive R session:
> mosaic.result=run_mosaic("Moroccan","example_data/",18:22,2,2)

To simulate admixture from the model and real data and then fit (without re-using the panels used to simulate) use
> Rscript mosaic.R simulated example_data/ -a 2 -n 4 -c 18:22 -k TRUE
or equivalently in an interactive R session:
> mosaic.result=run_mosaic("simulated","example_data/",18:22,2,2,TRUE)

User defined simulations can also be provided by specifying a vector of populations:
> Rscript mosaic.R simulated example_data/ -c 21:22 -n 2 -k "English Mandenka"

###################################################################################################################
A folder called RESULTS is required to hold log-files (foo.out) and results (foo.RData).  
A folder called PLOTS is required to hold the plots created by default by a MOSAIC run.

You'll need to install the R package MOSAIC_1.0.tar.gz. This is to once-off compile the C++ code that does the forward-backward, etc. 
###################################################################################################################

The main work is multiple rounds of EM->phase->thin where:
	EM runs a few iterations of the EM algorithm to update parameters. 
	phase does some rephasing (see below).
	thin finds the useful donor haplotypes at each gridpoint, based on the current EM estimations of parameters (see below).

########  PARAMETERS INFERRED   ###################################################################################
There are 4 sets of parameters inferred via EM:
	1. PI (prob. of switching between latent ancestries, including switch to same anc; LxL)
	2. rho (prob. of switching haps within each ancestry)
	3. Mu (copying matrix i.e. Mu[i,j] is  prob. of picking a hap from group i given latent ancestry j; dimension is KxL where K is #donorpops) 
	4. theta (error / mutation vector of length L; prob. of a difference b/w copied and copying haps at a locus)

Note that PI and rho will depend on grid granularity (GpcM); a finer grid means lower prob. of switching as we move from a gridpoint to the next one. 
###################################################################################################################

#######  CHANGING DEFAULTS   ######################################################################################
MOSAIC can be tuned to run differently by changing from the default values used. Try
> Rscript mosaic.R --help 
to view these arguments and their defaults. For example, turn re-phasing on or off (see below).

########  PHASE  ##################################################################################################
Phase is by default updated from initial phasing. You can turn this off by setting phase=FALSE
Re-phasing is done by two algorithms: "phase hunting" and MCMC. The hunting is called in each round of EM->phase->thin. 
Setting s.M=0 turns off MCMC phasing until after EM convergence and this is the default. 
###################################################################################################################

###################################################################################################################
thinning is an approximation that greatly speeds the code up when there are a large number of donor haplotypes. 
###################################################################################################################

############### After running and attaching the results from RESULTS try these plots: ################
(1) ord.Mu<-plot_Mu(Mu,alpha,NL) # to look at the copying matrix 
(2) plot_localanc(chrnos,g.loc,localanc) # cycle through all local ancestry plots (one plot per target diploid chromosome)
(3) plot_coanccurves(acoancs,dr,2,2) # plot some co-ancestry curves
###################################################################################################################

email michael.salter-townshend@ucd.ie for help
