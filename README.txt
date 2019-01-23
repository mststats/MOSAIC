######  OVERVIEW   ###############################################################################################
MOSAIC is designed to run on the linux command line. 
You can install MOSAIC by running the command
> R CMD INSTALL -l rlib MOSAIC_1.2.tar.gz
where rlib is the location where you'd like to install the R package (e.g. /usr/lib/R/site-library/)
or use 
> install.packages("MOSAIC_1.2.tar.gz")
within R

This is once-off and compiles all C++ code for the Hidden Markov Models underlying MOSAIC. 
To run MOSAIC on the command line it is also necessary to install "argparser" (available from CRAN). 

#######  HELP    #################################################################################################
Try 
> Rscript mosaic.R --help 
to see how to use MOSAIC. The only two arguments that must be provided are:
  (1) The name of the target admixed population 
  (2) The folder in which the data are stored (see below for details).
MOSAIC may also be used in an interactive R session. First load the package
> require(MOSAIC)
then use 
> run_mosaic(target,folder,chromosomes,A,n) 
# where A is the number of unseen mixing groups and n is the number of target haplotypes. Specifying n larger than
the number of haplotypes in the target file results in running MOSAIC on all of them. A defaults to 2 and n to 1000.

#######  INPUTS   ################################################################################################
There should be a folder with 4 types of input file:
1. phased haplotypes: "pop.genofile.chr" in the format #snps rows and #haps columns.
2. pop names: "sample.names" format unimportant apart from first column should have all the population names.
3. snp files: "snpfile.chr" #snps rowns and 6 columns of rsID, chr, distance, position, allele ?, allele ?. 
4. recombination map: "rates.chr" 3 rows of #sites, position, cumulative recombination rate (in centiMorgans). 

#######  PLOTS   #################################################################################################
In R, after loading the results (including localanc_foo file) of a MOSAIC run (stored in MOSAIC_RESULTS by default) use:
> plot_all_mosaic(pathout="MOSAIC_PLOTS/", EM)
to output default plots to the folder "MOSAIC_PLOTS/"
You can also use:
(1) > ord.Mu<-plot_Mu(Mu,alpha,NL) # to look at the copying matrix 
(2) > plot_coanccurves(acoancs,dr) # plot some co-ancestry curves that are used to infer event timings
(3) > plot_localanc(chrnos,g.loc,localanc) # cycles through all local ancestry plots (one plot per target diploid chromosome)
(4) > plot_loglike(extract_log(logfile)) # plots the model fit across iterations of thin/phase/EM (if EM is on)

#######  EXAMPLES  ###############################################################################################
The "example_data" folder contains example data for chromosomes 18 to 22 and a real-data example run of mosaic can be done using:
> Rscript mosaic.R Moroccan example_data/ -a 2 -n 4 -c 18:22
or equivalently in an interactive R session:
> mosaic.result=run_mosaic("Moroccan","example_data/",18:22,2,4)

User defined simulations can also be provided by specifying a vector of populations:
> Rscript mosaic.R simulated example_data/ -c 18:22 -n 2 -p "English Mandenka"
or equivalently in an interactive R session:
> mosaic.result=run_mosaic("simulated","example_data/",18:22,2,4,c("English","Mandenka"))
Note that additional groups will be used as the donor panels

##### OUTPUTS ####################################################################################################
A folder called MOSAIC_RESULTS is required to hold log-files (foo.out) and results (foo.RData).  
A folder called MOSAIC_PLOTS is required to hold the plots created by default by a MOSAIC run.
A folder called FREQS is required to hold the frequencies used to compute Fst statistics if required.
##################################################################################################################

The main work is multiple rounds of EM->phase->thin where:
	EM runs a few iterations of the EM algorithm to update parameters. 
	phase does some rephasing (see below).
	thin finds the useful donor haplotypes at each gridpoint, based on the estimated parameters (see below).

########  PARAMETERS INFERRED   ##################################################################################
There are 4 sets of parameters inferred via EM:
	1. PI (prob. of switching between latent ancestries, including switch to same anc; AxA)
	2. rho (prob. of switching haps within each ancestry)
	3. Mu (copying matrix; Mu[i,j] is  prob. of donor from group i given ancestry j; KxA where K is #donorpops) 
	4. theta (error / mutation; vector length A; prob. of a difference b/w copied and copying haps at a locus)

Note that PI and rho will depend on grid granularity (GpcM); finer grid => lower prob. of switching between gridpoints.
##################################################################################################################

#######  CHANGING DEFAULTS   #####################################################################################
MOSAIC can be tuned to run differently by changing from the default values used. Try
> Rscript mosaic.R --help 
to view these arguments and their defaults. For example, turn re-phasing on or off (see below).

########  PHASE  #################################################################################################
Phase is by default updated from initial phasing. You can turn this off by including -nophase in the command. 
Re-phasing is done by two algorithms: "phase hunting" and MCMC. The hunting is called in each round of EM->phase->thin. 
Setting s.M=0 turns off MCMC phasing until after EM convergence and this is the default. 
##################################################################################################################

##################################################################################################################
thinning is an approximation that greatly speeds the code up when there are a large number of donor haplotypes. 
##################################################################################################################

##################################################################################################################

email michael.salter-townshend@ucd.ie for help
