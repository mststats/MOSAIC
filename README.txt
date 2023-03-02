######  OVERVIEW   ###############################################################################################
MOSAIC is designed to run on the linux command line or within R. 
You can install MOSAIC by running the command
> R CMD INSTALL -l rlib MOSAIC_1.5.0.tar.gz
where rlib is the location where you'd like to install the R package (e.g. /usr/lib/R/site-library/)
or use 
> install.packages("MOSAIC_1.5.0.tar.gz")
within R

This is once-off and compiles all C++ code for the Hidden Markov Models underlying MOSAIC. 
To run MOSAIC on the command line it is also necessary to install "argparser" (available from CRAN). 

#######  HELP    #################################################################################################
There is a vignette with walkthrough examples available by calling
> vignette("MOSAIC")
within R or download it at www.maths.ucd.ie/~mst/MOSAIC/MOSAIC.pdf 

Try 
$ Rscript mosaic.R --help 
to see how to use MOSAIC. The only two arguments that must be provided are:
  (1) The name of the target admixed population 
  (2) The folder in which the data are stored (see below for details).
MOSAIC may also be used in an interactive R session. First load the package
> require(MOSAIC)
then use 
> run_mosaic(target,folder,chromosomes,a,n) 
# where "a" is the number of unseen mixing groups and n is the number of target individuals. Specifying n larger than
the number of individuals in the target file results in running MOSAIC on all of them. a defaults to 2 and n to 1000.

#######  INPUTS   ################################################################################################
There should be a folder with 4 types of input file:
1. phased haplotypes: "pop.genofile.chr" in the format #snps rows and #haps columns.
2. pop names: "sample.names" format unimportant apart from first column should have all the population names.
3. snp files: "snpfile.chr" #snps rows and 6 columns of rsID, chr, distance, position, allele ?, allele ?. 
4. recombination map: "rates.chr" 3 rows of #sites, position, cumulative recombination rate (in centiMorgans). 

#######  PLOTS   #################################################################################################
In R, after loading the results (including localanc_foo file) of a MOSAIC run (stored in MOSAIC_RESULTS by default) use:
> plot_all_mosaic(pathout="MOSAIC_RESULTS/",target) # note that this is already run automatically by default in run_mosaic
to output default plots to the folder "MOSAIC_RESULTS/"
You can also individually use:
(1) > ord.Mu=plot_Mu(Mu,alpha,NL) # to look at the copying matrix 
(2) > plot_coanccurves(acoancs,dr) # plot some co-ancestry curves that are used to infer event timings
(3) > plot_localanc(chrnos,g.loc,localanc) # cycles through all local ancestry plots (one plot per target diploid chromosome)
(4) > plot_loglike(extract_log(logfile)) # plots the model fit across iterations of thin/phase/EM (if EM is on)

#######  converting local ancestry to SNP positions ##############################################################
MOSAIC outputs local ancestry estimates along the genome to a file called localanc_$target_$a-way_$indfirst-$indlast_$chr1-$chrlast_etc.RData 
in the folder MOSAIC_RESULTS/ (by default). Once you load this file and the file of the same name but without the localanc_ prefix in R you
can convert to local ancestry at your SNP positions on the first chromosome you analysed using:
> ans=grid_to_pos(localanc[[1]],loci,g.loc[[1]]) # where loci are the SNP positions you'd like to map back to.


#######  EXAMPLES  ###############################################################################################
The "extdata" folder contains example data for chromosomes 18 to 22 and a real-data example run of mosaic can be done using:
$ Rscript mosaic.R Moroccan extdata/ -a 2 -n 2 -c 18:22
Note that the extdata folder is in MOSAIC/inst/ in the source code but at installation-library/MOSAIC/ after installation of MOSAIC
where installation-library is the location on your system where you have installed MOSAIC. 
Alternatively, in an interactive R session you can use:
> fpath=system.file("extdata", package="MOSAIC")
> mosaic.result=run_mosaic("Moroccan", fpath, 18:22, 2, 2)

User defined simulations can also be provided by specifying a vector of populations:
$ Rscript mosaic.R simulated extdata/ -c 18:22 -n 3 -p "English Mandenka"
or equivalently in an interactive R session:
> fpath=system.file("extdata", package="MOSAIC")
> mosaic.result=run_mosaic("simulated", fpath, 18:22, 2, 2, c("English","Mandenka"))
Note that additional groups will be used as the donor panels but can also be specified manually as follows:
> mosaic.result=run_mosaic("simulated", fpath, 18:22, 2, 2, c("English","Mandenka", "French", "Yoruba"))

##### OUTPUTS ####################################################################################################
A folder called MOSAIC_RESULTS storing log-files (foo.out) and results (foo.RData).  This also stores the PDF plots created by MOSAIC. 
A folder called FREQS within the results folder containing the frequencies used to compute Fst statistics if required.
##################################################################################################################

The main work is multiple rounds of EM->phase->thin where:
	EM runs a few iterations of the EM algorithm to update parameters. 
	phase does some rephasing (see below).
	thin finds the useful donor haplotypes at each gridpoint, based on the estimated parameters (see below).

########  PARAMETERS INFERRED   ##################################################################################
There are 4 sets of parameters inferred via EM:
	1. Pi (prob. of switching between latent ancestries, including switch to same anc; AxA)
	2. rho (prob. of switching haps within each ancestry)
	3. Mu (copying matrix; Mu[k,a] is  prob. of donor from group k given ancestry a; KxA where K is #donorpops) 
	4. theta (error / mutation; vector length A; prob. of a difference b/w copied and copying haps at a locus)

Note that Pi and rho will depend on grid granularity (GpcM); finer grid => lower prob. of switching between gridpoints.
##################################################################################################################

#######  CHANGING DEFAULTS   #####################################################################################
MOSAIC can be tuned to run differently by changing from the default values used. Try
$ Rscript mosaic.R --help 
to view these arguments and their defaults. For example, turn re-phasing on or off (see below).
You can also tell MOSAIC which of the above 4 sets of parameters to update by setting a flag --notheta on the command line
to stay with the initial estimates of theta for example. In the interactive mode you can set dotheta=FALSE as an argument 
to run_mosaic() within R. 

Specific starting values for the 4 sets of parameters can be set using init.PI, init.rho, init.Mu, and init.theta. 
These could for example be final estimates from a previous MOSAIC run. You can either supply these as arguments to
run_mosaic() in R:
> fpath=system.file("extdata", package="MOSAIC")
> mosaic.result=run_mosaic("Moroccan", fpath, 18:21,A=2,NUMI=4) # run MOSAIC on chromosomes 18 to 21 only to estimate parameters
> load("MOSAIC_RESULTS/Moroccan_2way_1-4_18-21_162_60_0.99_100.RData") # load the results. 
# Note that you could also use attach(mosaic.result) for this, if still in the same R session 
# or read parameters from the log file using paras=extract_paras(extract_log(EMlogfile)) 
> fpath=system.file("extdata", package="MOSAIC")
> mosaic.result22=run_mosaic("Moroccan", fpath, 22, A=2, NUMI=4, EM=F, init.PI=PI, init.rho=rho,init. Mu=Mu, init.theta=theta) 
# above line runs MOSAIC on chromosome 22 using the parameter estimates from the previous run as initial values and don't update via EM

On the command line the parameters can be provided as strings containing values separated by spaces. For example:
$ Rscript mosaic.R Moroccan extdata/ -a 2 -n 2 -c 22:22 --init.rho "0.09 0.09"  -p "English Mandenka"
or (more usefully) you can extract values from the log file of a previous run and use that directly:
$ logfile=MOSAIC_RESULTS/Moroccan_2way_1-4_18-21_162_60_2021_07_14_10:29:37_EMlog.out # replace with the log file generated by your MOSAIC run above
$ PI_command='cat\(unlist\(MOSAIC::extract_paras\(MOSAIC::extract_log\(\"$logfile\"\)\)\$PI\)\)'
$ PI=`eval Rscript -e $PI_command | cut -d " " -f 1-`
$ Mu_command='cat\(c\(MOSAIC::extract_paras\(MOSAIC::extract_log\(\"$logfile\"\)\)\$Mu\)\)'
$ Mu=`eval Rscript -e $Mu_command | cut -d " " -f 1-`
$ rho_command='cat\(c\(MOSAIC::extract_paras\(MOSAIC::extract_log\(\"$logfile\"\)\)\$rho\)\)'
$ rho=`eval Rscript -e $rho_command | cut -d " " -f 1-`
$ theta_command='cat\(c\(MOSAIC::extract_paras\(MOSAIC::extract_log\(\"$logfile\"\)\)\$theta\)\)'
$ theta=`eval Rscript -e $theta_command | cut -d " " -f 1-`
$ Rscript mosaic.R Moroccan extdata/ -c 22:22 -a 2 -n 4 --noEM --init.PI "$PI" --init.Mu "$Mu" --init.rho "$rho" --init.theta "$theta"

MOSAIC uses the ff package to store and manipulate "fast files". The default location for these when running in
an interactive R session is set by a call to tempdir(), whereas when running using Rscript on the command line it uses 
/dev/shm/ (which assumes a linux installation). These defaults can be changed by using the -f option for the command
line or the ffpath argument when using run_mosaic() within R. 

########  PHASE  #################################################################################################
Phase is by default updated from initial phasing. You can turn this off by including --nophase in the command. 
Re-phasing is done by two algorithms: "phase hunting" and MCMC. The hunting is called in each round of EM->phase->thin. 
Setting s.M=0 turns off MCMC phasing until after EM convergence and this is the default. 
##################################################################################################################

##################################################################################################################
Thinning is an approximation that greatly speeds the code up when there are a large number of donor haplotypes by 
using only the most likely to copied from ones (MOSAIC finds these for you). 
You can adjust the number used by changing prop.don and / or max.donors. 
##################################################################################################################

##################################################################################################################

email michael.salter-townshend@ucd.ie for help
