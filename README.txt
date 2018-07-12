###################################################################################################################
MOSAIC is designed to run on the linux command line, however interactive usage within R is also possible. 
To run on the command line use Rscript mosaic.R -t target -d example_data/ -A 2 -i 1 -n 2 -c 21:22
In an interactive session use run_mosaic(target,folder,chromosomes,A,n,NULL)
###################################################################################################################
There should be a folder with 4 types of input file:
1. phased haplotypes: "pop.genofile.chr" in the format #snps rows and #haps columns.
2. pop names: "sample.names" format unimportant apart from first column should have all the population names.
3. snp files: "snpfile.chr" #snps rowns and 6 columns of rsID, chr, position, distance, allele ?, allele ?. 
4. recombination map: "rates.chr" 3 rows of #sites, position, recombination rate. 
example_data contains example data for chromosomes 21 to 22 and an example run of mosaic can be done using this data via:
> Rscript mosaic.R Moroccan example_data/ 2 1 2 2 21:22
###################################################################################################################


###################################################################################################################
A folder called RESULTS is required to hold log-files (foo.out) and results (foo.RData).  
A folder called PLOTS is required to hold the plots created by default by a MOSAIC run.

You'll need to install the R package mosaicpackage_1.0.tar.gz. This is to once-off compile the C++ code that does the forward-backward, etc. 
###################################################################################################################

mosaic.R contains an example run script. It can be run from the command line using Rscript and requires 5 arguments: 
  (1) Target population (2) Folder where data is stored (3) Number of mixing groups (4) Index of first individual (5) Number of target haplotypes 
  Optional 6th argument is consecutive chromosome numbers in the form n:m
  Optional 7th argument is number of cores to use for parallel processing. 

  set (4) to 1 to run from first ind 
  set (5) to a number greater than or equal to double the number of individuals in the target pop to run all individuals at once

The main work is multiple rounds of EM->phase->thin where:
	EM runs a few iterations of the EM algorithm to update parameters. 
	phase does some rephasing (see below).
	thin finds the useful donor haplotypes at each gridpoint, based on the current EM estimations of parameters (see below).

Some other parameters that can be changed (you can set within mosaic.R): 
	1. GpcM (granularity of grid; higher is finer; #gridpoints on each chromosome=G[chr] will be GpcM*(length of chromosome in centiMorgans))
	MOSAIC currently uses 60 gidpoints per centiMorgan as a default.
###################################################################################################################

###################################################################################################################
Less important things to choose sometimes: 
	1. PHASE (logical indicating whether you want to rephase the data to avoid spurious ancestry switches)
	2. MC (number of cores to use when parallelized code is running)
	3. nl (max number of donor haps to use per donor pop)
	4. max.donors (max #haps to actually copy from at any gridpoint; defaults to #donors when set higher than this)
	5. prop.don (proportion of probability to cover; setting less than 1 means copy best hap, second best, etc. until max.donors or this is reached at each gridpoint)
	7. total (number of final EM iterations)
	8. s.total (number of EM iterations w/in rounds of EM->phase->thin)
	9. REPS (number of rounds of EM->phase->thin)
	10. s.M (controls number of MCMC iterations of rephasing to run in each round of EM->phase->thin; s.M*G[chr] is #iterations)
	11. M (controls number of MCMC iterations of rephasing to run at the end; M*G[chr] is #iterations)
Note that if you set nl or NUMA greater than available in the data they just revert to using all available. 
###################################################################################################################

###################################################################################################################
There are 4 sets of parameters inferred via EM:
	1. PI (prob. of switching between latent ancestries, including switch to same anc; LxL)
	2. rho (prob. of switching haps within each ancestry)
	3. Mu (copying matrix i.e. Mu[i,j] is  prob. of picking a hap from group i given latent ancestry j; dimension is KxL where K is #donorpops) 
	4. theta (error / mutation vector of length L; prob. of a difference b/w copied and copying haps at a locus)

Note that PI and rho will depend on GpcM; a finer grid means lower prob. of switching as we move from a gridpoint to the next one. 
###################################################################################################################

###################################################################################################################
Phase is by default updated from initial phasing. You can turn this off by setting PHASE=F
Re-phasing is done by two algorithms: "phase hunting" and MCMC. The hunting is called in each round of EM->phase->thin. 
Setting s.M=0 turns off MCMC phasing until after EM convergence and this is the default. 
###################################################################################################################

###################################################################################################################
thinning is an approximation that greatly speeds the code up when there are a large number of donor haplotypes. 
###################################################################################################################

############### After running and attaching the results from RESULTS try these plots: ################
(1) ord.Mu<-plot_Mu(Mu,MODE="copy",cexa=2) # to look at the copying matrix 
(2) ord.Mu<-plot_Mu(Mu,MODE="scaled",cexa=2) # to look at the copying matrix scaled by #donor haps in each panel. 
(3) plot_localanc(chrnos,g.loc,localanc) # cycle through all local ancestry plots (one plot per target diploid chromosome)
(4) plot_coanccurves(ind,acoancs,2,2) # plot some co-ancestry curves
###################################################################################################################

email michael.salter-townshend@ucd.ie for help
