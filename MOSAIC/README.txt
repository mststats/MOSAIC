There should be a folder with 4 types of input file:
1. phased haplotypes: "pop.genofile.chr" in the format #snps rows and #haps columns.
2. pop names: "sample.names" format unimportant apart from first column should have all the population names.
3. snp files: "snpfile.chr" #snps rowns and 6 columns of rsID, chr, position, distance, allele ?, allele ?. 
4. recombination map: "rates.chr" 3 rows of #sites, position, recombination rate. 


A folder called RESULTS needs to be created to hold log-files (foo.out) and results (foo.RData).  
Also folders called LOGS and PLOTS. 

You'll need to install the mosaicpackage_1.0.tar.gz. This is to once-off compile the C++ code that does the forward-backward, etc. 

run.R then contains an example script to try. It can be run from the command line using Rscript and requires 5 arguments: 
  (1) target pop (2) datasource (3) #ancestries (4) firstind (5) Number of target haplotypes 
  Optional 6th arguments in #cores to use for parallel processing. 

  set (4) to 1 to run from first ind 
  set (5) to a number greater than or equal to the number of inds in the target pop to run all inds
  You can also hard code these into run.R if you prefer. 

The main work is multiple rounds of EM->phase->thin where:
	EM runs a few iterations of the EM algorithm to update parameters. 
	phase does some rephasing (see below).
	thin finds the useful donor haplotypes at each gridpoint, based on the current EM estimations of parameters (see below).

Important things to choose (set within run.R): 
	1. chrnos (which chromosomes to run on; full genome is via e.g. chrnos=22:1)
	2. GpcM (granularity of grid; higher is finer; #gridpoints on each chromosome=G[chr] will be GpcM*(length of chromosome in centiMorgans))
	We're currently using 60 gidpoints per centiMorgan. 

Less important things to choose sometimes: 
	1. PHASE (logical indicating whether you want to rephase the data to avoid spurious ancestry switches)
	2. MC (number of cores to use when parallelized code is running)
	3. nl (max number of donor haps to use per donor pop)
	4. max.donors (max #haps to actually copy from at any gridpoint; defaults to #donors when set higher than this)
	5. prop.don (proportion of probability to cover; setting less than 1 means copy best hap, second best, etc. until max.donors or this is reached at each gridpoint)
	6. drop.threshold (used in proposing donor groups to drop either b/c they are not good at discriminating b/w the ancestries or they are just never copied from)
	7. total (number of final EM iterations)
	8. s.total (number of EM iterations w/in rounds of EM->phase->thin)
	9. REPS (number of rounds of EM->phase->thin)
	10. s.M (controls number of MCMC iterations of re-phasing to run in each round of EM->phase->thin; s.M*G[chr] is #iterations)
	11. M (controls number of MCMC iterations of re-phasing to run at the end; M*G[chr] is #iterations)
Note that if you set nl or NUMA greater than available in the data they just revert to using all available. 

4 sets of parameters inferred via EM:
	1. Q (prob. of switching between latent ancestries, including switch to same anc; LxL)
	2. rho (prob. of switching haps within each ancestry)
	3. Mu (copying matrix i.e. Mu[i,j] is  prob. of picking a hap from group i given latent ancestry j; dimension is KxL where K is #donorpops) 
	4. theta (error / mutation vector of length L; prob. of a difference b/w copied and copying haps at a locus)

Note that Q and rho will depend of GpcM; a finer grid means lower prob. of switching as we move from a gridpoint to the next one. 

Phase is by default updated from initial phasing. You can turn this off by setting PHASE=F
Re-phasing is done by two algorithms: "phase hunting" and MCMC. The hunting is called in each round of EM->phase->thin. 
Setting s.M=0 would turn off MCMC phasing until after EM convergence. 

thinning is an approximation that greatly speeds the code up when there are a large number of donor haplotypes. 

############### After running and loading the results try these plots: ################
(1) source("plot_funcs.R") # to load various plotting functions
(2) ord.Mu<-plot_Mu(Mu,MODE="copy",cexa=2) # to look at the copying matrix 
(3) ord.Mu<-plot_Mu(Mu,MODE="scaled",cexa=2) # to look at the copying matrix scaled by #donor haps in each panel. 
(3) PNG=F;MODE="DIP";showgradient=F;source("barplot_localanc.R") # cycle through all local ancestry plots (one plot per target diploid chromosome)
(4) source("coancestry.R") # to load functions required to create co-ancestry plots
(5) plot_coanccurves(ind,coancs,2,2) # plot some co-ancestry curves

