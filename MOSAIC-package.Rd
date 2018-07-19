\name{MOSAIC-package}
\alias{MOSAIC-package}
\alias{MOSAIC}
\docType{package}
\title{
  MOSAIC Organises Segments of Ancestry In Chromosomes.
}
\description{
MOSAIC is a tool for modelling multi-way admixture events using dense haplotype data. The Hidden Markov Model approach is scalable to thousands of samples and, unlike existing methods, does not require specification of surrogate reference populations for the mixing ancestral groups. Instead, potentially complex relationships between multiple observed reference groups and unseen ancestral groups are inferred wherein the reference groups may also be admixed with respect to the ancestral groups. It automatically infers mixing proportions, timings, and parameters related to recombination rate and mutation. Our method outputs estimates of the model parameters along with measurements of model fit and confidence based on the local ancestry estimates and distribution of ancestral chromosome segment lengths with uncertainty. Phasing errors that can potentially mask ancestry are detected and corrected within the algorithm. 
}
\details{
  MOSAIC is designed to run on the linux command line. 
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

There should be a folder with 4 types of input file:
1. phased haplotypes: "pop.genofile.chr" in the format #snps rows and #haps columns.
2. pop names: "sample.names" format unimportant apart from first column should have all the population names.
3. snp files: "snpfile.chr" #snps rowns and 6 columns of rsID, chr, position, distance, allele ?, allele ?. 
4. recombination map: "rates.chr" 3 rows of #sites, position, recombination rate. 

}
\author{
Michael Salter-Townshend, michael.salter-townshend@ucd.ie

Maintainer: michael.salter-townshend@ucd.ie <michael.salter-townshend@ucd.ie>
}
\references{
  "Fine-scale Inference of Ancestry Segments without Prior Knowledge of Admixing Groups" by Salter-Townshend, M. and Myers, S.R. 2018.
}
\keyword{admixture, local ancestry, population genetics, Hidden Markov Model}
\seealso{
  
}
\examples{
  \dontrun{
    mosaic.result=run_mosaic("Moroccan","example_data/",18:22,2,2) # real data example
    mosaic.result=run_mosaic("simulated","example_data/",18:22,2,2,c("English","Mandenka")) # simulated example
  }
}
