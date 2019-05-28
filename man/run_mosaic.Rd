\name{run_mosaic}
\alias{run_mosaic}
\title{run MOSAIC analysis; reads in data and performs full inferential algorithm to estimate parameters and local ancestry. 
}
\description{calls all other MOSAIC functions to read in data, compress, perform inference, and output results. 
}
\usage{run_mosaic(target, folder, chromosomes, A, NUMI, 
    pops = NULL, mask = NULL, PLOT = FALSE, doFst = TRUE, 
    PHASE = TRUE, gens = 0, ratios = NULL, EM = TRUE, 
    ffpath = "/dev/shm/", MC = 0, return.res = TRUE, 
    REPS = 0, GpcM = 60, nl = 1000, max.donors = 100, prop.don = 0.99, 
    doMu = TRUE, doPI = TRUE, dorho = TRUE, dotheta = TRUE, 
    firstind = 1, verbose = TRUE, Ne = 90000) 
}
\arguments{
  \item{target}{a string giving the name of admixed target population
}
  \item{folder}{a string giving the path to folder where inputs are stored
}
  \item{chromosomes}{vector of integers giving the chromosomes to use
}
  \item{A}{integer greater than or equal to 2 giving the number of ancestries to model
}
  \item{NUMI}{integer giving the number of admixed samples to use
}
  \item{pops}{vector of strings giving names of populations; first A are used to simulate admixture if target="simulated" is used. Defaults
    to NULL so that all available populations are used as donor panels
}
  \item{mask}{
}
  \item{PLOT}{
}
  \item{doFst}{
}
  \item{PHASE}{logical specifying whether to rephase targets or not
}
  \item{gens}{
}
  \item{ratios}{
}
  \item{EM}{
}
  \item{ffpath}{
}
  \item{MC}{
}
  \item{return.res}{
}
%  \item{REPS}{integer number of rounds of thin + phase + EM
%}
%  \item{GpcM}{integer granularity of genetic distance grid as number of gridpoints per centiMorgan
%}
%  \item{nl}{
%}
%  \item{max.donors}{
%}
%  \item{prop.don}{
%}
%  \item{doMu}{
%}
%  \item{doPI}{
%}
%  \item{dorho}{
%}
%  \item{dotheta}{
%}
%  \item{firstind}{
%}
%  \item{verbose}{
%}
%  \item{Ne}{
%}
  \item{...}{additional optional arguments
}
}
\value{a list containing inferred parameters, local ancestry estimates, etc. 
}
\author{Michael Salter-Townshend
}
%\seealso{}
%\keyword{}

\examples{
# infer 2-way local ancestry in 2 Moroccans using example data on chromosomes 18,19,20,21,22 
mosaic.result=run_mosaic("Moroccan","example_data/",18:22,2,2)
# simulate and then infer 2-way local ancestry between English and Mandenkans in 2 individuals using example data on chromosomes 18,19,20,21,22 
mosaic.result=run_mosaic("simulated","example_data/",18:22,2,2,c("English","Mandenka"))
}
