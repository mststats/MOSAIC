\name{run_mosaic}
\alias{run_mosaic}
\title{run MOSAIC analysis; reads in data and performs full inferential algorithm to estimate parameters and local ancestry. 
}
\description{calls all other MOSAIC functions to read in data, compress, perform inference, and output results. 
}
\usage{run_mosaic(target, datasource, chrnos, A, NUMI, 
pops = NULL, mask = NULL, PLOT = TRUE, doFst = TRUE, 
PHASE = TRUE, gens = NULL, ratios = NULL, EM = TRUE, 
ffpath = tempdir(), MC = 0, return.res = TRUE, 
REPS = 0, GpcM = 60, nl = 1000, max.donors = 100, prop.don = 0.99, 
doMu = TRUE, doPI = TRUE, dorho = TRUE, dotheta = TRUE, 
firstind = 1, verbose = TRUE, Ne = 90000, MODE="DIP", singlePI=FALSE,
init.rho=NULL, init.theta=NULL, init.Mu=NULL, init.PI=NULL) 
}
\arguments{
  \item{target}{a string giving the name of admixed target population.
  }
  \item{datasource}{a string giving the path to folder where inputs are stored.
  }
  \item{chrnos}{vector of integers giving the chromosomes to use.
  }
  \item{A}{integer greater than or equal to 2 giving the number of ancestries to model.
  }
  \item{NUMI}{integer giving the number of admixed samples to use.
  }
  \item{pops}{vector of strings giving names of populations; first A are used to simulate admixture if target="simulated" is used. 
    Defaults to NULL so that all available populations are used as donor populations when target is not simulated. 
  }
  \item{mask}{vector of strings giving names of populations not to use. 
    Defaults to NULL so that no populations are dropped from the analysis. 
  }
  \item{PLOT}{logical specifying whether MOSAIC should save default plots to summarise results in folder MOSAIC_PLOTS/.
  }
  \item{doFst}{logical to set whether or not to calculate Fst between each pair of ancestries (using estimated local ancestry). 
    Defaults to TRUE and run times are sped up by using FALSE. 
  }
  \item{PHASE}{logical specifying whether to rephase targets or not. Defaults to TRUE. 
  }
  \item{gens}{integer number of generations ago from which to simulate admixture (if target is simulated). 
  }
  \item{ratios}{vector of proportions; when simulating data (if target is simulated), what vector of proportions should be used for the ancestries? 
    Default NULL causes MOSAIC to choose a random vector that is close to equal proportions of 1/A. 
  }
  \item{EM}{logical to set whether or not to perform EM inference on model parameters. Defaults to TRUE. 
  }
  \item{ffpath}{string path to where temporary read / write files are stored. Default is "/dev/shm"
  }
  \item{MC}{integer number of CPU cores to hog for the MOSAIC inference. Defaults to use half of all available cores. 
  }
  \item{return.res}{logical to set whether a list containing all inferred parameters, local ancestry estimates, Fst values, is returned. Defaults to TRUE. 
  }
  \item{REPS}{integer number of rounds of thin + phase + EM
  }
  \item{GpcM}{integer granularity of genetic distance grid as number of gridpoints per centiMorgan. Defaults to 60. 
  }
  \item{nl}{integer maximum number of individuals to use from each reference population. Defaults to 1000 and where this is larger 
    than the number of individuals in a population all are used. 
  }
  \item{max.donors}{integer specifying the maximum number of top ranked most useful haplotypes to consider copying from for each target individual
    at each gridpoint. Defaults to 100. The ranked top 100 are found in the thinning part of the algorithm and this leads to a large speedup in
    fitting the two-layer HMM with a very small decrease in accuracy of local ancestry estimation. 
  }
  \item{prop.don}{proportion specifying the maximum cumulative probability of copying from for each target individual
    at each gridpoint. Defaults to 0.99. Top ranked max.donors or until copying cumulative probability from top down reaches prop.don
    are considered at each gridpoint. 
  }
  \item{doMu}{logical specifying whether to update the copying matrix Mu parameters. 
    Defaults to TRUE but can be turned off to force MOSAIC to fit a particular version of a model e.g. using parameters learnt previously. 
  }
  \item{doPI}{logical specifying whether to update the ancestry transition matrix parameters in the HMM. 
    Defaults to TRUE but can be turned off to force MOSAIC to fit a particular version of a model. 
  }
  \item{dorho}{logical specifying whether to update the recombination rate within an ancestry parameter in the HMM. 
    Defaults to TRUE but can be turned off to force MOSAIC to fit a particular version of a model. 
  }
  \item{dotheta}{logical specifying whether to update the miscopying rate (due to mutation + error) parameter (or emission probability in the HMM). 
    Defaults to TRUE but can be turned off to force MOSAIC to fit a particular version of a model. 
  }
  \item{firstind}{integer index of first target individual to analyse from the input dataset. Use with NUMI to analyses batches of individuals
    or do them one at a time. Defaults to 1. 
  }
  \item{verbose}{logical specifying whether to print algorithm updates to screen during the inferential routine. 
    Defaults to TRUE. 
  } 
  \item{Ne}{integer specifying effective population size for the species. Only usage is in estimating initial value of
    theta above which will be re-estimated by default as part of the algorithm. 
  }
  \item{MODE}{either "DIP" or "HAP" to run coancestry curve estimation on diploid or haploid genomes. Note that if "HAP" then PHASE should be FALSE.
    Defaults to "DIP".
  }
  \item{singlePI}{should MOSAIC be constrained to assume that all admixed genomes experienced the same admixture event (timings and proportions)?
    Defaults to FALSE. 
  }
  \item{init.rho}{User specified values for recombination rate within latent ancestry. Could be output from previous MOSAIC run. 
  }
  \item{init.theta}{User specified values for error / mutation rate within latent ancestry. Could be output from previous MOSAIC run. 
  }
  \item{init.Mu}{User specified values for copying rates from panels within each ancestry. Could be output from previous MOSAIC run. 
  }
  \item{init.PI}{User specified values for switching rates between ancestries. Could be output from previous MOSAIC run. Note this overrides ratios and gens values. 
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
\dontrun{mosaic.result=run_mosaic("Moroccan","example_data/",18:22,A=2,NUMI=4)}
# simulate and then infer 2-way local ancestry between English and Mandenkans in 2 individuals using example data on chromosomes 18,19,20,21,22 
\dontrun{mosaic.result=run_mosaic("simulated","example_data/",18:22,A=2,NUMI=3,c("English","Mandenka"))}
}
