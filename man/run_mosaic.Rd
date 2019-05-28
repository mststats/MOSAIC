\name{run_mosaic}
\alias{run_mosaic}
\title{run MOSAIC analysis; reads in data and performs full inferential algorithm to estimate parameters and local ancestry. 
}
\description{calls all other MOSAIC functions to read in data, compress, perform inference, output results, and create plots. 
}
\usage{run_mosaic(target,folder,chromosomes,a,n)
}
\arguments{
  \item{target}{name of admixed target population
}
  \item{folder}{path to folder where inputs are stored
}
  \item{chromosomes}{vector of chromosomes to use
}
  \item{a}{number of ancestries to model
}
  \item{n}{number of admixed samples to use
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
