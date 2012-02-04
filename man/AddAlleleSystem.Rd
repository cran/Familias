\name{AddAlleleSystem}
\alias{AddAlleleSystem}
\title{
Add an Allele system (locus) to the Familias list of loci. 
}
\description{Add an allele system (locus) to the Familias list of loci. 
The various properties of the system is specified, such as allele frequencies, 
mutation rates, etc. 
}
\usage{
AddAlleleSystem(frequencies, 
mutationRateFemale = 0, mutationRateMale = 0, 
mutationModelFemale = "stable", mutationModelMale = "stable", 
mutationRangeFemale = 0.1, mutationRangeMale = 0.1, 
silentFrequency = NULL)
}
\arguments{
  \item{frequencies}{A vector containing the allele frequencies of the system. Note that the 
specified frequencies may add up to less than 1: It is then assumed that the system contains 
other alleles one has chosen not to specify. 
}
  \item{mutationRateFemale}{The female mutation rate. 
}
  \item{mutationRateMale}{The male mutation rate. 
}
  \item{mutationModelFemale}{The female mutation "model" (see below).  
}
  \item{mutationModelMale}{The male mutation "model" (see below). 
}
  \item{mutationRangeFemale}{The female mutation "range" (see below). 
}
  \item{mutationRangeMale}{The male mutation "range" (see below). 
}
  \item{silentFrequency}{
One may specify a frequency for silent (i.e., unobserved in data) alleles. 
This will then influence the probability of observing homozygote data. (It is legal, though 
hardly meaningful, to specify the silent allele in data: It is then specified with the index
after the other alleles specified in \code{frequencies}). 
}
}
\details{The mutation model parameters regulate the probabilities of which allele a mutation 
ends up with, given that a mutation occurs. They must be either "equal" (meaning equal probabilities
for all alleles), "frequencies" (meaning that the probability is proportional to each allele's frequency), 
"step" (meaning that the probability is highest for alleles close to the unmutated allele in the list 
of alleles), or "stable" (same as "step", except that the probabilities have been adapted so that repeated 
generations of mutations will not change the allele frequencies from those given in \code{frequencies}). 
The mutation range parameters are meaningful only for the "step" and "stable" models. They then regulate
how fast the probability falls with the distance from the unmutated allele. 
}
\value{
The index of the allele system in the Familias list of allele systems. 
}
\references{
For more details, please go to www.familias.name
}
\author{
Petter Mostat <mostad@chalmers.se> 
}
\seealso{
\code{\link{EditAlleleSystem}}, \code{\link{RemoveAlleleSystem}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPedigree()
AddAlleleSystem(c(0.2, 0.3, 0.5))
AddDNAObservation(1, 1, 1, 2)
GetProbabilities()
}
