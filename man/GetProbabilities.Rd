\name{GetProbabilities}
\alias{GetProbabilities}
\title{
Compute prior, likelihood, and posterior probabilities for Pedigrees.  
}
\description{
Compute prior, likelihood, and posterior probabilities for the pedigrees 
stored in Familias given the DNA data stored. Input parameters can be 
used to regulate the prior probability for the pedigrees, and a kinship 
parameter regulates kinship. Note that, before computations, pedigrees that 
are equivalent to earlier pedigrees in the Familias list of pedigrees are removed. 
A list indicating which pedigrees have been kept is returned by the function. 
}
\usage{
GetProbabilities(kinship = 0, generationsParameter = 1, 
inbreedingParameter = 1, partnerParameter = 1, 
onlyprior = FALSE, maxGenerations = NULL)
}
\arguments{
  \item{kinship}{
A non-negative number (usually close to zero) specifying the kinship parameters
used when calculating likelihoods. 
}
  \item{generationsParameter}{
A positive parameter influencing the prior for the pedigrees. The default
1 means that the number of generations in the pedigree does not influence the 
prior.
}
  \item{inbreedingParameter}{
A positive parameter influencing the prior for the pedigrees. The default
1 means that the amount of inbreeding in the pedigree does not influence the 
prior.
}
  \item{partnerParameter}{
A positive parameter influencing the prior for the pedigrees. The default
1 means that the number of partners each parent has 
in the pedigree does not influence the prior.
}
  \item{onlyprior}{If TRUE, only the prior probabilities for the pedigrees
will be computed. 
}
  \item{maxGenerations}{
If given, pedigrees containing more than \code{maxGenerations} generations 
will get prior probability zero. 
}
}
\value{
  \item{posterior}{A vector with lengt equal to the number of pedigrees, 
   specifying the posterior probability for each pedigree.}
  \item{prior}{A vector with length equal to the number of pedigrees, 
   specifying the prior probability for each pedigree.}
  \item{likelihoods}{A vector with length equal to the number of pedigrees, 
   specifying the likelihood of each pedigree.}
  \item{likelihoodsPerSystem}{A matrix with rows corresponding to the 
   different allele systems, and columns corresponding to the pedigrees, 
   specifying the likelihood for each pedigree at each allele system.}
  \item{pedigreeKept}{Before the probability calculations, pedigrees that are 
      equivalent to earlier pedigrees in the Familias list of pedigrees are removed. 
      Thus the results above are for this pruned list of pedigrees. \code{pedigreeKept} is 
      a logical vector indicating which pedigrees from the old list have been kept.}
}
\references{
For more details about the parameters and the computations, 
see www.familias.name. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\examples{
NewFamilias()
AddPerson(FALSE)
AddPerson(TRUE)
AddPerson(TRUE)
AddFixedRelation(1, 3)
AddPedigree()
AddRelation(2, 3, 1)
AddPedigree()
AddAlleleSystem(c(0.2, 0.3, 0.5))
AddAlleleSystem(c(0.1, 0.1, 0.1))
AddDNAObservation(1, 1, 1)
AddDNAObservation(1, 2, 1, 2)
AddDNAObservation(2, 1, 2)
AddDNAObservation(2, 2, 1, 2)
AddDNAObservation(3, 1, 1, 2)
AddDNAObservation(3, 2, 1, 2)
GetProbabilities()
}
