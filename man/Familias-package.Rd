\name{Familias-package}
\alias{Familias-package}
\alias{Familias}
\docType{package}
\title{Familias: Inferring paternity and indentification based on DNA data
}
\description{
Familias has existed for a number of years as a Windows program for calculating
probabilities in connection with using DNA data to infer family relationships. 
The code is now released in an R package. 
}
\details{
\tabular{ll}{
Package: \tab Familias\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-02-17\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package represents a bare-bones implementation of an interface to the 
core Familias functions, which are programmed in C++. The package itself
functions as a kind of database, where information about persons, 
pedigrees, allele systems and observations for persons are entered stepwise. 
In the end, probability calculations are made. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\references{
For more information, see www.familias.name
}
\keyword{ package }
\examples{
# Minimal computations for a paternity case: 
NewFamilias()
# Mother, putative father and child is added: 
AddPerson(FALSE)
AddPerson(TRUE)
AddPerson(TRUE)
# Mother is specified as paret to child: 
AddFixedRelation(1, 3)
# A first pedigree is added: 
AddPedigree()
# In this pedigree, it is specified that the putative father is the true father
AddRelation(2, 3, 1)
# A second pedigree is added, where the putative father is not the true father
AddPedigree()
# A first allele system is added, with three possible alleles
AddAlleleSystem(c(0.2, 0.3, 0.5))
# A second allele system is added, with three possible alleles
AddAlleleSystem(c(0.1, 0.1, 0.1))
# DNA observations are added for the mother, father, and child: 
AddDNAObservation(1, 1, 1)
AddDNAObservation(1, 2, 1, 2)
AddDNAObservation(2, 1, 2)
AddDNAObservation(2, 2, 1, 2)
AddDNAObservation(3, 1, 1, 2)
AddDNAObservation(3, 2, 1, 2)
# The result is calculated
GetProbabilities()
}
