\name{AddDNAObservation}
\alias{AddDNAObservation}
\title{
Add a DNA observation for a specified person and a specified allele system. 
}
\description{
Add a DNA observation for a specified person and a specified allele system. 
}
\usage{
AddDNAObservation(indexperson, indexAlleleSystem, 
indexAllele1, indexAllele2 = indexAllele1)
}
\arguments{
  \item{indexperson}{
The index of the person in the Familias list of persons. 
}
  \item{indexAlleleSystem}{
The index of the allele system in the Familias list of allele systems. 
}
  \item{indexAllele1}{
The index of one of the observed alleles in the list of alleles for this system.
}
  \item{indexAllele2}{
The index of the other observed alleles in the list of alleles for this system.
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{RemoveDNAObservation}}
}
\examples{
NewFamilias()
AddPerson(FALSE)
AddAlleleSystem(c(0.1, 0.1, 0.1))
AddDNAObservation(1, 1, 1, 3)
}
