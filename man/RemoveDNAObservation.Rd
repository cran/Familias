\name{RemoveDNAObservation}
\alias{RemoveDNAObservation}
\title{
Remove a DNA observation for a specified person and a specified allele system. 
}
\description{
Remove a DNA observation for a specified person and a specified allele system. 
}
\usage{
RemoveDNAObservation(indexperson, indexAlleleSystem)
}
\arguments{
  \item{indexperson}{
The index of the person in the Familias list of persons. 
}
  \item{indexAlleleSystem}{
The index of the allele system in the Familias list of allele systems. 
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{AddDNAObservation}}
}
\examples{
NewFamilias()
AddPerson(FALSE)
AddAlleleSystem(c(0.1, 0.1, 0.1))
AddDNAObservation(1, 1, 1, 3)
RemoveDNAObservation(1, 1)
}
