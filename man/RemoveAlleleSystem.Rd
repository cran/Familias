\name{RemoveAlleleSystem}
\alias{RemoveAlleleSystem}
\title{
Remove Allele System (locus) from the Familias list of loci. 
}
\description{
Remove an allele system (locus) from the Familias list of loci, given 
the index of the system in the list. 
}
\usage{
RemoveAlleleSystem(index)
}
\arguments{
  \item{index}{
The index of the system in the Familias list of allele systems. 
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{AddAlleleSystem}}, \code{\link{EditAlleleSystem}}
}
\examples{
NewFamilias()
AddPerson(FALSE)
AddPedigree()
AddAlleleSystem(c(0.1, 0.1, 0.1))
RemoveAlleleSystem(1)
}
