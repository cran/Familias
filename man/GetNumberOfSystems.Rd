\name{GetNumberOfSystems}
\alias{GetNumberOfSystems}
\title{
Get Number of Allele Systems (loci) stored in Familias. 
}
\description{
Get the number of allele systems (loci) that are stored in Familias. 
}
\usage{
GetNumberOfSystems()
}
\value{
The number of allele systems stored in Familias. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{AddAlleleSystem}}
}
\examples{
NewFamilias()
AddAlleleSystem(c(0.1, 0.1, 0.1))
AddAlleleSystem(c(0.2, 0.2, 0.2, 0.2))
GetNumberOfSystems()
}
