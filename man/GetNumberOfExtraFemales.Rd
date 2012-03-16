\name{GetNumberOfExtraFemales}
\alias{GetNumberOfExtraFemales}
\title{
Get Number of Extra Females in a Famlias Pedigree
}
\description{
Get the number of females in a Familias pedigree that have been 
added to the persons that appear in all the pedigrees. 
}
\usage{
GetNumberOfExtraFemales(pedigree)
}
\arguments{
  \item{pedigree}{
The index of the pedigree in the Familias list of pedigrees. 
}
}
\value{
The number of extra females in the pedigree. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{GetNumberOfExtraMales}}
}
\examples{
NewFamilias()
AddPedigree(3,3)
GetNumberOfExtraFemales(1)
}
