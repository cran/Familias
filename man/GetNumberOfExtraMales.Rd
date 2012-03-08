\name{GetNumberOfExtraMales}
\alias{GetNumberOfExtraMales}
\title{
Get Number of Extra Males in a Famlias Pedigree
}
\description{
Get the number of males in a Familias pedigree that have been 
added to the persons that appear in all the pedigrees. 
}
\usage{
GetNumberOfExtraMales(pedigree)
}
\arguments{
  \item{pedigree}{
The index of the pedigree in the Familias list of pedigrees. 
}
}
\value{
The number of extra males in the pedigree. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{GetNumberOfExtraFemales}}
}
\examples{
NewFamilias()
AddPedigree(3,3)
GetNumberOfExtraMales(1)
}
