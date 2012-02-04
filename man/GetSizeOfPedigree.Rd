\name{GetSizeOfPedigree}
\alias{GetSizeOfPedigree}
\title{
Get the number of persons in a Familias pedigree. 
}
\description{
Get the number of persons in a Familias pedigree, given its index. 
}
\usage{
GetSizeOfPedigree(index)
}
\arguments{
  \item{index}{
The index of the pedigree in the Familias list of pedigrees. 
}
}
\value{
The number of persons in the pedigree. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{GetPedigree}}
}
\examples{
NewFamilias()
AddPerson(FALSE)
AddPerson(FALSE)
AddPedigree(2,2)
GetSizeOfPedigree(1)
}
