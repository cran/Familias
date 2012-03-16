\name{RemovePedigree}
\alias{RemovePedigree}
\title{
Remove a pedigree from Familias list of pedigrees. 
}
\description{
Remove a pedigree from the Familias list of pedigrees. The input is 
the index of the pedigree in the list. 
}
\usage{
RemovePedigree(index)
}
\arguments{
  \item{index}{
The index of the pedigre in the Familias list of pedigrees.
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se>
}
\seealso{
\code{\link{AddPedigree}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(TRUE)
AddPedigree()
RemovePedigree(1)
}
