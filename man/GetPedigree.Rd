\name{GetPedigree}
\alias{GetPedigree}
\title{Get the relationships in a Familias pedigree. 
}
\description{
Get a logical matrix describing the parent-child relationships in a 
pedigree, given its index in the Familias list of pedigrees. 
}
\usage{
GetPedigree(index)
}
\arguments{
  \item{index}{
The index in the Familias list of pedigrees. 
}
}
\value{A square logical matrix of size equal to the number of 
persons in the pedigree. The value in row i and column j is TRUE 
if and only if person i is a parent of person j.
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{GetSizeOfPedigree}}
}
\examples{
NewFamilias()
AddPerson(FALSE)
AddPerson(FALSE)
AddFixedRelation(1,2)
AddPedigree(1,1)
GetPedigree(1)
}
