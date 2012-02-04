\name{RemoveRelation}
\alias{RemoveRelation}
\title{
Remove a parent-child relation in a Familias pedigree. 
}
\description{
Remove a parent-child relation in a Familias pedigree. 
}
\usage{
RemoveRelation(parentindex, childindex, pedigree)
}
\arguments{
  \item{parentindex}{
The index of the parent in the list of persons of this pedigree. 
}
  \item{childindex}{
The index of the child in the list of persons of this pedigree. 
}
  \item{pedigree}{
The index of the pedigree in the Familias list of pedigrees. 
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se>
}
\seealso{
\code{\link{AddRelation}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(FALSE)
AddPedigree(1,1)
AddRelation(1,3,1)
RemoveRelation(1,3,1)
}
