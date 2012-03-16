\name{AddRelation}
\alias{AddRelation}
\title{
Add a parent-child relation to a Familias pedigree. 
}
\description{
Add a parent-child relation to a Familias pedigree. 
}
\usage{
AddRelation(parentindex, childindex, pedigree)
}
\arguments{
  \item{parentindex}{
The index of the parent in the list of persons for this pedigree. 
}
  \item{childindex}{
The index of the child in the list of persons for this pedigree. 
}
  \item{pedigree}{
The index of this pedigree in the Familias list of pedigrees. 
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{RemoveRelation}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(TRUE)
AddPedigree()
AddRelation(1, 2, 1)
}
