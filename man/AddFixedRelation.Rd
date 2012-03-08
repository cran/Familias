\name{AddFixedRelation}
\alias{AddFixedRelation}
\title{Add a Fixed Relation between Familias persons. 
}
\description{
Adds a parent-child relationship between two persons in Familias' list of 
persons. The relationship is added to all registered pedigrees. For some
of the pedigrees, the relationship may be incompatible with already registered
relationships. These pedigrees are then removed from the list of pedigrees,
and a vector is returned with TRUE for the pedigrees that are kept and FALSE 
for those that are removed. 
}
\usage{
AddFixedRelation(parentindex, childindex)
}
\arguments{
  \item{parentindex}{
The index of the parent in Familias' list of persons. 
}
  \item{childindex}{
The index of the child in Familias' list of persons. 
}
}
\value{
A logical vector is returned of the same length as the (old) number of 
pedigrees. A side effect is that pedigrees with relationships incompatible 
with the given new relationship are removed from the list of pedigrees. 
The vector is TRUE for those pedigrees that are kept, and FALSE for those 
that are removed. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{RemoveFixedRelation}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(TRUE)
AddFixedRelation(1, 2) #Adding a father-son relationship
}
