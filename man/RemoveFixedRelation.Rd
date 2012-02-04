\name{RemoveFixedRelation}
\alias{RemoveFixedRelation}
\title{
Remove a parent-child Relation
}
\description{
Remove a parent-child relation from the list of relations that all 
pedigrees contain. The indices of the parent and the child in the 
list of persons is the input. 
}
\usage{
RemoveFixedRelation(parentindex, childindex)
}
\arguments{
  \item{parentindex}{
The index of the parent in the list of persons stored in Familias. 
}
  \item{childindex}{
The index of the child in the list of persons stored in Familias. 
}
}
\value{None
}
\author{
Petter mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{AddFixedRelation}}
}
\examples{
NewFamilias()
AddPerson(1)
AddPerson(1)
AddFixedRelation(1,2)
RemoveFixedRelation(1,2)
}
