\name{GetParents}
\alias{GetParents}
\title{
Finds the parents of each person in a Familias pedigree. 
}
\description{For a pedigree in the Familias list of pedigrees, 
two vectors are returned, containing a list of each persons mother 
respective father within the pedigree.  
}
\usage{
GetParents(index)
}
\arguments{
  \item{index}{
The index in the Familias list of pedigrees. 
}
}
\value{
A list with two components, called mother and father. 
Each has the same length as the number of persons in this pedigree, 
and for each lists the index of the persons mother, respective father. 
If the person has no such parent in the pedigree, NA is returned. 
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
AddPerson(TRUE)
AddFixedRelation(1,2)
AddPedigree(1,1)
GetParents(1)
}
