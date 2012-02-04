\name{AddExtraPerson}
\alias{AddExtraPerson}
\title{
Add extra person to a Familias pedigree. 
}
\description{
Add an extra person (male or female) to a specific Familias pedigree. 
}
\usage{
AddExtraPerson(male, pedigree)
}
\arguments{
  \item{male}{
Should be TRUE for male or FALSE for female. 
}
  \item{pedigree}{
The index in the Familias list of pedigrees of the pedigree to which 
the person should be added. 
}
}
\value{
None
}
\note{If an extra male is added, his index will be after any other extra males 
for this pedigree, but before extra females. Thus the indices of all extra females
will be increased by 1.
} 
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{RemoveExtraPerson}}, \code{\link{AddRelation}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(FALSE)
AddPedigree(1,1)
AddExtraPerson(TRUE, 1)
#The person in pedigree 1 with index 5 is female. 
}
