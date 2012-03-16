\name{AddPedigree}
\alias{AddPedigree}
\title{Add a Pedigree to the List of Pedigrees
}
\description{
Add a pedigree to the list of pedigrees stored in Familias. 
Initially, the pedigree contains only the relations fixed with 
\code{\link{AddFixedRelation}}. A number of extra females and males
may be specified. Additional relationships may then be specified 
using \code{\link{AddRelation}}. 
}
\usage{
AddPedigree(nExtraFemales = 0, nExtraMales = 0)
}
\arguments{
  \item{nExtraFemales}{
The number of extra females (who may later be related to others in the 
pedigree using \code{\link{AddRelation}}). The extra females are given 
indexes after the persons that appear in all pedigrees, but before the 
extra males. 
}
  \item{nExtraMales}{
The number of extra males (who may later be related to others in the 
pedigree using \code{\link{AddRelation}}). The extra males are given 
indexes after all other persons in the pedigree. 
}
}
\value{
The index of the new pedigree in the list of pedigrees is returned. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{RemovePedigree}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(FALSE)
AddPerson(FALSE)
AddPedigree()
}
