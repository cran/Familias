\name{GeneratePedigrees}
\alias{GeneratePedigrees}
\title{
Generate possible pedigrees based on information in Familias. 
}
\description{
Based on the list of persons stored in Familias, all possible pedigrees are generated. 
Equivalent pedigrees among existing ones are first removed. Then the additional pedigrees are
generated and added to the Familias list of pedigrees. Additional persons without DNA data
may be specified, so that more complex pedigrees may be generated. 
}
\usage{
GeneratePedigrees(nExtraFemales = 0, nExtraMales = 0)
}
\arguments{
  \item{nExtraFemales}{
The number of extra females, beyond those in the Familias list of persons, that the algorithm can use to create pedigrees. 
}
  \item{nExtraMales}{
The number of extra males, beyond those in the Familias list of persons, that the algorithm can use to create pedigrees. 
}
}
\value{
The output is a logical vector, with the same length as the old number of pedigrees. It indicates which 
pedigrees have been removed from the old list, because they were equivalent to pedigrees earlier in the list. 
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\note{
Unless the number of persons involved is quite limited, this function can take a LONG time, and may seem to not terminate.
In particular, the number of extra males or females specified should be very low, if different from zero.  
}
\seealso{
\code{\link{GetNumberOfPedigrees}}, \code{\link{GetPedigree}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(TRUE)
GeneratePedigrees()
GetNumberOfPedigrees()
}
