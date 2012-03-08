\name{RemovePerson}
\alias{RemovePerson}
\title{Remove a Person from Familias
}
\description{
Removes a person from Familias' list of persons, based on the index of 
the person in Familias' list. 
}
\usage{
RemovePerson(index)
}
\arguments{
  \item{index}{
The index of the person to be removed in Familias' list of persons. 
}
}
\value{None. 
}
\author{
Petter mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{AddPerson}}
}
\examples{
NewFamilias()
AddPerson(TRUE) # Adds a male. 
AddPerson(FALSE) # Adds a female. 
RemovePerson(1) # Removes the male. The female gets index 1. 
}
