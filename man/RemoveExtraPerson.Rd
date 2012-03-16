\name{RemoveExtraPerson}
\alias{RemoveExtraPerson}
\title{
Remove an extra person from a Familias pedigree. 
}
\description{
Remove an extra person from a Familias pedigree. All pedigrees contain the 
persons generated with \code{\link{AddPerson}}, but each pedigree may 
contain extra persons. Such persons may be removed with this command. 
}
\usage{
RemoveExtraPerson(person, pedigree)
}
\arguments{
  \item{person}{
The index of the person to be removed in the list of persons in the pedigree. 
}
  \item{pedigree}{
The index of the pedigree in the Familias list of pedigrees. 
}
}
\value{
None
}
\author{
Petter Mostad
}
\seealso{
\code{\link{AddExtraPerson}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPerson(FALSE)
AddPerson(FALSE)
AddPedigree(1, 1)
RemoveExtraPerson(4, 1)
}
