\name{AddPerson}
\alias{AddPerson}
\title{
Add a Person to Familias
}
\description{
Adds a person to the list of persons registered in Familias. 
The index of the person in Familias' list of persons is returned. 
}
\usage{
AddPerson(male, isChild = FALSE, yearOfBirth = NULL)
}
\arguments{
  \item{male}{
Should be either TRUE or FALSE. 
}
  \item{isChild}{
Should be either TRUE or FALSE. When a person is designated as a "child", he or she may not have children in pedigrees. 
}
  \item{yearOfBirth}{
Optionally, the year of birth of a person may be given. When pedigrees are 
generated or constructed, it is checked that they are compatible with year of 
birth information. 
}
}
\value{
The index in Familias' list of persons is returned.
}
\author{
Petter Mostad <mostad@chalmers.se> 
}
\seealso{
\code{\link{RemovePerson}}
}
\examples{
AddPerson(FALSE) #Adding an adult female
AddPerson(TRUE, TRUE, 1964) #Adding a male child born 1964
}
