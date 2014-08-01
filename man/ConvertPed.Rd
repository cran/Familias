\name{ConvertPed}
\alias{ConvertPed}
\title{
Converts pedigree data frame to kinship2 pedigree object
}
\description{
A data frame describing a pedigree is converted to a
kinship2 pedigree object to be used by say \code{\link{FamiliasPosterior}}.
Typically the input is created using paramlink.
}
\usage{
ConvertPed(ped, persons = NULL)
}

\arguments{
  \item{ped}{
A data frame with columns ID, FID (father id), MID (mother id),sex (1=male, 2=female), AFF (not used)
possible followed by marker data, two columns for each marker
}
  \item{persons}{
A character vector indicating names of persons in the case ordered as \code{ped[,1].
In other words,
persons[1] corresponds tp ped[1,1] etc.}
}
}
\value{
\item{ped }{An object of class pedigree}
\item{datamatrix }{A data frame. One line for each person, one column for each allele.}
}

\author{
Thore Egeland <Thore.Egeland@gmail.com>
}
\seealso{
\code{\link{FamiliasPosterior}}
}

\examples{
persons <- c("AF","mother", "daughter")
ped1=swapSex(nuclearPed(1),3)
locus1=marker(ped1,1,c(1,1),2,c(2,2),3,c(1,2))
ped1=addMarker(ped1,locus1)
# A plot function such as the following can now be used: 
# plot(ped1,marker=1,id.labels=persons,title="ped1") 
ped1=as.data.frame(ped1) 
ConvertPed(ped1,persons)
}
