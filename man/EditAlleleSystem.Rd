\name{EditAlleleSystem}
\alias{EditAlleleSystem}
\title{
Edit an Allele system (locus) in the Familias list of loci. 
}
\description{Edit an allele system (locus) in the Familias list of loci. 
Various properties of the system, such as number of alleles, frequencies, 
mutation rates, etc, can be changed. A vector can be given specifying the 
correspondence between the old and new alleles, so that the DNA data of 
persons can be updated correctly. 
}
\usage{
EditAlleleSystem(indexSystem, frequencies, 
correspondence = 1:length(frequencies), 
mutationRateFemale = 0, mutationRateMale = 0, 
mutationModelFemale = "stable", mutationModelMale = "stable", 
mutationRangeFemale = 0.1, mutationRangeMale = 0.1, 
silentFrequency = NULL)
}
\arguments{
  \item{indexSystem}{
The index of the system in the Familias list of allele systems. 
}
  \item{frequencies}{A vector specifying the updated allele frequencies. 
}
  \item{correspondence}{
A vector, of the same length as the old number of alleles in the system, 
specifying for each allele the corresponding index in the new system. This
makes sure DNA data of persons can be updated correctly. 
}
  \item{mutationRateFemale}{
The female mutation rate (see \code{\link{AddAlleleSystem}}). 
}
  \item{mutationRateMale}{
The male mutation rate (see \code{\link{AddAlleleSystem}}). 
}
  \item{mutationModelFemale}{
The female mutation "model" (see \code{\link{AddAlleleSystem}}). 
}
  \item{mutationModelMale}{
The male mutation "model" (see \code{\link{AddAlleleSystem}}). 
}
  \item{mutationRangeFemale}{
The female mutation "range" (see \code{\link{AddAlleleSystem}}). 
}
  \item{mutationRangeMale}{
The male mutation "range" (see \code{\link{AddAlleleSystem}}). 
}
  \item{silentFrequency}{
The silent allele frequency (see \code{\link{AddAlleleSystem}}). 
}
}
\value{
None
}
\author{
Petter Mostad <mostad@chalmers.se>
}
\seealso{
\code{\link{AddAlleleSystem}}, \code{\link{RemoveAlleleSystem}}
}
\examples{
NewFamilias()
AddPerson(TRUE)
AddPedigree()
AddAlleleSystem(c(0.2, 0.3, 0.5), mutationRateMale = 0.01)
EditAlleleSystem(1, c(0.3, 0.3, 0.4), mutationRateMale = 0.01)
#Non-default parameters must be re-specified when editing! 
}
