\name{Familias-package}
\alias{Familias-package}
\alias{Familias}
\docType{package}
\title{Familias: Inferring paternity and indentification based on DNA data
}
\description{
Familias has existed for a number of years as a Windows program for calculating
probabilities in connection with using DNA data to infer family relationships. 
The code is now released in an R package. 
}
\details{
\tabular{ll}{
Package: \tab Familias\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-07-17\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package represents an implementation of an interface to the 
core Familias functions, which are programmed in C++. 
Note that this version, 2.0, implements 
a completely different interface to the core functions compared to the previous versions. 
}
\author{
Petter Mostad <mostad@chalmers.se> and Thore Egeland <Thore.Egeland@gmail.com>
}
\references{
For more information, see www.familias.name
}
\keyword{ package }
\examples{
persons <- c("mother", "daughter", "AF")
ped1 <- FamiliasPedigree(id=persons, dadid=c(NA, "AF", NA), 
momid=c(NA, "mother", NA), 
sex=c("female", "female", "male"))
ped2 <- FamiliasPedigree(id=c(persons, "TF"), dadid=c(NA, "TF", NA, NA), 
momid=c(NA, "mother", NA, NA), 
sex=c("female", "female", "male", "male"))
ped3 <- FamiliasPedigree(id=c(persons, "TF", "gf", "gm"), 
dadid = c(NA, "TF", "gf", "gf", NA, NA), 
momid=c(NA, "mother", "gm", "gm", NA, NA), 
sex=c("female", "female", "male", "male", "male", "female"))
par(mfrow=c(3,1))
plot(ped1);title("ped1, i.e., AF is father")
plot(ped2);title("ped2, i.e., AF is unrelated")
plot(ped3);title("ped3, i.e., AF is uncle")
par(mfrow=c(1,1))
mypedigrees <- list(isFather = ped1, unrelated=ped2, isUncle = ped3)
locus1 <- FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("A", "B", "C", "D"))
locus2 <- FamiliasLocus(c(0.2, 0.3, 0.5), "loc2", c(17, 18, 19), 
femaleMutationRate = 0.05)
myloci <- list(locus1, locus2)
datamatrix <- data.frame(locus1.1=c("A", "A", "A"), locus1.2=c("B", "B", "C"), 
locus2.1=c(17, 19, 19), locus2.2=c(18, 18, 18))
rownames(datamatrix) <- persons
FamiliasPosterior(mypedigrees, myloci, datamatrix)
}
