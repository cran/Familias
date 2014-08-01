\name{FamiliasPosterior}
\alias{FamiliasPosterior}
\title{
Calculates posterior probabilities and likelihoods for pedigrees
}
\description{
The calculations of the windows version of Familias 2.0
are made available in an R environment.
}
\usage{
FamiliasPosterior(pedigrees, loci, datamatrix, prior, ref = 1, 
kinship = 0, simplifyMutations = FALSE)
}

\arguments{
  \item{pedigrees}{
An object of type 'pedigree' or 'FamiliasPedigree', or a list of such objects
}
  \item{loci}{
A FamiliasLocus object or a list of such objects
}
  \item{datamatrix}{
A data frame. The row names must be the names of the persons you have data for. The columns contain
the alleles, two columns for each marker, in the same order used in the loci list. 
}
  \item{prior}{
A vector of non-negative probabilities summing to 1. As default a flat prior is used. 
}
  \item{ref}{
The index in the list of pedigrees of the pedigree which should be used as reference when computing likelihood ratios. 
The default value is 1. 
}
  \item{kinship}{
A real in [0,1], commonly denoted theta in forensics included to model sub-population effects as departures
from Hardy-Weinberg equilibrium. The default value is zero.
}
  \item{simplifyMutations}{
In pedigrees with several generations multistep mutations may happen. If the probability of mutating to an allele 
depends on which allele it mutates from, exact likelihood computations must keep track of all the possible values of 
mutated alleles in such multistage mutations, and this may slow down computations considerably. Instead, one may 
use in computations for multistage mutations a single allele that is not among those observed in the data. 
When this approach gives exact results it is always used; 
in other cases one may choose to use it as an approximation by setting 
simplifyMutations to TRUE. The properties of the extra allele is the weighted average (by population frequencies) 
of the alleles not observed in the data. 
}
}
\value{
\item{posterior}{probabilities for each pedigree}
\item{prior}{prior returned for convenience and backward compatibility}
\item{LR}{Likelihood ratios using the pedigree indicated with the ref parameter (default value 1) as basis}
\item{LRperMarker}{Likelihood ratios per marker using the pedigree indicated with the ref parameter (default value 1) as basis}
\item{likelihoods}{for each pedigree}
\item{likelihoodsPerSystem}{likelihoods for each locus and each pedigree}
}

\author{
Petter Mostad <mostad@chalmers.se> and Thore Egeland <Thore.Egeland@gmail.com>
}

\examples{
#Example 1
#Data is available for "mother", "daughter" and "AF" for two loci (systems). 
#Three pedigrees are defined, with "mother" being the mother of "daughter" 
#in all cases. "AF" may be the father (ped1), unrelated (ped1) or
#uncle (ped3). The posterior probabilities for the pedigrees are calculated 
#and likelihoods are also given so that
#likelihood ratios can be computed. Compared to the windows version of Familias 2.0, 
#it is easy to plot pedigrees and define arbitrary priors for the three alternative pedigrees.
#The below implementation uses the R package kinship2 to define and plot pedigrees

persons <- c("mother", "daughter", "AF")
ped1 <- pedigree(id=persons, dadid=c(NA, "AF", NA), momid=c(NA, "mother", NA), 
sex=c("female", "female", "male"))
ped2 <- pedigree(id=c(persons, "TF"), dadid=c(NA, "TF", NA, NA), momid=c(NA, "mother", NA, NA), 
sex=c("female", "female", "male", "male"))
ped3 <- pedigree(id=c(persons, "TF", "gf", "gm"), dadid = c(NA, "TF", "gf", "gf", NA, NA), 
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
locus2 <- FamiliasLocus(c(0.2, 0.3, 0.5), "loc2", c(17, 18, 19), femaleMutationRate = 0.05)
myloci <- list(locus1, locus2)
datamatrix <- data.frame(locus1.1=c("A", "A", "A"), locus1.2=c("B", "B", "C"), 
locus2.1=c(17, 19, 19), locus2.2=c(18, 18, 18))
rownames(datamatrix) <- persons
result=FamiliasPosterior(mypedigrees, myloci, datamatrix,ref=2)

#Example 2. Example 1 continued.
#Below commands (from scratch) are based on paramlink to prepare pedigrees
persons <- c("AF","mother", "daughter")
ped1=swapSex(nuclearPed(1),3)
ped1=as.data.frame(ped1)
ped1.c=ConvertPed(ped1,persons)$ped
plot(ped1.c,main="ped1")

ped1=swapSex(nuclearPed(1),3)
AF=singleton(4)
ped2=rbind(as.data.frame(ped1),as.data.frame(AF))
persons <- c("TF","mother", "daughter","AF")
ped2.c=ConvertPed(ped2,persons)$ped
plot(ped2.c,main="ped2")

ped=nuclearPed(2)
ped=addOffspring(ped,father=4,noff=1,sex=2)
persons=c("gf","gm","AF","TF","mother","daughter")
plot(ped,id.labels=persons,main="ped3")
ped=as.data.frame(ped)
ped3.c=ConvertPed(ped,persons)$ped
mypedigrees <- list(isFather = ped1.c, unrelated=ped2.c, isUncle = ped3.c)

locus1 <- FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("A", "B", "C", "D"))
locus2 <- FamiliasLocus(c(0.2, 0.3, 0.5), "loc2", c(17, 18, 19), femaleMutationRate = 0.05)
myloci <- list(locus1, locus2)
datamatrix <- data.frame(locus1.1=c("A", "A", "A"), locus1.2=c("B", "B", "C"), 
locus2.1=c(17, 19, 19), locus2.2=c(18, 18, 18))
rownames(datamatrix) <- c("mother", "daughter", "AF")
result=FamiliasPosterior(mypedigrees, myloci, datamatrix, ref=2)

#Example 3. Example 2 continued
#Below commands (from scratch) are based on paramlink to prepare pedigrees and datamatrix.
#The handling of datamatrix is artificial in *this* example; alternatively genotypes can be provided
#in other, perhaps more natural ways.
persons <- c("AF","mother", "daughter") #Names in paramlink are integers 1:length(persons)={1,2,3}
ped1=swapSex(nuclearPed(1),3)
dat <- data.frame(locus1.1=c("A", "A", "A"), locus1.2=c("C","B", "B"), 
locus2.1=c(19,17, 19), locus2.2=c(18, 18, 18)) #Consistent with persons above
ped1=cbind(as.data.frame(ped1),dat)
foo=ConvertPed(ped1,persons)
ped1.c=foo$ped
datamatrix=foo$datamatrix

ped1=swapSex(nuclearPed(1),3)
AF=singleton(4)
ped2=rbind(as.data.frame(ped1),as.data.frame(AF))
persons <- c("TF","mother", "daughter","AF")
ped2.c=ConvertPed(ped2,persons)$ped

ped=nuclearPed(2)
ped=addOffspring(ped,father=4,noff=1,sex=2)
persons=c("gf","gm","AF","TF","mother","daughter")
ped=as.data.frame(ped)
ped3.c=ConvertPed(ped,persons)$ped
mypedigrees <- list(isFather = ped1.c, unrelated=ped2.c, isUncle = ped3.c)

locus1 <- FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("A", "B", "C", "D"))
locus2 <- FamiliasLocus(c(0.2, 0.3, 0.5), "loc2", c(17, 18, 19), femaleMutationRate = 0.05)
myloci <- list(locus1, locus2)

result=FamiliasPosterior(mypedigrees, myloci, datamatrix, ref=2)

#Example 4: Using FamiliasPedigree
# When using non-stationary mutation models, results change when persons with no data and 
# only one child are added to pedigrees. To avoid adding such "spurious" persions, 
# FamiliasPedigree should be used instead of the pedigree function. An example where 
# data is available only for a person and his alledged father: 
persons <- c("person", "AF")
sex <- c("male", "male")
ped1 <- FamiliasPedigree(id=persons, dadid=c(NA, NA), momid=c(NA, NA), sex=sex)
ped2 <- FamiliasPedigree(id=persons, dadid=c("AF", NA), momid=c(NA, NA), sex=sex)
mypedigrees <- list(unrelated = ped1, isFather=ped2)
locus1 <- FamiliasLocus(c(0.1, 0.2, 0.3, 0.4), "locus1", c("A", "B", "C", "D"),
   maleMutationModel = "equal", maleMutationRate = 0.005)
locus2 <- FamiliasLocus(c(0.2, 0.3, 0.5), "locus2", c(17, 18, 19),
   maleMutationModel = "equal", maleMutationRate = 0.005) 
myloci <- list(locus1, locus2)
datamatrix <- data.frame(locus1.1=c("A", "A"), locus1.2=c("B", "B"), 
locus2.1=c(17, 19), locus2.2=c(18, 18))
rownames(datamatrix) <- persons
result <- FamiliasPosterior(mypedigrees, myloci, datamatrix)

#Example 5: User-specified mutation matrices
persons <- c("son", "mother", "AF")
sex <- c("male", "female", "male")
ped1 <- FamiliasPedigree(id=persons, dadid=c(NA, NA, NA), momid=c("mother", NA, NA), sex=sex)
ped2 <- FamiliasPedigree(id=persons, dadid=c("AF", NA, NA), momid=c("mother", NA, NA), sex=sex)
mypedigrees <- list(unrelated = ped1, isFather=ped2)
locus1 <- FamiliasLocus(c(0.1, 0.2, 0.3, 0.4), "locus1", 
   maleMutationModel = "custom", maleMutationMatrix = matrix(
   c(0.99, 0.005, 0.003, 0.002, 
     0.004, 0.99, 0.004, 0.002, 
     0.002, 0.004, 0.99, 0.004, 
     0.002, 0.003, 0.005, 0.99), 
   4, 4, byrow=TRUE), femaleMutationModel = "custom", femaleMutationMatrix = matrix(
   c(0.99, 0.005, 0.003, 0.002, 
     0.004, 0.99, 0.004, 0.002, 
     0.002, 0.004, 0.99, 0.004, 
     0.002, 0.003, 0.005, 0.99), 
   4, 4, byrow=TRUE))
datamatrix <- data.frame(locus1.1=c("A", "A", "C"), locus1.2=c("B", "A", "C"))
rownames(datamatrix) <- persons
result <- FamiliasPosterior(mypedigrees, locus1, datamatrix)


#Example 6: Using kinship: 
persons <- c("AF","mother", "daughter")
ped1=swapSex(nuclearPed(1),3)
locus1=marker(ped1,1,c(1,1),2,c(2,2),3,c(1,2))
ped1=addMarker(ped1,locus1)
plot(ped1,marker=1,id.labels=persons,main="ped1") #plot function of paramlink
ped1=as.data.frame(ped1) 
foo=ConvertPed(ped1,persons)
ped1.c=foo$ped
datamatrix=foo$datamatrix
locus1 <- FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", allelenames= 1:4)
FamiliasPosterior(ped1.c, locus1, datamatrix) #likelihood=0.1^2*0.2^2=4e-04
th=0.03;FamiliasPosterior(ped1.c, locus1, datamatrix,kinship=th) #Formula below
((th+(1-th)*0.1)*(1-th)*0.1*(th*0.2+(1-th)*0.2^2))/((1+2*th)*(1+th)) #OK

# Example 7: Exercise S11 Silent allele
# See http://arken.umb.no/~theg/alcala/familiasExercises.pdf
ped1=nuclearPed(1)
persons <- c("AF","mother", "child")
locus1=marker(ped1,1,c("B","B"),2,c("A","A"),3,c("A","A"),alleles=c("A","B","C"))
ped1=as.data.frame(addMarker(ped1,locus1))
foo=ConvertPed(ped1,persons)
ped1.c=foo$ped
datamatrix=foo$datamatrix

ped2=nuclearPed(1)
ped2=relabel(ped2,old=c(1,2,3),new=c(4,2,3))
AF=singleton(1)
ped2=rbind(as.data.frame(AF),as.data.frame(ped2))
persons <- c("AF","mother", "child","TF")
ped2.c=ConvertPed(ped2,persons)$ped
mypedigrees=list(ped1.c,ped2.c)
locus1 <- FamiliasLocus(frequencies=c(0.1, 0.1,0.75, 0.05), name="locus1", 
allelenames= c("A","B","C", "silent"))
result=FamiliasPosterior(mypedigrees, locus1, datamatrix, ref=2) #LR=1.36 OK
}
