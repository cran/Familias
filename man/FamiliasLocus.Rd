\name{FamiliasLocus}
\alias{FamiliasLocus}

\title{
Creates an object with information on a locus, including its mutation matrices. 
}
\description{
The user provides input needed to define a locus (also called system or forensic marker) to be used for pedigree calculations.
The input is checked and if no errors are found a list with class \code{FamiliasLocus} is returned containing
the information.
}
\usage{
FamiliasLocus(frequencies, name, allelenames, femaleMutationModel = "stable", 
maleMutationModel = "stable", femaleMutationRate = 0, maleMutationRate = 0, 
femaleMutationRange = 0.5, maleMutationRange = 0.5, 
femaleMutationMatrix, maleMutationMatrix)
}

\arguments{
  \item{frequencies}{
A vector containing the allele frequencies of the system, including a possible silent allele. The frequencies must sum to 1. 
}
  \item{name}{
Characters like "D3S1358", used to identify the locus (system). 
}
  \item{allelenames}{
Names of the alleles, like \code{15} or "A". Note that thel last allele may be called 'silent' (or 'Silent'). 
It is then treated as a silent allele in subsequent likelihood calculations. 
}
  \item{femaleMutationModel}{
The female mutation model: It may be 'equal', 'frequencies', 'step', 'stable', or 'custom'.
}
  \item{maleMutationModel}{
The male mutation model: It may be 'equal', 'frequencies', 'step', 'stable', or 'custom'.
}
  \item{femaleMutationRate}{
The female mutation rate (not used when the female mutation model is 'custom'). 
}
  \item{maleMutationRate}{
The male mutation rate (not used when the male mutation model is 'custom'). 
}
  \item{femaleMutationRange}{
The female mutation 'range' (only used when the female mutation model is 'step' or 'stable').
}
  \item{maleMutationRange}{
The male mutation 'range' (only used when the male mutation model is 'step' or 'stable').
}
  \item{femaleMutationMatrix}{
Only used when the female mutation model is 'custom'. 
}
  \item{maleMutationMatrix}{
Only used when the male mutation model is 'custom'. 
}
}
\details{
The probabilities for when and how mutations happen can be specified in mutation matrices, where 
the row corresponding to an allele indicates the probabilities that the allele is transferred as 
the allele indicated by the column. Mutation matrices may be specified directly by setting the 
value of the mutation model parameter to 'custom'. Otherwise, the values of this parameter, 
together with values of the rate and range parameters, specify how the mutation matrix is constructed. 
The value 'equal' means equal probability for mutating to any allele, given that a mutation happens. 
It is referred to as "Equal probability (simple and fast)" in Familias 2.0. The value 'frequencies' 
means that the probability to mutate to an allele is proportional to the allele's frequency. It is 
referred to as "Probability proportional to frequency (stable)" in Familias 2.0. 
The value 'step' means that the probability is highest for alleles close to the unmutated allele 
in the list of alleles. How fast the probability decreases with range is regulated with the range parameter. 
The mutation model is referred to as "Probability decreasing with range (equal)" in Familias 2.0. 
Finally, the value 'stable' is similar to 'step', except that the probabilities have been modified so that repeated
generations of mutations will not change the allele frequencies from those given in frequencies. The mutation model 
is referred to as "Probability decreasing with range (stable)" in Familias 2.0.
}
\value{
A list of class \code{FamiliasLocus} containing 
\item{locusname}{The name of the locus}
\item{alleles}{The frequencies of the alleles. The names of the alleles are included as the vector names.}
\item{femaleMutationType}{A string specifying the type of the female mutations.}
\item{femaleMutationMatrix}{The mutation matrix used for female transfer.}
\item{maleMutationType}{A string specifying the type of the male mutations.}
\item{maleMutationMatrix}{The mutation matrix used for male transfer.}
\item{simpleMutationMatrices}{Indicates whether the probability of mutating to an 
allele is always independent of which allele the mutation happens from. 
If this is true, some likelihood computations can be done faster. 
}
}
\author{
Petter Mostad mostad@chalmers.se and Thore Egeland Thore.Egeland@gmail.com
}


\examples{
#A very simple specification
FamiliasLocus(1:4/10)
#A simple specification with more information but without mutations or silent alleles
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("A", "B", "C", "D"))
#Modified to include a silent frequency
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.3, 0.1), name="locus1", 
allelenames= c("08", "09", "10", "11", "silent"))
#Mutation rates added to default mutation model
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("08", "09", "10", "11"),femaleMutationRate=0.001,maleMutationRate=0.005)
#Mutation matrices specified directly 
femaleMM <- matrix(c(0.99, 0.005, 0.003, 0.002, 0.005, 0.99, 0.005, 0, 
0, 0.005, 0.99, 0.005, 0.002, 0.003, 0.005, 0.99), 4, 4, byrow=TRUE)
maleMM <- matrix(c(0.99, 0.005, 0.003, 0.002, 0.005, 0.99, 0.005, 0, 
0, 0.005, 0.99, 0.005, 0.002, 0.003, 0.005, 0.99), 4, 4, byrow=TRUE)
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("08", "09", "10", "11"), 
femaleMutationModel = "custom", femaleMutationMatrix = femaleMM, 
maleMutationModel = "custom", maleMutationMatrix = maleMM)
}