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
FamiliasLocus(frequencies, name, allelenames, 
			   MutationModel = "Stepwise", 
                           MutationRate  = 0, 
                           MutationRange = 0.5, 
                           MutationRate2 = 0, 
			   MutationMatrix, 
			   Stabilization = "None", 
			   StabilizationFactor = 0.1, 
			   femaleMutationModel, 
			   femaleMutationRate, 
			   femaleMutationRange, 
			   femaleMutationRate2, 
			   femaleMutationMatrix, 
			   maleMutationModel, 
			   maleMutationRate, 
			   maleMutationRange, 
			   maleMutationRate2, 
			   maleMutationMatrix)
}

\arguments{
  \item{frequencies}{
The first input of FamiliasLocus may be either a vector containing allele frequencies, or a previously created FamiliasLocus object. 
In the first case, the vector may include a possible silent allele; that it is silent is indicated in the allelenames. 
The frequencies must sum to 1. 
In the second case, the new object will be identical to the old object in terms of frequencies, names of alleles, and name of locus, 
whereas the mutation matrices will be based on the input given to the function. 
}
  \item{name}{
Characters like 'D3S1358', used to identify the locus (marker). The default is to use the name of the frequencies argument to this function. 
}
  \item{allelenames}{
Names of the alleles, like \code{15} or 'A'. Note that thel last allele may be called 'Silent' (or 'silent'). 
It is then treated as a silent allele in subsequent likelihood calculations. The default is to use the names attribute of 
the frequencies, if it exists; otherwise the default is to use consecutive integers, starting at 1.
Note that if the 'Stepwise' mutation model is used, allele names (except for a silent allele) must be integers, 
with microvariants named as for example 15.2. 
}
  \item{MutationModel}{
The mutation model, used to create the mutation matrix. It may be 'Equal', 'Proportional', 'Stepwise', or 'Custom', see Details. 
}
  \item{MutationRate}{
The mutation rate; for the 'Stepwise' model the rate of integer-step mutations. It is not used when the MutationModel is 'Custom'.
}
  \item{MutationRange}{
Only used when the MutationModel is 'Stepwise'. It then indicates the relative probability of mutating n+1 steps versus 
mutating n steps. 
}
  \item{MutationRate2}{
Only used when the MutationModel is 'Stepwise'. It then indicates the rate of non-integer-step mutations, e.g., mutations 
from an allele with an integer name to alleles with decimal names indicating microvariants. 
}
  \item{MutationMatrix}{
Only used when the MutationModel is 'Custom'. It then directly specifies the mutation matrix. 
}
  \item{Stabilization}{
This parameter may currently have values 'None' or 'Proportional'. If the value is 'Proportional', 
the male and female mutation matrices derived from on the other mutation parameters are stabilized, in the sense
that allele frequencies after one or more generation of mutations will be equal to the original allele frequencies. 
The 'Proportional' method will keep the relative proportions of different mutation probabilities from one allele to other alleles. 
}
  \item{StabilizationFactor}{
Not used when stabilization is 'None'. Otherwise it indicates a lower bound for the factor an element of the mutation matrix 
may be changed with in a stabilization. 
}
  \item{femaleMutationModel}{
Specifies a separate female value for MutationModel; defaults to MutationModel. 
}
  \item{femaleMutationRate}{
Specifies a separate female value for MutationRate; defaults to MutationRate. 
}
  \item{femaleMutationRange}{
Specifies a separate female value for MutationRange; defaults to MutationRange. 
}
  \item{femaleMutationRate2}{
Specifies a separate female value for MutationRate2; defaults to MutationRate2. 
}
  \item{femaleMutationMatrix}{
Specifies a separate female value for MutationMatrix; defaults to MutationMatrix. 
}
  \item{maleMutationModel}{
Specifies a separate male value for MutationModel; defaults to MutationModel. 
}
  \item{maleMutationRate}{
Specifies a separate male value for MutationRate; defaults to MutationRate. 
}
  \item{maleMutationRange}{
Specifies a separate male value for MutationRange; defaults to MutationRange. 
}
  \item{maleMutationRate2}{
Specifies a separate male value for MutationRate2; defaults to MutationRate2. 
}
  \item{maleMutationMatrix}{
Specifies a separate male value for MutationMatrix; defaults to MutationMatrix. 
}
}
\details{
The probabilities for when and how mutations happen can be specified in mutation matrices, where 
the row corresponding to an allele indicates the probabilities that the allele is transferred as 
the allele indicated by the column. Mutation matrices may be specified directly 
in the MutationMatrix parameters by setting the 
value of the MutationModel parameter to 'Custom'. Otherwise they are computed based on the values 
of the MutationModel, MutationRate, MutationRate2, and MutationRange parameters. 
If MutationModel is 'Equal', there is an equal probability of mutating to any non-silent allele, given that a mutation happens.   
This model is referred to as 'Equal probability (simple and fast)' in Familias 2.0. 
If MutationModel is 'Proportional', the probability of mutating to any non-silent allele is proportional to its frequency. 
It is referred to as 'Probability proportional to frequency (stable)' in Familias 2.0.
If MutationModel is 'Stepwise', it is required that the names of all non-silent alleles are positive integers,
indicating the number of sequence repetitions of an STR marker, or decimal numbers with a single decimal, such as '15.2', 
indicating a microvariant. Mutations are then divided into two types: Those that 
add or subtract an integer to the allele, and those that add or subtract some fractional amount. 
The rate of these two types of mutations are given separately as MutationRate and MutationRate2, respectively. 
Relative probabilities of different mutaitons of the first type are specified using the MutationRange parameter. 
The model with only integer alleles is referred to as 'Probability decreasing with range (equal)' in Familias 2.0, 
while the more general model is called 'Extended stepwise' in Familias 3.0. 
Note that the probability of mutations to or from silent alleles is set to zero in all models except the 'Custom' model. 

The 'Proportional' setting for the stabilization parameter may in principle be used 
together with all the other settings of the mutation parameters, thereby changing both the male and 
female mutation matrices. When the 'Proportional' setting is used together with the 'Stepwise' MutationModel 
and all allele names are integer, the resulting model is referred to as 
'Probability decreasing with range (stable)' in Familias 2.0. 
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
If this is true, some likelihood computations can be done faster.}
\item{Stabilization}{The stabilization method used.}
}
\author{
Petter Mostad mostad@chalmers.se and Thore Egeland Thore.Egeland@gmail.com
}

\examples{
#Simple examples
FamiliasLocus(1:4/10)
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("A", "B", "C", "D"))

#Modified to include a silent frequency
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.3, 0.1), name="locus1", 
allelenames= c("8", "9", "10", "11", "silent"))

#Mutation rates added
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("8", "9", "10", "11"), femaleMutationRate=0.001, maleMutationRate=0.005)

#Mutation matrices specified directly 
MM <- matrix(c(0.99, 0.005, 0.003, 0.002, 0.005, 0.99, 0.005, 0, 
0, 0.005, 0.99, 0.005, 0.002, 0.003, 0.005, 0.99), 4, 4, byrow=TRUE)
FamiliasLocus(frequencies=c(0.1, 0.2, 0.3, 0.4), name="locus1", 
allelenames= c("08", "09", "10", "11"), MutationModel = "Custom", 
MutationMatrix = MM)

#A locus is first created, and then edited
loc <- FamiliasLocus(c(0.2, 0.5, 0.3))
FamiliasLocus(loc, maleMutationRate = 0.001)

#A locus using standard Norwegian frequencies is created
data(NorwegianFrequencies)
FamiliasLocus(NorwegianFrequencies$TH01)
}