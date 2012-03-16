\name{preFamilias}
\alias{preFamilias}
\title{
Prepares input for GetProbabilities
}
\description{
The function assumes that a \code{\link{linkdat}} object has been created using the R-package
paramlink. This could be done within R or by reading files on Merlin format using
the \code{\link{linkdat}} function of paramlink. The required calls to the \code{\link{Familias}} functions are then
performed. This function is normally called once for each pedigree. The parameter \code{new=TRUE}
for the first call and thereafter \code{new=FALSE}. The function has only been tested for the case
when individuals, markers and alleles are given integer names 1,2,...
}
\usage{
preFamilias(x, new = TRUE, mutationRateFemale = 0, mutationRateMale = 0, mutationModelFemale = "stable", mutationModelMale = "stable", mutationRangeFemale = 0.1, mutationRangeMale = 0.1, silentFrequency = NULL)
}

\arguments{
  \item{x}{
A \link{linkdat} object.
}
  \item{new}{
A logical variable set to TRUE first time function is called and 0 otherwise.
}
  \item{mutationRateFemale}{
See \code{\link{AddAlleleSystem}}
}
  \item{mutationRateMale}{
See \code{\link{AddAlleleSystem}}
}
  \item{mutationModelFemale}{
See \code{\link{AddAlleleSystem}}
}
  \item{mutationModelMale}{
See \code{\link{AddAlleleSystem}}
}
  \item{mutationRangeFemale}{
See \code{\link{AddAlleleSystem}}
}
  \item{mutationRangeMale}{
See \code{\link{AddAlleleSystem}}
}
  \item{silentFrequency}{
See \code{\link{AddAlleleSystem}}
}
}

\value{
\item{ped}{The number of the pedigree that gas been defined}
}
\references{
For more details about the parameters and the computations, 
see www.familias.name. 
}
\author{
Thore Egeland <Thore.Egeland@umb.no> 
}

\examples{
#Example 1. Same example as in the documentation of GetProbabilities.
require(paramlink)
require(Familias)
x1=swapSex(nuclearPed(1),1)
afreq1=c(0.2,0.3,0.5)
m1=marker(x1,1,c(1,1),2,c(2,2),3,c(1,2),afreq=afreq1,alleles=1:length(afreq1))
x1=addMarker(x1,m1)
afreq2=c(0.1,0.1,0.1,0.7)
m2=marker(x1,1,c(1,2),2,c(1,2),3,c(1,2),afreq=afreq2,alleles=1:length(afreq2))
x1=addMarker(x1,m2)
tt=c("System 1 freqs: 0.2,0.3,0.5. System 2 freqs: 0.1,0.1,0.1,0.7")
plot(x1,marker=1:2,title=tt)
preFamilias(x1,new=TRUE)
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res)

#Alternatively, the input can be read from files on Merlin format as below
write.linkdat(x1,prefix="Familias",what=c("ped","freq","dat","map")) #Skip if merlin files exist
x1=linkdat("Familias.ped",freq="Familias.freq",dat="Familias.dat",map="Familias.map") 
preFamilias(x1,new=TRUE)
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
GetProbabilities()


# Example 2. (Exercise S1)
# This and some following examples are based on Exercises in http://arken.umb.no/~theg/familiasExercises.pdf

x1=swapSex(nuclearPed(1),1)
afreq1=c(0.05,0.05,0.9)
m1=marker(x1,1,c(2,2),2,c(1,1),3,c(1,2),afreq=afreq1,alleles=1:length(afreq1))
x1=addMarker(x1,m1)
afreq2=c(0.1,0.1,0.1,0.7)
m2=marker(x1,1,c(2,3),2,c(1,1),3,c(1,3),afreq=afreq2,alleles=1:length(afreq2))
x1=addMarker(x1,m2)
tt=c("System 1 freqs: 0.2,0.3,0.5. System 2 freqs: 0.1,0.1,0.1,0.7")
plot(x1,marker=1:2,title=tt)
preFamilias(x1,new=TRUE)
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) #Answer in 1a-d

res=GetProbabilities(kinship=0.02)
postFamilias(res,ref=2) #Answer 1h

# Example 3 (Exercise S2. Mutation)
x1=swapSex(nuclearPed(1),1)
afreq1=c(0.072,0.082,0.212,0.292,0.222,0.097,0.02,0.003)
m1=marker(x1,2,c(14,15),3,c(16,17),afreq=afreq1,alleles=14:(14+length(afreq1)-1))
x1=addMarker(x1,m1)
plot(x1,marker=1,title="Exercise 2")
preFamilias(x1,new=TRUE,mutationRateFemale = 0.007, mutationRateMale = 0.007, mutationModelFemale = "equal", 
mutationModelMale = "equal")
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) # 

# Example 4 (Exercise S7. More mutations)
x1=swapSex(nuclearPed(1),1)
afreq1=c(0.072,0.082,0.212,0.292,0.222,0.097,0.02,0.003)
m1=marker(x1,2,c(14,15),3,c(16,17),afreq=afreq1,alleles=14:(14+length(afreq1)-1))
x1=addMarker(x1,m1)
plot(x1,marker=1,title="Exercise 2")
#a)
preFamilias(x1,new=TRUE,mutationRateFemale = 0.005, mutationRateMale = 0.005, mutationModelFemale = "equal", 
mutationModelMale = "equal")
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) # 

#b)
preFamilias(x1,new=TRUE,mutationRateFemale = 0.005, mutationRateMale = 0.005, mutationModelFemale = "frequencies", 
mutationModelMale = "frequencies")
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) # 

#c)
preFamilias(x1,new=TRUE,mutationRateFemale = 0.005, mutationRateMale = 0.005, mutationModelFemale = "step", 
mutationModelMale = "step",mutationRangeFemale = 0.5, mutationRangeMale = 0.5)
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) # 

#d)
preFamilias(x1,new=TRUE,mutationRateFemale = 0.005, mutationRateMale = 0.005, mutationModelFemale = "stable", 
mutationModelMale = "stable",mutationRangeFemale = 0.5, mutationRangeMale = 0.5)
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) # 

# Example 5 (Exercise S11. Silent alleles)
x1=nuclearPed(1)
afreq1=c(0.1,0.1,0.75)
m=marker(x1,1,c(2,2),2,c(1,1),3,c(1,1),afreq=afreq1,alleles=1:length(afreq1)) #Warning, ignore
nn=c("AF","mother","child")
x1=addMarker(x1,m)
plot(x1,marker=1,id.labels=nn,title="Figure 5. Silent allele?")
preFamilias(x1,new=TRUE,silentFrequency=0.05)
x2=x1
x2$pedigree[3,2]=0
preFamilias(x2,new=FALSE)
res=GetProbabilities()
postFamilias(res,ref=2) # 

# Example 6. Simulations
x = halfCousinPed(1);plot(x)
x=swapSex(x,c(4,6,8))
plot(x)
LR=NULL
for (i in 1:10){
y=SNPsim(x,N=100,available=c(1,8:9),seed=i)
preFamilias(y)
AddPedigree()
res=GetProbabilities()
post=postFamilias(res,ref=2)$LR[1]
LR=c(LR,post)
}

}
