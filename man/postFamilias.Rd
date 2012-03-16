\name{postFamilias}
\alias{postFamilias}

\title{
Presents LRs
}
\description{
Postprocessing of output from \code{\link{GetProbabilities}}. Additional ouput compared to 
\code{\link{GetProbabilities}}: likelihood ratios (LRs) for each
marker and combined. logarithms are used to avoid numerical problems.
}
\usage{
postFamilias(res, ref = 1)
}
\arguments{
  \item{res}{
A \link{linkdat} object.
}
\item{ref}{
An integer determining the reference pedigree (denominator) in the LR.
}
}

\value{
  \item{LRperMarker}{LR for each marker}
   \item{LR}{LR based on all marker}
\item{res }{Copy of output from \code{\link{GetProbabilities}}}
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
}
