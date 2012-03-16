postFamilias <-
function(res,ref=1){
noOfPedigrees=length(res$pedigreeKept)
if(ref<1|ref>noOfPedigrees) stop("Impossible reference pedigree")
if(noOfPedigrees==1) return(res)
else{
x=res$likelihoodsPerSystem
LReach= x/x[,ref]
logLik=log(res$likelihoodsPerSystem)
LR=exp(apply(logLik-logLik[,ref],2,sum))
list(LRperMarker=LReach,LR=LR,res=res)
}
}
