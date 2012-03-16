preFamilias <-
function(x,new=TRUE, mutationRateFemale = 0, mutationRateMale = 0, mutationModelFemale = "stable", 
mutationModelMale = "stable", mutationRangeFemale = 0.1, mutationRangeMale = 0.1, silentFrequency = NULL) {
if(new){
 NewFamilias()
 for (i in x$orig.ids)
 AddPerson(x$pedigree[i,4]==1)
 for (i in 1:x$nMark)
 AddAlleleSystem(attr(x$markerdata[i][[1]],which="afreq"),mutationRateFemale = mutationRateFemale, 
mutationRateMale = mutationRateMale, mutationModelFemale = mutationModelFemale, 
mutationModelMale = mutationModelMale, mutationRangeFemale = mutationRangeFemale, mutationRangeMale = mutationRangeMale, 
silentFrequency = silentFrequency) 
 for (i in 1:x$nMark)
 for (j in x$available)
 AddDNAObservation(j, i, x$markerdata[[i]][j,1],x$markerdata[[i]][j,2]) 
}

ped=AddPedigree() 
if(max(x$ped[,2:3])>0){
for (i in x$orig.ids){
Father=x$pedigree[i,2]
Mother=x$pedigree[i,3]
Child=x$pedigree[i,1]
if(Father>0) AddRelation(Father, Child, ped)
if(Mother>0) AddRelation(Mother, Child, ped)
}
}
ped
}
