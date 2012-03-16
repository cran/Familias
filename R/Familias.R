#This file should contain the R interface functions: 

GetNumberOfSystems <- function()
{
   result <- .C("GetNumberOfSystems", nsystems = integer(1))
   result$nsystems
}

GetProbabilities <- function(
                   kinship = 0, 
                   generationsParameter = 1,  
		   inbreedingParameter = 1, 
		   partnerParameter = 1,
                   onlyprior = FALSE, 
		   maxGenerations = NULL) 
{
if (kinship < 0 | generationsParameter < 0 | inbreedingParameter < 0 | partnerParameter < 0) 
   cat("ERROR: The parameters cannot be negative.\n")
else {
nsystems <- GetNumberOfSystems()
npedigrees <- GetNumberOfPedigrees()
if (is.null(maxGenerations)) maxGenerations <- -1
result <- .C("GetProbabilities", 
       	  as.double(generationsParameter), 
          as.integer(maxGenerations), 
    	  as.double(inbreedingParameter), 
    	  as.double(partnerParameter), 
    	  as.integer(!onlyprior),
	  as.double(kinship), 
	  redundant = integer(npedigrees), 
	  probabilities = double(npedigrees), 
	  likelihoods = double(nsystems*npedigrees), 
	  error = integer(1))
if (result$error==1) 
   cat("ERROR: Wrong input.\n")
else if (result$error==2)
   cat("ERROR: All pedigrees have probability zero.\n")
else
   prior <- result$probabilities 
   likelihoodsPerSystem <- matrix(result$likelihoods, 
     GetNumberOfSystems(), GetNumberOfPedigrees())
   likelihoods <- apply(likelihoodsPerSystem, 2, prod)
   posterior <- prior * likelihoods
   posterior <- posterior/sum(posterior)
   pedigreeUnique <- !as.logical(result$redundant)
   list(posterior = posterior, 
        prior = prior,  
        likelihoods = likelihoods,
        likelihoodsPerSystem = likelihoodsPerSystem, 
	pedigreeKept = pedigreeUnique)
}
}

RemoveDNAObservation <- function(indexperson, 
			indexAlleleSystem)
{
result <- .C("RemoveDNAObservation", as.integer(indexperson - 1), 
			   as.integer(indexAlleleSystem - 1), 
			   error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
}

AddDNAObservation <- function(indexperson, 
		  	indexAlleleSystem, 
			indexAllele1, 
			indexAllele2 = indexAllele1)
{
result <- .C("AddDNAObservation", as.integer(indexperson - 1), 
			as.integer(indexAlleleSystem - 1), 
			as.integer(indexAllele1 - 1), 
                        as.integer(indexAllele2 - 1), 
			error = integer(1))		  
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
} 

EditAlleleSystem <- function(indexSystem, 
                         frequencies,
                         correspondence = 1:length(frequencies), 
			 mutationRateFemale = 0,
			 mutationRateMale = 0,
			 mutationModelFemale = "stable", 
			 mutationModelMale = "stable", 
			 mutationRangeFemale = 0.1, 
			 mutationRangeMale = 0.1, 
			 silentFrequency = NULL)
{
nAlleles <- length(frequencies)   
if (!is.null(silentFrequency)) {
   frequencies <- c(frequencies, silentFrequency)
   correspondence <- c(correspondence, length(correspondence) + 1)
   nAlleles <- nAlleles + 1
}
if (any(frequencies<0))
   cat("ERROR: Allele frequencies cannot be negative.\n")
else if (sum(frequencies)>1)
   cat("ERROR: Allele frequencies cannot sum to more than 1.\n")
else if (mutationRateFemale < 0 | mutationRateMale < 0)
   cat("ERROR: Mutation rates cannot be negative.\n")
else if ((mutationModelFemale!="equal" & mutationModelFemale!="frequencies" & mutationModelFemale!="step" & mutationModelFemale!="stable") |
         (mutationModelMale  !="equal" & mutationModelMale  !="frequencies" & mutationModelMale  !="step" & mutationModelMale  !="stable"))
   cat("ERROR: Mutation models must be either \"equal\", \"frequencies\", \"step\", or \"stable\".\n")
else if (mutationRangeFemale < 0 | mutationRangeMale < 0)
   cat("ERROR: Mutation ranges cannot be negative.\n")
else if (length(correspondence)!=length(frequencies))
   cat("ERROR: The correspondence vector must be of the same length as the frequencies vector.\n")
else {
   if (mutationModelFemale=="equal") mutationModelFemale <- 0
   else if (mutationModelFemale=="frequencies") mutationModelFemale <- 1
   else if (mutationModelFemale=="step") mutationModelFemale <- 2
   else if (mutationModelFemale=="stable") mutationModelFemale <- 3

   if (mutationModelMale=="equal") mutationModelMale <- 0
   else if (mutationModelMale=="frequencies") mutationModelMale <- 1
   else if (mutationModelMale=="step") mutationModelMale <- 2
   else if (mutationModelMale=="stable") mutationModelMale <- 3

   correspondence[is.na(correspondence)] <- 0 

   nPossibilities <- nAlleles #Dummy not really used

result <- .C("EditAlleleSystem", as.integer(indexSystem - 1), 
		       as.integer(nAlleles), 
                       as.double(mutationRateFemale), 
		       as.double(mutationRateMale), 
                       as.integer(mutationModelFemale), 
                       as.integer(mutationModelMale), 
                       as.integer(nPossibilities), 
                       as.double(mutationRangeFemale), 
                       as.double(mutationRangeMale), 
                       as.double(frequencies), 
                       as.integer(correspondence - 1), 
                       as.integer(!is.null(silentFrequency)),
		       error = integer(1))
if (result$error>1) 
   cat("ERROR: Wrong input.\n")	
}		  
} 

RemoveAlleleSystem <- function(index)
{
result <- .C("RemoveAlleleSystem", as.integer(index - 1), 
			  error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
}

AddAlleleSystem <- function(frequencies, 
		        mutationRateFemale = 0, 
			mutationRateMale = 0, 
			mutationModelFemale = "stable", 
			mutationModelMale = "stable", 
			mutationRangeFemale = 0.1, 
			mutationRangeMale = 0.1, 
			silentFrequency = NULL)
{
nAlleles <- length(frequencies)   
if (!is.null(silentFrequency)) {
   frequencies <- c(frequencies, silentFrequency)
   nAlleles <- nAlleles + 1
}
if (any(frequencies<0))
   cat("ERROR: Allele frequencies cannot be negative.\n")
else if (sum(frequencies)>1)
   cat("ERROR: Allele frequencies cannot sum to more than 1.\n")
else if (mutationRateFemale < 0 | mutationRateMale < 0)
   cat("ERROR: Mutation rates cannot be negative.\n")
else if ((mutationModelFemale!="equal" & mutationModelFemale!="frequencies" & mutationModelFemale!="step" & mutationModelFemale!="stable") |
         (mutationModelMale  !="equal" & mutationModelMale  !="frequencies" & mutationModelMale  !="step" & mutationModelMale  !="stable"))
   cat("ERROR: Mutation models must be either \"equal\", \"frequencies\", \"step\", or \"stable\".\n")
else if (mutationRangeFemale < 0 | mutationRangeMale < 0)
   cat("ERROR: Mutation ranges cannot be negative.\n")
else {
   if (mutationModelFemale=="equal") mutationModelFemale <- 0
   else if (mutationModelFemale=="frequencies") mutationModelFemale <- 1
   else if (mutationModelFemale=="step") mutationModelFemale <- 2
   else if (mutationModelFemale=="stable") mutationModelFemale <- 3

   if (mutationModelMale=="equal") mutationModelMale <- 0
   else if (mutationModelMale=="frequencies") mutationModelMale <- 1
   else if (mutationModelMale=="step") mutationModelMale <- 2
   else if (mutationModelMale=="stable") mutationModelMale <- 3

   nPossibilities <- nAlleles #Dummy not really used
   result <- .C("AddAlleleSystem", as.integer(nAlleles), 
    		      as.double(mutationRateFemale), 
		      as.double(mutationRateMale), 
		      as.integer(mutationModelFemale), 
		      as.integer(mutationModelMale),
		      as.integer(nPossibilities), 
		      as.double(mutationRangeFemale), 
		      as.double(mutationRangeMale), 
		      as.double(frequencies), 
		      as.integer(!is.null(silentFrequency)), 
		      index = integer(1), 
		      error = integer(1))
   if (result$error>0) 
      cat("ERROR: Wrong input.\n")			  
   else
      result$index + 1	
   }	     
}

RemoveRelation <- function(parentindex, 
		       childindex, 
		       pedigree)
{
result <- .C("RemoveRelation", as.integer(parentindex-1), 
		     as.integer(childindex-1), 
		     as.integer(pedigree-1), 
		     error = integer(1))
if (result$error==1) 
   cat("ERROR: Wrong input.\n")	
else if (result$error==2)
   cat("ERROR: Trying to remove a fixed relation.\n")		  
}

AddRelation <- function(parentindex, 
		    childindex, 
		    pedigree)
{
result <- .C("AddRelation", as.integer(parentindex-1), 
		  as.integer(childindex-1), 
		  as.integer(pedigree-1), 
		  error = integer(1))
if (result$error==1) 
   cat("ERROR: Wrong input.\n")	
else if (result$error==2)
   cat("ERROR: Illegal relation based on Year-of-birth or is-Child data.\n")
else if (result$error==3)
   cat("ERROR: Cycle in the pedigree or duplicate parent.\n")
} 

AddExtraPerson <- function(male, 
		pedigree) 
{
result <- .C("AddExtraPerson", as.integer(male), 
		     as.integer(pedigree - 1), 
		     error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
} 

RemoveExtraPerson <- function(person, 
			  pedigree)
{
result <- .C("RemoveExtraPerson", as.integer(person - 1), 
			as.integer(pedigree - 1), 
			error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
} 

GetPedigree <- function(index)
{
size <- GetSizeOfPedigree(index)
result <- .C("GetPedigree", as.integer(index-1), 
		  matrix = integer(size*size), 
		  error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else
   matrix(as.logical(result$matrix), size, size)
} 

GetParents <- function(index)
{
size <- GetSizeOfPedigree(index)
result <- .C("GetParents", as.integer(index-1), 
		 mother = integer(size), 
		 father = integer(size), 
		 error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else {
   result$mother[result$mother==-1] <- NA
   result$father[result$father==-1] <- NA
   list(mother = result$mother+1, father = result$father+1)
}
}

GetSizeOfPedigree <- function(index)
{
result <- .C("GetSizeOfPedigree", as.integer(index - 1), 
			size = integer(1),
			error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else
   result$size
}

GetNumberOfExtraFemales <- function(pedigree)
{
result <- .C("GetNumberOfExtraFemales", as.integer(pedigree-1), 
			       number = integer(1), 
			       error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else
   result$number
}

GetNumberOfExtraMales <- function(pedigree)
{
result <- .C("GetNumberOfExtraMales", as.integer(pedigree-1), 
			    number = integer(1), 
			    error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else 
   result$number
} 

RemovePedigree <- function(index)
{
result <- .C("RemovePedigree", as.integer(index - 1), 
		     error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
} 

AddPedigree <- function(nExtraFemales = 0, 
		   nExtraMales = 0)
{
result <- .C("AddPedigree", as.integer(nExtraFemales), 
		  as.integer(nExtraMales), 
		  index = integer(1), 
		  error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else
   result$index + 1
}

GetNumberOfPedigrees <- function()
{
result <- .C("GetNumberOfPedigrees", number = integer(1))
result$number
}

GeneratePedigrees <- function(nExtraFemales = 0, 
			  nExtraMales = 0)
{
npedigrees <- GetNumberOfPedigrees()
result <- .C("GeneratePedigrees", as.integer(nExtraFemales), 
			as.integer(nExtraMales), 
			removed = integer(npedigrees), 
			error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
else
   !as.logical(result$removed)
} 

RemoveFixedRelation <- function(parentindex, 
			    childindex)
{
result <- .C("RemoveFixedRelation", as.integer(parentindex-1), 
			  as.integer(childindex-1), 
			  error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")			  
}

AddFixedRelation <- function(parentindex, 
			 childindex) 
{
npedigrees <- GetNumberOfPedigrees()
result <- .C("AddFixedRelation", as.integer(parentindex-1), 
		       as.integer(childindex-1), 
		       removed = integer(npedigrees), 
		       error = integer(1))
if (result$error==1) 
   cat("ERROR: Wrong input.\n")	
else if (result$error==2)
   cat("ERROR: Illegal relation based on Year-of-birth or is-Child data.\n")
else if (result$error==3)
   cat("ERROR: Illegal relation: A child cannot be its own parent!\n")
else if (result$error==4)
   cat("ERROR: Illegal relation: Child is ancestor of parent.\n")
else if (result$error==5)
   cat("ERROR: Illegal relation: Child already has such parent.\n")
else
   !as.logical(result$removed) 
}

RemovePerson <- function(index) 
{
result <- .C("RemovePerson", as.integer(index - 1), 
		   error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")	
}

AddPerson <- function(male, 
                      isChild = FALSE, 
		      yearOfBirth=NULL)
{
if (is.null(yearOfBirth)) yearOfBirth = -1
result <- .C("AddPerson", 
                as.integer(male), 
		as.integer(yearOfBirth), 
		as.integer(isChild), 
		index = integer(1), 
		error = integer(1))
if (result$error>0) 
   cat("ERROR: Wrong input.\n")	
else 
   result$index + 1
} 

NewFamilias <- function() 
{
.C("NewFamilias")
   cat("All information stored in Familias was removed.\n")
}