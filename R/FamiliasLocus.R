FamiliasLocus <- function (frequencies, name, allelenames, 
          femaleMutationModel  = "stable", maleMutationModel = "stable", 
          femaleMutationRate   = 0,        maleMutationRate  = 0, 
          femaleMutationRange  = 0.5,      maleMutationRange = 0.5, 
	  femaleMutationMatrix,            maleMutationMatrix) 
{
   if (missing(frequencies) || !is.numeric(frequencies) || any(frequencies<=0))
      stop("frequencies must be a vector of positive numbers")
   if (round(sum(frequencies), 6)!=1)
      stop("The frequencies must sum to 1")
   if (missing(name)) 
      name <- "locus"
   if (missing(allelenames))
      allelenames <- LETTERS[1:length(frequencies)]
   else
      if (length(allelenames)!=length(frequencies))
         stop("The number of allele names must be the same as the number of frequencies")   
   if (anyDuplicated(allelenames))
      stop("There cannot be duplicates among the allele names")
   if (any(allelenames[-length(allelenames)]=="silent") || any(allelenames[-length(allelenames)]=="Silent"))
      stop("Only the last allele can be specified as silent")
   names(frequencies) <- allelenames
   nAlleles <- length(frequencies)

   if ((femaleMutationModel != "equal" & femaleMutationModel != 
        "frequencies" & femaleMutationModel != "step" & femaleMutationModel != 
        "stable" & femaleMutationModel != "custom") | (maleMutationModel != "equal" & 
	maleMutationModel != "frequencies" & maleMutationModel != "step" & 
        maleMutationModel != "stable" & maleMutationModel != "custom")) 
      stop("ERROR: Mutation models must be either 'equal', 'frequencies', 'step', 'stable', or 'custom'")
   if (femaleMutationRate < 0 | maleMutationRate < 0) 
      stop("ERROR: Mutation rates cannot be negative")
   if (femaleMutationRange < 0 | maleMutationRange < 0) 
      stop("ERROR: Mutation ranges cannot be negative")  

   if ((femaleMutationRate > 0 & (femaleMutationModel == "step" | femaleMutationModel == "stable")) | 
       (maleMutationRate   > 0 & (maleMutationModel   == "step" | maleMutationModel   == "stable"  ))) {
       if (any(sort(allelenames)!=allelenames))
          stop(paste("ERROR: When mutation models 'step' or 'stable' are used with non-zero mutation rates,", 
                     "the resulting mutation matrix depends on the ordering of the alleles. To avoid errors and confusion,", 
                     "it is in these cases required that alleles are listed in lexicographical order."))
   }

# Parameter indicating if the probability of mutating to an allele is always independent of which allele it is mutating from: 
simpleMutationMatrices <- TRUE

# CHECKING OR PRODUCING FEMALE MUTATION MATRIX: 
   if (femaleMutationModel == "custom") {
      if (missing(femaleMutationMatrix)) 
         stop("When the female mutation model is 'custom' the female mutation matrix must be specified")
      if (!is.matrix(femaleMutationMatrix) | dim(femaleMutationMatrix)[1] != nAlleles | dim(femaleMutationMatrix)[2] != nAlleles)
         stop("The female mutation matrix must be of a dimension corresponding to the vector of frequencies")
      if (any(as.vector(femaleMutationMatrix)<0))
         stop("The female mutation matrix cannot have negative entries")
      if (any(round(apply(femaleMutationMatrix, 1, sum), 6) != 1))
         stop("The rows in the female mutation matrix must sum to 1")
      for (j in 1:nAlleles) {
         v <- femaleMutationMatrix[-j, j]
         if (any(round(v - v[1], 6)!=0)) simpleMutationMatrices <- FALSE
      }     
      femaleMutationType <- "A 'custom' specified mutation matrix"
   } else femaleMutationMatrix <- matrix(0, nAlleles, nAlleles)
   if (femaleMutationModel == "equal") {
      for (i in 1:nAlleles)
      	  for (j in 1:nAlleles) 
              if (i==j) femaleMutationMatrix[i,j] <- 1 - femaleMutationRate
              else femaleMutationMatrix[i,j] <- femaleMutationRate/(nAlleles - 1)
      if (femaleMutationRate > 0) 
         femaleMutationType <- paste("An 'equal' mutation model with mutation rate",femaleMutationRate)
      else
         femaleMutationType <- "No mutations"
   } else if (femaleMutationModel == "frequencies") {
      alpha <- femaleMutationRate/sum(frequencies*(1-frequencies))
      for (i in 1:nAlleles)
         for (j in 1:nAlleles)
            if (i==j) femaleMutationMatrix[i,j] <- 1 - alpha + alpha*frequencies[j]
            else femaleMutationMatrix[i,j] <- alpha*frequencies[j]
      if (femaleMutationRate > 0) 
         femaleMutationType <- paste("A 'frequencies' mutation model with expected mutation rate", femaleMutationRate)
      else
         femaleMutationType <- "No mutations"
   } else if (femaleMutationModel == "step") {
      for (i in 1:nAlleles) {
         k <- femaleMutationRate*(1-femaleMutationRange)/femaleMutationRange/
              (2-femaleMutationRange^(i-1)-femaleMutationRange^(nAlleles-i))
         for (j in 1:nAlleles)
            if (i==j) femaleMutationMatrix[i,j] <- 1 - femaleMutationRate
            else femaleMutationMatrix[i,j] <- k*femaleMutationRange^abs(i-j)
      }
      if (femaleMutationRate > 0) {
         simpleMutationMatrices <- FALSE
         femaleMutationType <- paste("A 'step' mutation model with expected mutation rate", femaleMutationRate, 
                         "and range parameter", femaleMutationRange) 
      } else femaleMutationType <- "No mutations"
   } else if (femaleMutationModel == "stable") {
      myconst <- femaleMutationRate*(1-femaleMutationRange)^2/2/
                 femaleMutationRange/(nAlleles - nAlleles*femaleMutationRange - 1 +
                 femaleMutationRange^nAlleles)
      for (i in 1:nAlleles)
         for (j in 1:nAlleles)
            if (i==j) femaleMutationMatrix[i,j] <- 1 - myconst*femaleMutationRange*(2 - 
                      femaleMutationRange^(i-1) - femaleMutationRange^(nAlleles-i))/
                      frequencies[i]/(1-femaleMutationRange)
            else  femaleMutationMatrix[i,j] <- myconst*femaleMutationRange^abs(i-j)/frequencies[i]
      if (femaleMutationRate > 0) {
         simpleMutationMatrices <- FALSE
         femaleMutationType <- paste("A 'stable' mutation model with expected mutation rate", femaleMutationRate, 
                      "and range parameter", femaleMutationRange)
      } else femaleMutationType <- "No mutations"
   }
   rownames(femaleMutationMatrix) <- allelenames
   colnames(femaleMutationMatrix) <- allelenames   

# CHECKING OR PRODUCING MALE MUTATION MATRIX: 
   if (maleMutationModel == "custom") {
      if (missing(maleMutationMatrix)) 
         stop("When the male mutation model is 'custom' the male mutation matrix must be specified")
      if (!is.matrix(maleMutationMatrix) | dim(maleMutationMatrix)[1] != nAlleles | dim(maleMutationMatrix)[2] != nAlleles)
         stop("The male mutation matrix must be of a dimension corresponding to the vector of frequencies")
      if (any(as.vector(maleMutationMatrix)<0))
         stop("The male mutation matrix cannot have negative entries")
      if (any(round(apply(maleMutationMatrix, 1, sum), 6) != 1))
         stop("The rows in the male mutation matrix must sum to 1")
      for (j in 1:nAlleles) {
         v <- maleMutationMatrix[-j, j]
         if (any(round(v - v[1], 6)!=0)) simpleMutationMatrices <- FALSE
      }     
      maleMutationType <- "A 'custom' specified mutation matrix"
   } else maleMutationMatrix <- matrix(0, nAlleles, nAlleles)
   if (maleMutationModel == "equal") {
      for (i in 1:nAlleles)
      	  for (j in 1:nAlleles) 
              if (i==j) maleMutationMatrix[i,j] <- 1 - maleMutationRate
              else maleMutationMatrix[i,j] <- maleMutationRate/(nAlleles - 1)
      if (maleMutationRate > 0) 
         maleMutationType <- paste("An 'equal' mutation model with mutation rate", maleMutationRate)
      else
         maleMutationType <- "No mutations"
   } else if (maleMutationModel == "frequencies") {
      alpha <- maleMutationRate/sum(frequencies*(1-frequencies))
      for (i in 1:nAlleles)
         for (j in 1:nAlleles)
            if (i==j) maleMutationMatrix[i,j] <- 1 - alpha + alpha*frequencies[j]
            else maleMutationMatrix[i,j] <- alpha*frequencies[j]
      if (maleMutationRate > 0) 
         maleMutationType <- paste("A 'frequencies' mutation model with expected mutation rate", maleMutationRate)
      else
         maleMutationType <- "No mutations"
   } else if (maleMutationModel == "step") {
      for (i in 1:nAlleles) {
         k <- maleMutationRate*(1-maleMutationRange)/maleMutationRange/
              (2-maleMutationRange^(i-1)-maleMutationRange^(nAlleles-i))
         for (j in 1:nAlleles)
            if (i==j) maleMutationMatrix[i,j] <- 1 - maleMutationRate
            else maleMutationMatrix[i,j] <- k*maleMutationRange^abs(i-j)
      }
      if (maleMutationRate > 0) {
         simpleMutationMatrices <- FALSE
         maleMutationType <- paste("A 'step' mutation model with expected mutation rate", maleMutationRate, 
                         "and range parameter", maleMutationRange) 
      } else maleMutationType <- "No mutations"
   } else if (maleMutationModel == "stable") {
      myconst <- maleMutationRate*(1-maleMutationRange)^2/2/
                 maleMutationRange/(nAlleles - nAlleles*maleMutationRange - 1 +
                 maleMutationRange^nAlleles)
      for (i in 1:nAlleles)
         for (j in 1:nAlleles)
            if (i==j) maleMutationMatrix[i,j] <- 1 - myconst*maleMutationRange*(2 - 
                      maleMutationRange^(i-1) - maleMutationRange^(nAlleles-i))/
                      frequencies[i]/(1-maleMutationRange)
            else  maleMutationMatrix[i,j] <- myconst*maleMutationRange^abs(i-j)/frequencies[i]
      if (maleMutationRate > 0) {
         simpleMutationMatrices <- FALSE
         maleMutationType <- paste("A 'stable' mutation model with expected mutation rate", maleMutationRate, 
                      "and range parameter", maleMutationRange)
      } else maleMutationType <- "No mutations"
   }
   rownames(maleMutationMatrix) <- allelenames
   colnames(maleMutationMatrix) <- allelenames   

   result <- list(locusname = name, 
                  alleles = frequencies, 
                  femaleMutationType = femaleMutationType, 
                  femaleMutationMatrix = femaleMutationMatrix, 
                  maleMutationType = maleMutationType, 
                  maleMutationMatrix = maleMutationMatrix, 
                  simpleMutationMatrices = simpleMutationMatrices) 
   class(result) <- "FamiliasLocus"
   result
}
