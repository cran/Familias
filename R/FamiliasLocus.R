FamiliasLocus <- function (frequencies, name, allelenames, 
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
{
    if (missing(frequencies)) 
        stop("The first argument must be a list of frequencies or a FamiliasLocus object")
    if(class(frequencies)=="FamiliasLocus") {
        if (!missing(name) || !missing(allelenames))
            stop("Only mutation parameters can be edited with the FamiliasLocus function")
        x <- frequencies
        frequencies <- x$alleles
        name        <- x$locusname
        allelenames <- names(frequencies)
    } else {
        if (!is.numeric(frequencies) || any(frequencies <= 0)) 
            stop("frequencies must be a vector of positive numbers")
        if (round(sum(frequencies), 6) != 1) 
            stop("The frequencies must sum to 1")
        if (missing(name)) 
            name <- deparse(substitute(frequencies))
        if (missing(allelenames)) {
            if (is.null(names(frequencies)))  
                allelenames <- as.character(1:length(frequencies))
            else 
                allelenames <- names(frequencies)
            }
        else if (length(allelenames) != length(frequencies)) 
            stop("The number of allele names must be the same as the number of frequencies")
        if (anyDuplicated(allelenames)) 
            stop("There cannot be duplicates among the allele names")
        if (any(allelenames[-length(allelenames)] == "silent") || 
            any(allelenames[-length(allelenames)] == "Silent")) 
            stop("Only the last allele can be specified as silent")
        names(frequencies) <- allelenames
    }
    nAlleles <- length(frequencies)
    hasSilentAllele <- (allelenames[nAlleles] == "silent" || allelenames[nAlleles] == "Silent")
    nAll <- nAlleles - hasSilentAllele
    freq <- frequencies[1:nAll]

    if (missing(femaleMutationModel)) femaleMutationModel <- MutationModel
    if (missing(femaleMutationRate)) femaleMutationRate <- MutationRate
    if (missing(femaleMutationRange)) femaleMutationRange <- MutationRange
    if (missing(femaleMutationRate2)) femaleMutationRate2 <- MutationRate2
    if (missing(femaleMutationMatrix) && !missing(MutationMatrix)) femaleMutationMatrix <- MutationMatrix
    if (missing(maleMutationModel)) maleMutationModel <- MutationModel
    if (missing(maleMutationRate)) maleMutationRate <- MutationRate
    if (missing(maleMutationRange)) maleMutationRange <- MutationRange
    if (missing(maleMutationRate2)) maleMutationRate2 <- MutationRate2
    if (missing(maleMutationMatrix) && !missing(MutationMatrix)) maleMutationMatrix <- MutationMatrix

    if ((femaleMutationModel != "Equal"        & femaleMutationModel != "equal" &
         femaleMutationModel != "Proportional" & femaleMutationModel != "proportional" & 
         femaleMutationModel != "Stepwise"     & femaleMutationModel != "stepwise" & 
         femaleMutationModel != "Custom"       & femaleMutationModel != "custom") | 
        (maleMutationModel   != "Equal"        & maleMutationModel   != "equal" &
         maleMutationModel   != "Proportional" & maleMutationModel   != "proportional" &
         maleMutationModel   != "Stepwise"     & maleMutationModel   != "stepwise" & 
         maleMutationModel   != "Custom"       & maleMutationModel   != "custom")) 
        stop("ERROR: Mutation models must be either 'Equal', 'Proportional', 'Stepwise', or 'Custom'")
    if (femaleMutationRate < 0 | maleMutationRate < 0 | femaleMutationRate > 1 | maleMutationRate > 1) 
        stop("ERROR: Mutation rates must be between 0 and 1.")
    if (femaleMutationRate2 < 0 | maleMutationRate2 < 0 | femaleMutationRate2 > 1 | maleMutationRate2 > 1) 
        stop("ERROR: Mutation rates must be between 0 and 1.")
    if (femaleMutationRange <= 0 | maleMutationRange <= 0) 
        stop("ERROR: Mutation ranges must be positive.")

    if (femaleMutationModel == "Custom" | femaleMutationModel == "custom") {
        if (missing(femaleMutationMatrix)) 
            stop("When the female mutation model is 'Custom' the female mutation matrix must be specified")
        if (!is.matrix(femaleMutationMatrix) | dim(femaleMutationMatrix)[1] != 
            nAlleles | dim(femaleMutationMatrix)[2] != nAlleles) 
            stop("The female mutation matrix must be of a dimension corresponding to the vector of frequencies")
        if (any(as.vector(femaleMutationMatrix) < 0)) 
            stop("The female mutation matrix cannot have negative entries")
        if (any(round(apply(femaleMutationMatrix, 1, sum), 6) != 1)) 
            stop("The rows in the female mutation matrix must sum to 1")
        femaleMutationType <- "A 'Custom' specified mutation matrix"
    }
    else 
    {
        femaleMutationMatrix <- matrix(0, nAlleles, nAlleles)
        diag(femaleMutationMatrix) <- 1
        if (femaleMutationRate == 0 ) 
            femaleMutationType <- "No mutations"
        else if (femaleMutationModel == "Equal" | femaleMutationModel == "equal") {
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                femaleMutationMatrix[i, j] <- 1 - femaleMutationRate
            else femaleMutationMatrix[i, j] <- femaleMutationRate/(nAll - 1)
            femaleMutationType <- paste("An 'Equal' mutation model with mutation rate", 
                femaleMutationRate)
        }
        else if (femaleMutationModel == "Proportional" | femaleMutationModel == "proportional") {
            sumfreq <- sum(freq)
            frq     <- freq/sumfreq
            alpha <- femaleMutationRate/sumfreq/sum(frq * (1 - frq))
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                femaleMutationMatrix[i, j] <- 1 - alpha + alpha * frq[j]
            else femaleMutationMatrix[i, j] <- alpha * frq[j]
            femaleMutationType <- paste("A 'Proportional' mutation model with expected mutation rate", 
                femaleMutationRate)
        }
        else if (femaleMutationModel == "Stepwise" | femaleMutationModel == "stepwise") {
            numfreq <- as.numeric(names(freq))
            if (any(is.na(numfreq)))
                stop("The 'Stepwise' mutation model requires all non-silent alleles to have numerical names.")
            if (any(round(numfreq, 1)!=numfreq))
                stop("Microvariants must be named as a decimal number with one decimal.")
            microgroup <- (numfreq - round(numfreq))*10
            for (i in 1:nAll) {
                microcompats <- (microgroup == microgroup[i])
                for (j in 1:nAll) {
                    if (i==j) {
                        if (all(microcompats)) femaleMutationMatrix[i,j] <- 1 - femaleMutationRate
                        else if (sum(microcompats)==1) femaleMutationMatrix[i,j] <- 1 - femaleMutationRate2
                        else femaleMutationMatrix[i,j] <- 1 - femaleMutationRate - femaleMutationRate2
                    } else if (microcompats[j])
                        femaleMutationMatrix[i,j] <- femaleMutationRange^abs(numfreq[i]-numfreq[j])
                    else
                        femaleMutationMatrix[i,j] <- femaleMutationRate2/(nAll-sum(microcompats))
                }
                microcompats[i] <- FALSE
                if (any(microcompats))
                    femaleMutationMatrix[i,microcompats] <- femaleMutationMatrix[i,microcompats]/
						            sum(femaleMutationMatrix[i,microcompats])*
						            femaleMutationRate
            }
            if (all(microgroup==0))
                femaleMutationType <- paste("A 'Stepwise' mutation model with mutation rate", 
                femaleMutationRate, "and mutation range", femaleMutationRange)
            else 
                femaleMutationType <- paste("A 'Stepwise' mutation model with mutation rate ", 
                femaleMutationRate, ", range ", femaleMutationRange, ", and fractional mutation rate ", 
		femaleMutationRate2, sep="")
        }
    }

    if (maleMutationModel == "Custom" | maleMutationModel == "custom") {
        if (missing(maleMutationMatrix)) 
            stop("When the male mutation model is 'Custom' the male mutation matrix must be specified")
        if (!is.matrix(maleMutationMatrix) | dim(maleMutationMatrix)[1] != 
            nAlleles | dim(maleMutationMatrix)[2] != nAlleles) 
            stop("The male mutation matrix must be of a dimension corresponding to the vector of frequencies")
        if (any(as.vector(maleMutationMatrix) < 0)) 
            stop("The male mutation matrix cannot have negative entries")
        if (any(round(apply(maleMutationMatrix, 1, sum), 6) != 1)) 
            stop("The rows in the male mutation matrix must sum to 1")
        maleMutationType <- "A 'Custom' specified mutation matrix"
    }
    else 
    {
        maleMutationMatrix <- matrix(0, nAlleles, nAlleles)
        diag(maleMutationMatrix) <- 1
        if (maleMutationRate == 0 ) 
            maleMutationType <- "No mutations"
        else if (maleMutationModel == "Equal" | maleMutationModel == "equal") {
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                maleMutationMatrix[i, j] <- 1 - maleMutationRate
            else maleMutationMatrix[i, j] <- maleMutationRate/(nAll - 1)
            maleMutationType <- paste("An 'Equal' mutation model with mutation rate", 
                maleMutationRate)
        }
        else if (maleMutationModel == "Proportional" | maleMutationModel == "proportional") {
            sumfreq <- sum(freq)
            frq     <- freq/sumfreq
            alpha <- maleMutationRate/sumfreq/sum(frq * (1 - frq))
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                maleMutationMatrix[i, j] <- 1 - alpha + alpha * frq[j]
            else maleMutationMatrix[i, j] <- alpha * frq[j]
            maleMutationType <- paste("A 'Proportional' mutation model with expected mutation rate", 
                maleMutationRate)
        }
        else if (maleMutationModel == "Stepwise" | maleMutationModel == "stepwise") {
            numfreq <- as.numeric(names(freq))
            if (any(is.na(numfreq)))
                stop("The 'Stepwise' mutation model requires all non-silent alleles to have numerical names.")
            if (any(round(numfreq, 1)!=numfreq))
                stop("Microvariants must be named as a decimal number with one decimal.")
            microgroup <- (numfreq - round(numfreq))*10
            for (i in 1:nAll) {
                microcompats <- (microgroup == microgroup[i])
                for (j in 1:nAll) {
                    if (i==j) {
                        if (all(microcompats)) maleMutationMatrix[i,j] <- 1 - maleMutationRate
                        else if (sum(microcompats)==1) maleMutationMatrix[i,j] <- 1 - maleMutationRate2
                        else maleMutationMatrix[i,j] <- 1 - maleMutationRate - maleMutationRate2
                    } else if (microcompats[j])
                        maleMutationMatrix[i,j] <- maleMutationRange^abs(numfreq[i]-numfreq[j])
                    else
                        maleMutationMatrix[i,j] <- maleMutationRate2/(nAll-sum(microcompats))
                }
                microcompats[i] <- FALSE
                if (any(microcompats))
                    maleMutationMatrix[i,microcompats] <- maleMutationMatrix[i,microcompats]/
						          sum(maleMutationMatrix[i,microcompats])*
						          maleMutationRate
            }
            if (all(microgroup==0))
                maleMutationType <- paste("A 'Stepwise' mutation model with mutation rate", 
                maleMutationRate, "and mutation range", maleMutationRange)
            else 
                maleMutationType <- paste("A 'Stepwise' mutation model with mutation rate ", 
                maleMutationRate, ", range ", maleMutationRange, ", and fractional mutation rate ", 
		maleMutationRate2, sep="")
        }
    }
    originalMaleMutationMatrix <- maleMutationMatrix
    originalFemaleMutationMatrix <- femaleMutationMatrix

    # do the stabilization 
    if (Stabilization == "Proportional" | Stabilization == "proportional") {
        stat <- eigen(t(femaleMutationMatrix))$vectors[,1]
        stat <- stat/sum(stat)
        dig  <- diag(femaleMutationMatrix)
        Dvec <- (1 - sum(dig*frequencies))/(1 - sum(dig*stat))*stat/frequencies
        femaleMutationMatrix <- diag(Dvec)%*%femaleMutationMatrix - diag(Dvec - 1)
        stat <- eigen(t(maleMutationMatrix))$vectors[,1]
        stat <- stat/sum(stat)
        dig  <- diag(maleMutationMatrix)
        Dvec <- (1 - sum(dig*frequencies))/(1 - sum(dig*stat))*stat/frequencies
        maleMutationMatrix <- diag(Dvec)%*%maleMutationMatrix - diag(Dvec - 1)
    } else if (Stabilization != "None" & Stabilization != "none")
        stop("The stabilization parameter must be 'None', 'Proportional', or 'Optimal'.")

    if ((Stabilization == "Proportional" | Stabilization == "proportional") & 
        (any(maleMutationMatrix < originalMaleMutationMatrix * StabilizationFactor) | 
         any(femaleMutationMatrix < originalFemaleMutationMatrix * StabilizationFactor)))
        stop("Stabilization algorithm failed")

    rownames(femaleMutationMatrix) <- names(frequencies)
    colnames(femaleMutationMatrix) <- names(frequencies)
    rownames(maleMutationMatrix) <- names(frequencies)
    colnames(maleMutationMatrix) <- names(frequencies)

    simpleMutationMatrices <- TRUE
    for (j in 1:nAlleles) {
        v <- femaleMutationMatrix[-j, j]
        if (any(round(v - v[1], 6) != 0)) 
            simpleMutationMatrices <- FALSE
    }
    for (j in 1:nAlleles) {
        v <- maleMutationMatrix[-j, j]
        if (any(round(v - v[1], 6) != 0)) 
            simpleMutationMatrices <- FALSE
    }
    result <- list(locusname = name, alleles = frequencies, femaleMutationType = femaleMutationType, 
        femaleMutationMatrix = femaleMutationMatrix, maleMutationType = maleMutationType, 
        maleMutationMatrix = maleMutationMatrix, simpleMutationMatrices = simpleMutationMatrices, 
	Stabilization = Stabilization)
    class(result) <- "FamiliasLocus"
    result
}

