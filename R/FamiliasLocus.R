FamiliasLocus <- function (frequencies, allelenames, name, 
			   MutationModel = "Stepwise", 
                     MutationRate  = 0, 
                     MutationRange = 0.5, 
                     MutationRate2 = 0, 
			   MutationMatrix, 
			   Stabilization = "None", 
			   MaxStabilizedMutrate = 1,
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
        stop("The first argument must be a list of frequencies or a FamiliasLocus object.")
    if(class(frequencies)=="FamiliasLocus") {
        if (!missing(name) || !missing(allelenames))
            stop("Only mutation parameters can be edited with the FamiliasLocus function.")
        x <- frequencies
        frequencies <- x$alleles
        name        <- x$locusname
        allelenames <- names(frequencies)
        if (missing(MutationModel) & missing(MutationRate) & 
            missing(MutationRange) & missing(MutationRate2) & 
            missing(MutationMatrix) & missing(femaleMutationModel) & 
            missing(femaleMutationRate) & missing(femaleMutationRange) & 
            missing(femaleMutationRate2) & missing(femaleMutationMatrix) & 
            missing(maleMutationModel) & missing(maleMutationRate) & 
            missing(maleMutationRange) & missing(maleMutationRate2) & 
            missing(maleMutationMatrix)) {
            if (missing(Stabilization)) 
                stop("When a FamiliasLocus object is given, at least one mutation parameter or 'Stabilization' must be specified.")
            MutationModel <- "Custom"
            maleMutationMatrix <- x$maleMutationMatrix
            femaleMutationMatrix <- x$femaleMutationMatrix
        }
    } else {
        if (!is.numeric(frequencies) || any(frequencies <= 0)) 
            stop("frequencies must be a vector of positive numbers.")
        if (round(sum(frequencies), 6) != 1) 
            stop("The frequencies must sum to 1.")
        if (missing(name)) 
            name <- deparse(substitute(frequencies))
        if (missing(allelenames)) {
            if (is.null(names(frequencies)))  
                allelenames <- as.character(1:length(frequencies))
            else 
                allelenames <- names(frequencies)
            }
        else if (length(allelenames) != length(frequencies)) 
            stop("The number of allele names must be the same as the number of frequencies.")
        if (anyDuplicated(allelenames)) 
            stop("There cannot be duplicates among the allele names.")
        if (any(allelenames[-length(allelenames)] == "silent") || 
            any(allelenames[-length(allelenames)] == "Silent")) 
            stop("Only the last allele can be specified as silent.")
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
        stop("ERROR: Mutation models must be either 'Equal', 'Proportional', 'Stepwise', or 'Custom'.")
    if (femaleMutationRate < 0 | maleMutationRate < 0 | femaleMutationRate > 1 | maleMutationRate > 1) 
        stop("ERROR: Mutation rates must be between 0 and 1.")
    if (femaleMutationRate2 < 0 | maleMutationRate2 < 0 | femaleMutationRate2 > 1 | maleMutationRate2 > 1) 
        stop("ERROR: Mutation rates must be between 0 and 1.")
    if (femaleMutationRange <= 0 | maleMutationRange <= 0) 
        stop("ERROR: Mutation ranges must be positive.")

    if (femaleMutationModel == "Custom" | femaleMutationModel == "custom") {
        if (missing(femaleMutationMatrix)) 
            stop("When the female mutation model is 'Custom' the female mutation matrix must be specified.")
        if (!is.matrix(femaleMutationMatrix) | dim(femaleMutationMatrix)[1] != 
            nAlleles | dim(femaleMutationMatrix)[2] != nAlleles) 
            stop("The female mutation matrix must be of a dimension corresponding to the vector of frequencies.")
        if (any(as.vector(femaleMutationMatrix) < 0)) 
            stop("The female mutation matrix cannot have negative entries.")
        if (any(round(apply(femaleMutationMatrix, 1, sum), 6) != 1)) 
            stop("The rows in the female mutation matrix must sum to 1.")
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
            stop("When the male mutation model is 'Custom' the male mutation matrix must be specified.")
        if (!is.matrix(maleMutationMatrix) | dim(maleMutationMatrix)[1] != 
            nAlleles | dim(maleMutationMatrix)[2] != nAlleles) 
            stop("The male mutation matrix must be of a dimension corresponding to the vector of frequencies.")
        if (any(as.vector(maleMutationMatrix) < 0)) 
            stop("The male mutation matrix cannot have negative entries.")
        if (any(round(apply(maleMutationMatrix, 1, sum), 6) != 1)) 
            stop("The rows in the male mutation matrix must sum to 1.")
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

    # do the stabilization 
    Stabilization <- toupper(Stabilization)
    if (Stabilization != "NONE" & Stabilization != "DP" & 
        Stabilization != "RM" & Stabilization != "PM") 
        stop("The Stabilization parameter must be 'None', 'DP', 'RM', or 'PM'.")
    if (Stabilization != "NONE") {
        if (hasSilentAllele & all(maleMutationMatrix[1:nAll,nAlleles]==0) & 
                              all(maleMutationMatrix[nAlleles,1:nAll]==0) & 
                              all(femaleMutationMatrix[1:nAll,nAlleles]==0) & 
                              all(femaleMutationMatrix[nAlleles,1:nAll]==0)) {
            res <- stabilize(femaleMutationMatrix[1:nAll,1:nAll], frequencies[1:nAll]/sum(frequencies[1:nAll]), 
                             Stabilization, MaxStabilizedMutrate) 
            print(paste("Female mutation matrix f ratio:", res$fratio))
            print(paste("Female mutation matrix max specific mutation rate:", res$maxrate))
            femaleMutationMatrix[1:nAll,1:nAll] <- res$stabilized
            res <- stabilize(maleMutationMatrix[1:nAll,1:nAll], frequencies[1:nAll]/sum(frequencies[1:nAll]), 
                             Stabilization, MaxStabilizedMutrate) 
            print(paste("Male mutation matrix f ratio:", res$fratio))
            print(paste("Male mutation matrix max specific mutation rate:", res$maxrate))
            maleMutationMatrix[1:nAll,1:nAll] <- res$stabilized
        } else {
            res <- stabilize(femaleMutationMatrix, frequencies, 
                             Stabilization, MaxStabilizedMutrate) 
            print(paste("Female mutation matrix f ratio:", res$fratio))
            print(paste("Female mutation matrix max specific mutation rate:", res$maxrate))
            femaleMutationMatrix <- res$stabilized
            res <- stabilize(maleMutationMatrix, frequencies, 
                             Stabilization, MaxStabilizedMutrate) 
            print(paste("Male mutation matrix f ratio:", res$fratio))
            print(paste("Male mutation matrix max specific mutation rate:", res$maxrate))
            maleMutationMatrix <- res$stabilized
        }
    }
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

stabilize <- function(M,pe,stabilizationMethod="DP",t=1){
  #library('Rsolnp')
  R = 1-sum(diag(M)*pe)
  n = dim(M)[1]
  m = n^2
  xM = as.vector(M)
  tol = 1e-10
  if (all(abs(M - diag(n))<tol)) {
    P = M
  } else if (stabilizationMethod == "DP"){
    if (any(xM==0))
      stop("DP stabilization not possible unless all mutation matrix elements are positive.") 
    C = matrix(0,3*n-1,m)
    for (i in 1:n){
      C[seq(1,n),seq(n*(i-1)+1,n*i)] = diag(n)
      C[n+i,seq(n*(i-1)+1,n*i)] = pe
      C[2*n,seq(n*(n-1)+1,m)] = 0
      C[2*n+i-1,n*(i-1)+i] = 1
    }
    b = c(rep(1,n),pe[-n],diag(M))
    xP0 = solnp(pars=xM,
                fun=function(x) max(sum(x/xM),sum(xM/x)),
                eqfun = function(x) C%*%x,
                eqB=b,
                LB=rep(0,m),
                UB=rep(1,m),
                control=list("trace"=FALSE))
    xP = solnp(pars=xP0$pars,
               fun=function(x) max(abs(log(x)-log(xM))),
               eqfun=function(x) C%*%x,
               eqB=b,
               LB=rep(0,m),
               UB=rep(1,m),
               control=list("trace"=FALSE))
    if (xP$convergence!=0)
      warning("The optimization algorithm has not converged.")
    P = matrix(xP$pars,n,n)
  } else if (stabilizationMethod == "RM"){
    if (any(xM==0))
      stop("RM stabilization not possible unless all mutation matrix elements are positive.") 
    C = matrix(0,2*n,m)
    for (i in 1:n){
      C[seq(1,n),seq(n*(i-1)+1,n*i)] = diag(n)
      C[n+i,seq(n*(i-1)+1,n*i)] = pe
      C[2*n,seq(n*(n-1)+1,m)] = 0
      C[2*n,n*(i-1)+i] = pe[i]
    }
    b = c(rep(1,n),pe[-n],1-R)
    xP0 = solnp(pars=xM,
                fun=function(x) max(sum(x/xM),sum(xM/x)),
                eqfun = function(x) C%*%x,
                eqB=b,
                LB=as.vector((1-t)*diag(n)),
                UB=rep(1,m),
                control=list("trace"=FALSE))
    xP = solnp(pars=xP0$pars,
               fun=function(x) max(abs(log(x)-log(xM))),
               eqfun=function(x) C%*%x,
               eqB=b,
               LB=as.vector((1-t)*diag(n)),
               UB=rep(1,m),
               control=list("trace"=FALSE))
    if (xP$convergence!=0)
      warning("The optimization algorithm has not converged")  
    P = matrix(xP$pars,n,n)
  } else if (stabilizationMethod == "PM"){
    # No optimization needed here, the stabilization is unique (if it exists).
    v <- eigen(t(M))$vectors[,1]
    v <- v/sum(v)
    if (sum(diag(M)*v)<1) 
       d = R*v/((1-sum(diag(M)*v))*pe)
    else 
       d <- rep(1, n)
    if (any(d*(1-diag(M))>t)){
      stop("PM stabilization doesn't exist for these input parameters.")
    }
    P = diag(d)%*%(M-diag(n)) + diag(n)    
  } else 
    stop("Stabilization method must be either \"DP\",\"RM\" or \"PM\".")

  if (max(abs(pe%*%P-pe))>tol){
    stop("The proposed stabilization doesn't have the desired stationary distribution.")
  }
  if (max(abs(P%*%rep(1,n)-rep(1,n)))>tol){
    stop("The proposed stabilization isn't a mutation matrix.")
  }
  if (min(P)<0){
    stop("The proposed stabilization has negative elements.")
  } 

  fratio = max(max(P[M>0]/M[M>0]),max(M[M>0]/P[M>0]))
  maxrate = 1 - min(diag(P))
  return(list(stabilized=P,fratio=fratio,maxrate=maxrate))  
}

