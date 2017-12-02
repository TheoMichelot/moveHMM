
#' Matrix of all probabilities
#'
#' Used in functions \code{\link{viterbi}}, \code{\link{logAlpha}}, \code{\link{logBeta}}.
#'
#' @param data Object \code{moveData}.
#' @param nbStates Number of states of the HMM.
#' @param stepDist Name of the distribution of the step lengths.
#' @param angleDist Name of the distribution of the turning angles.
#' Set to "none" if the angle distribution should not be estimated.
#' @param stepPar Parameters of the step length distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of \code{stepDist}),
#' and one column for each state.
#' @param anglePar Parameters of the turning angle distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of \code{angleDist}),
#' and one column for each state. Default: \code{NULL} ; if the turning angles distribution
#' is not estimated.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#' Default: \code{FALSE}.
#' @param knownStates Vector of values of the state process which are known prior to fitting the
#' model (if any). Default: NULL (states are not known). This should be a vector with length the number
#' of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
#' the state is not known.
#'
#' @return Matrix of all probabilities.
#'
#' @examples
#' \dontrun{
#' stepPar <- c(1,10,1,5,0.2,0.3)
#' anglePar <- c(0,pi,0.5,2)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                  anglePar=anglePar,nbCovs=2,zeroInflation=TRUE)
#' P <- allProbs(data=data,nbStates=2,stepDist=stepDist,angleDist=angleDist,
#'                stepPar=matrix(stepPar,ncol=2,byrow=TRUE),anglePar=matrix(anglePar,ncol=2,
#'                byrow=TRUE),zeroInflation=TRUE)
#' }

allProbs <- function(data,nbStates,stepDist,angleDist,stepPar,anglePar=NULL,zeroInflation=FALSE,
                     knownStates=NULL)
{
    stepFun <- paste("d",stepDist,sep="")
    if(angleDist!="none") angleFun <- paste("d",angleDist,sep="")

    nbObs <- length(data$step)
    allProbs <- matrix(1,nrow=nbObs,ncol=nbStates)
    stepInd <- which(!is.na(data$step))
    if(angleDist!="none") angleInd <- which(!is.na(data$angle))

    sp <- stepPar

    for(state in 1:nbStates) {
        stepPar <- sp[,state]
        stepProb <- rep(1,nbObs)
        angleProb <- rep(1,nbObs)

        # Constitute the lists of state-dependent parameters for the step and angle
        stepArgs <- list(data$step[stepInd])
        if(angleDist!="none") angleArgs <- list(data$angle[angleInd])

        if(zeroInflation) {
            zeromass <- stepPar[length(stepPar)]
            stepPar <- stepPar[-length(stepPar)]
        }

        for(j in 1:length(stepPar))
            stepArgs[[j+1]] <- stepPar[j]

        # conversion between mean/sd and shape/scale if necessary
        if(stepDist=="gamma") {
            shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
            scale <- stepArgs[[3]]^2/stepArgs[[2]]
            stepArgs[[2]] <- shape
            stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        if(zeroInflation) {
            stepProb[stepInd] <- ifelse(data$step[stepInd]==0,
                                        zeromass, # if step==0
                                        (1-zeromass)*do.call(stepFun,stepArgs)) # if step != 0
        }
        else stepProb[stepInd] <- do.call(stepFun,stepArgs)

        if(angleDist!="none") {
            for(j in 1:nrow(anglePar))
                angleArgs[[j+1]] <- anglePar[j,state]

            angleProb[angleInd] <- do.call(angleFun,angleArgs)

            allProbs[,state] <- stepProb*angleProb
        }
        else allProbs[,state] <- stepProb # model step length only
    }

    # if some states are known a priori
    if(!is.null(knownStates)) {
        for(i in which(!is.na(knownStates))) {
            # set all probabilities to zero except for known state
            prob <- allProbs[i,knownStates[i]]
            allProbs[i,] <- 0
            allProbs[i,knownStates[i]] <- prob
        }
    }

    return(allProbs)
}
