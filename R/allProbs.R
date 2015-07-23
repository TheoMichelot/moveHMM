
#' Matrix of all probabilities
#'
#' Used in the computation of the log-likelihood.
#'
#' @param data Data relative to one observed individual, as a list. Required fields :
#' step, angle.
#' @param nbStates Number of states of the HMM.
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values. Defaults to "NULL"
#' if the turning angles distributions is not estimated.
#' @param stepPar Parameters of the step length distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of stepDist),
#' and one column for each state.
#' @param anglePar Parameters of the turning angle distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of angleDist),
#' and one column for each state. Defaults to NULL if the turning angles distributions
#' is not estimated.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return Matrix of all probabilities.
#'
#' @examples
#' stepPar <- matrix(c(1,10,
#'                    1,5,
#'                    0.2,0.3),nrow=3,byrow=TRUE)
#' anglePar <- matrix(c(0,pi,0.5,2),nrow=2,byrow=TRUE)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' stepParVec <- as.vector(t(stepPar))
#' angleParVec <- as.vector(t(anglePar))
#' data <- simData(5,2,stepDist,angleDist,stepParVec,angleParVec,nbCovs=2,zeroInflation=TRUE)
#' P <- allProbs(data[[1]],2,stepDist,angleDist,stepPar,anglePar,TRUE)

allProbs <- function(data,nbStates,stepDist=c("gamma","weibull","exp"),
                     angleDist=c("NULL","vm","wrpcauchy"),stepPar,anglePar=NULL,
                     zeroInflation=FALSE)
{
  stepDist <- match.arg(stepDist)
  stepFun <- paste("d",stepDist,sep="")
  angleDist <- match.arg(angleDist)
  if(angleDist!="NULL") angleFun <- paste("d",angleDist,sep="")

  nbObs <- length(data$step)
  allProbs <- matrix(1,nrow=nbObs,ncol=nbStates)
  stepInd <- which(!is.na(data$step))
  if(angleDist!="NULL") angleInd <- which(!is.na(data$angle))

  sp <- stepPar

  for(state in 1:nbStates) {
    stepPar <- sp
    stepProb <- rep(1,nbObs)
    angleProb <- rep(1,nbObs)

    # Constitute the lists of state-dependent parameters for the step and angle
    stepArgs <- list(data$step[stepInd])
    if(angleDist!="NULL") angleArgs <- list(data$angle[angleInd])

    if(zeroInflation) {
      zeromass <- stepPar[nrow(stepPar),state]
      stepPar <- stepPar[-nrow(stepPar),]
    }
    if(nrow(stepPar)==1) stepArgs[[2]] <- stepPar[state]
    else {
      for(j in 1:nrow(stepPar))
        stepArgs[[j+1]] <- stepPar[j,state]
    }
    # conversion between mean/sd and shape/scale if necessary
    if(stepFun=="dweibull" | stepFun=="dgamma") {
      shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
      scale <- stepArgs[[3]]^2/stepArgs[[2]]
      stepArgs[[2]] <- shape
      if(stepFun=="dgamma") stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      else stepArgs[[3]] <- scale # dweibull expects scale
    }
    if(zeroInflation) {
      stepProb[stepInd] <- ifelse(data$step[stepInd]==0,
                                  zeromass, # if step==0
                                  (1-zeromass)*do.call(stepFun,stepArgs)) # if step != 0
    }
    else stepProb[stepInd] <- do.call(stepFun,stepArgs)

    if(angleDist!="NULL") {
      if(nrow(anglePar)==1) angleArgs[[2]] <- anglePar[state]
      else {
        for(j in 1:nrow(anglePar))
          angleArgs[[j+1]] <- anglePar[j,state]
      }
      angleProb[angleInd] <- do.call(angleFun,angleArgs)

      allProbs[,state] <- stepProb*angleProb
    }
    else allProbs[,state] <- stepProb # model step length only
  }

  return(allProbs)
}
