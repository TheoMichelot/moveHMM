
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
#'
#' @return Matrix of all probabilities.
#'
#' @examples
#' stepPar <- matrix(c(1,1,10,5),nrow=2) # mean1, sd1, mean2, sd2
#' anglePar <- matrix(c(0,0.5,pi,2),nrow=2) # mean1, k1, mean2, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar,0.2,2)
#' P <- allProbs(data[[1]],2,stepDist,angleDist,stepPar,anglePar)

allProbs <- function(data,nbStates,stepDist=c("gamma","weibull","exp"),
                     angleDist=c("NULL","vm","wrpcauchy"),stepPar,anglePar=NULL)
{
  stepDist <- match.arg(stepDist)
  stepFun <- paste("d",stepDist,sep="")
  angleDist <- match.arg(angleDist)
  if(angleDist!="NULL") angleFun <- paste("d",angleDist,sep="")

  nbObs <- length(data$step)
  allProbs <- matrix(1,nrow=nbObs,ncol=nbStates)
  stepInd <- which(!is.na(data$step))
  if(angleDist!="NULL") angleInd <- which(!is.na(data$angle))

  for(i in 1:nbStates) {
    stepProb <- rep(1,nbObs)
    angleProb <- rep(1,nbObs)

    # Constitute the lists of state-dependent parameters for the step and angle
    stepArgs <- list(data$step[stepInd])
    if(angleDist!="NULL") angleArgs <- list(data$angle[angleInd])
    if(nrow(stepPar)==1) stepArgs[[2]] <- stepPar[i]
    else {
      for(j in 1:nrow(stepPar))
        stepArgs[[j+1]] <- stepPar[j,i]
    }
    # conversion between mean/sd and shape/scale if necessary
    if(stepFun=="dweibull" | stepFun=="dgamma") {
      shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
      scale <- stepArgs[[3]]^2/stepArgs[[2]]
      stepArgs[[2]] <- shape
      if(stepFun=="dgamma") stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      else stepArgs[[3]] <- scale # dweibull expects scale
    }
    stepProb[stepInd] <- do.call(stepFun,stepArgs)

    if(angleDist!="NULL") {
      if(nrow(anglePar)==1) angleArgs[[2]] <- anglePar[i]
      else {
        for(j in 1:nrow(anglePar))
          angleArgs[[j+1]] <- anglePar[j,i]
      }

      angleProb[angleInd] <- do.call(angleFun,angleArgs)
      allProbs[,i] <- stepProb*angleProb
    }
    else allProbs[,i] <- stepProb # model step length only
  }
  return(allProbs)
}
