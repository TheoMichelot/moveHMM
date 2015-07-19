
#' Simulation tool
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepFun Character string, among "gamma", "weibull", and "exp".
#' Name of the distribution from which to draw the step length values.
#' @param angleFun Character string, among "vm", and "wrpcauchy".
#' Name of the distribution from which to draw the turning angle values.
#' @param stepPar Parameters of the step length distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the function stepFun),
#' and one column for each state.
#' @param anglePar Parameters of the turning angle distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the function angleFun),
#' and one column for each state.
#' @param zeroInflation Probability of step length being zero (0 by default).
#' @param nbCov Number of covariates to simulate (0 by default).
#'
#' @return An object moveData
#' @examples
#' stepPar <- matrix(c(1,1,10,5),nrow=2) # mean1, sd1, mean2, sd2
#' anglePar <- matrix(c(0,0.5,pi,2),nrow=2) # mean1, k1, mean2, k2
#' stepFun <- "gamma"
#' angleFun <- "vm"
#' data <- simData(5,2,stepFun,angleFun,stepPar,anglePar,0.2,2)
#'
#' stepPar <- matrix(c(1,1,10,5),nrow=2) # mean1, sd1, mean2, sd2
#' anglePar <- matrix(c(0,0.5,pi,0.7),nrow=2) # mean1, k1, mean2, k2
#' stepFun <- "weibull"
#' angleFun <- "wrpcauchy"
#' data <- simData(5,2,stepFun,angleFun,stepPar,anglePar)
simData <- function(nbAnimals,nbStates,stepFun=c("gamma","weibull","exp"),
                    angleFun=c("vm","wrpcauchy"),stepPar,anglePar,zeroInflation=0,nbCov=0)
{
  stepFun <- match.arg(stepFun)
  stepFun <- paste("r",stepFun,sep="")
  angleFun <- match.arg(angleFun)
  angleFun <- paste("r",angleFun,sep="")

  data <- list()

  # generate regression parameters for transition probabilities
  beta <- matrix(rnorm((nbCov+1)*nbStates*(nbStates-1)),nrow=nbCov+1,ncol=nbStates*(nbStates-1))

  for(i in 1:nbAnimals) {
    nbObs <- sample(100:1500,1) # number of observation chosen in [100,1500]

    # initial state distribution
    delta <- rep(1,nbStates)/nbStates

    # generate covariate values
    covs <- NULL
    if(nbCov==1) covs <- rnorm(nbObs)
    if(nbCov>1) {
      covs <- rnorm(nbObs)
      for(j in 2:nbCov) covs <- cbind(covs,rnorm(nbObs))
    }

    # generate states sequence Z
    Z<-rep(NA,nbObs)
    Z[1]<-sample(1:nbStates,size=1,prob=delta)
    for (k in 2:nbObs){
      gamma <- diag(nbStates)
      g <- beta[1,]
      if(nbCov==1) g <- g + beta[2,]*covs[k]
      if(nbCov>1) {
        for(j in 2:(nbCov+1))
          g <- g + beta[j,]*covs[k,j-1]
      }
      gamma[!gamma]<-exp(g)
      gamma <- t(gamma)
      gamma<-gamma/apply(gamma,1,sum)

      Z[k]<-sample(1:nbStates,size=1,prob=gamma[Z[k-1],])
    }

    # simulate movement path X
    X <- matrix(NA,nrow=nbObs,ncol=2)
    step <- rep(NA,nbObs)
    angle <- rep(NA,nbObs)
    X[1,] <- c(0,0) # initial position of the animal
    phi <- 0
    for(k in 1:(nbObs-1)) {

      # Constitute the lists of state-dependent parameters for the step and angle
      stepArgs <- list(1); angleArgs <- list(1) # first argument = 1 (one random draw)
      if(nrow(stepPar)==1) stepArgs[[2]] <- stepPar[Z[k]]
      else {
        for(j in 1:nrow(stepPar))
          stepArgs[[j+1]] <- stepPar[j,Z[k]]
      }
      if(nrow(anglePar)==1) angleArgs[[2]] <- anglePar[Z[k]]
      else {
        for(j in 1:nrow(anglePar))
          angleArgs[[j+1]] <- anglePar[j,Z[k]]
      }
      angle[k] <- do.call(angleFun,angleArgs)
      phi<-phi+angle[k]

      # conversion between mean/sd and shape/scale if necessary
      if(stepFun=="rweibull" | stepFun=="rgamma") {
        shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
        scale <- stepArgs[[3]]^2/stepArgs[[2]]
        stepArgs[[2]] <- shape
        stepArgs[[3]] <- scale
      }
      if(runif(1)>zeroInflation)
        step[k] <- do.call(stepFun,stepArgs)
      else
        step[k] <- 0

      m <- step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
      X[k+1,] <- X[k,] + m
      angle[k] <- angle[k]-pi # angle between -pi and pi
    }
    angle[1] <- NA # the first angle value is arbitrary

    data[[i]] <- list(ID=as.character(i),x=X[,1],y=X[,2],step=step,angle=angle,covs=covs)
  }

  return(moveData(data))
}
