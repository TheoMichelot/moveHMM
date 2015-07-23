
#' Simulation tool
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepDist Name of the distribution from which to draw the step length values.
#' @param angleDist Name of the distribution from which to draw the turning angle values.
#' @param stepPar Parameters of the step length distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of stepDist),
#' and one column for each state.
#' @param anglePar Parameters of the turning angle distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of angleDist),
#' and one column for each state.
#' @param zeroInflation Probability of step length being zero (0 by default).
#' @param nbCovs Number of covariates to simulate (0 by default).
#'
#' @return An object moveData
#' @examples
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' anglePar <- c(0,pi,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar,0.2,2)
#'
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(0,pi,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "weibull"
#' angleDist <- "wrpcauchy"
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar)

simData <- function(nbAnimals,nbStates,stepDist=c("gamma","weibull","exp"),
                    angleDist=c("vm","wrpcauchy"),stepPar,anglePar,nbCovs=0)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  stepFun <- paste("r",stepDist,sep="")
  angleDist <- match.arg(angleDist)
  angleFun <- paste("r",angleDist,sep="")

  if(nbAnimals<1) stop("nbAnimals should be at least 1.")
  if(nbStates<1) stop("nbStates should be at least 1.")
  p <- parDef(stepDist,angleDist,nbStates,TRUE,TRUE)
  if(length(stepPar)!=p$parSize[1]*nbStates | length(anglePar)!=p$parSize[2]*nbStates)
    stop("Wrong number of parameters")
  stepBounds <- p$bounds[1:(p$parSize[1]*nbStates),]
  angleBounds <- p$bounds[(p$parSize[1]*nbStates+1):nrow(p$bounds),]
  if(length(which(stepPar<stepBounds[,1] | stepPar>stepBounds[,2]))>0 |
       length(which(anglePar<angleBounds[,1] | anglePar>angleBounds[,2]))>0)
    stop("Check the parameters bounds.")

  # generate regression parameters for transition probabilities
  beta <- matrix(rnorm((nbCovs+1)*nbStates*(nbStates-1)),nrow=nbCovs+1,ncol=nbStates*(nbStates-1))
  # initial state distribution
  delta <- rep(1,nbStates)/nbStates

  # format the parameters
  wpar <- n2w(c(stepPar,anglePar),p$bounds,beta,delta,nbStates)
  par <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs)
  zeroInflation <- par$stepPar[nrow(par$stepPar),]
  stepPar <- par$stepPar[-(nrow(par$stepPar)),]
  anglePar <- par$anglePar

  data <- list()
  for(i in 1:nbAnimals) {
    nbObs <- sample(1000:1500,1) # number of observations chosen in [1000,1500]

    # generate covariate values
    covs <- NULL
    if(nbCovs==1) covs <- rnorm(nbObs)
    if(nbCovs>1) {
      covs <- rnorm(nbObs)
      for(j in 2:nbCovs) covs <- cbind(covs,rnorm(nbObs))
    }

    # generate states sequence Z
    Z<-rep(NA,nbObs)
    Z[1]<-sample(1:nbStates,size=1,prob=delta)
    for (k in 2:nbObs){
      gamma <- diag(nbStates)
      g <- beta[1,]
      if(nbCovs==1) g <- g + beta[2,]*covs[k]
      if(nbCovs>1) {
        for(j in 2:(nbCovs+1))
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
      angle[k] <- do.call(angleFun,angleArgs)-pi # angle between -pi and pi
      phi<-phi+angle[k]

      # conversion between mean/sd and shape/scale if necessary
      if(stepFun=="rweibull" | stepFun=="rgamma") {
        shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
        scale <- stepArgs[[3]]^2/stepArgs[[2]]
        stepArgs[[2]] <- shape
        if(stepFun=="rgamma") stepArgs[[3]] <- 1/scale # rgamma expects rate=1/scale
        else stepArgs[[3]] <- scale # rweibull expects scale
      }

      if(runif(1)>zeroInflation[Z[k]])
        step[k] <- do.call(stepFun,stepArgs)
      else
        step[k] <- 0

      m <- step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
      X[k+1,] <- X[k,] + m
    }
    angle[1] <- NA # the first angle value is arbitrary

    data[[i]] <- list(ID=as.character(i),x=X[,1],y=X[,2],step=step,angle=angle,covs=covs)
  }

  return(moveData(data))
}
