
#' Simulation tool
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepPar Parameters of the step length gamma distribution. Must be provided in a
#' matrix of dimensions 2*nbStates. First row : mean ; second row : standard deviation.
#' @param anglePar Parameters of the turning angle von Mises distribution. Must be provided
#' in a matrix of dimensions 2*nbStates. First row : mean ; second row : concentration power.
#' @param nbCov Number of covariates to simulate (0 by default).
#'
#' @return An object moveData
#' @examples
#' stepPar <- matrix(c(1,1,10,5),nrow=2) # mean1, sd1, mean2, sd2
#' anglePar <- matrix(c(0,0.5,pi,2),nrow=2) # mean1, k1, mean2, k2
#' data <- simData(5,2,1,stepPar,anglePar)
simData <- function(nbAnimals,nbStates,stepPar,anglePar,nbCov=0)
{
  data <- list()

  # generate regression parameters for transition probabilities
  beta <- matrix(rnorm((nbCov+1)*nbStates*(nbStates-1)),nrow=nbCov+1,ncol=nbStates*(nbStates-1))

  for(i in 1:nbAnimals) {
    nbObs <- sample(100:150,1) # number of observation chosen in [100,150]

    # initial state distribution
    delta <- rep(1,nbStates)/nbStates

    # generate covariate values
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
      angle[k] <- rvm(1,mean=anglePar[1,Z[k]],k=anglePar[2,Z[k]])
      phi<-phi+angle[k]
      step[k]<-rgamma(1,shape=stepPar[1,Z[k]]^2/stepPar[2,Z[k]]^2,scale=stepPar[2,Z[k]]^2/stepPar[1,Z[k]])
      m<-step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
      X[k+1,]<-X[k,]+m
      angle[k] <- angle[k]-pi # angle between -pi and pi
    }
    angle[1] <- NA

    data[[i]] <- list(ID=as.character(i),x=X[,1],y=X[,2],step=step,angle=angle,covs=covs)
  }

  return(moveData(data))
}
