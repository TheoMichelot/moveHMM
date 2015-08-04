
#' Simulation tool
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepDist Name of the distribution from which to draw the step length values.
#' @param angleDist Name of the distribution from which to draw the turning angle values.
#' @param stepPar Parameters of the step length distribution.
#' @param anglePar Parameters of the turning angle distribution.
#' @param nbCovs Number of covariates to simulate (0 by default).
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#' @param obsPerAnimal Bounds of the number of observations per animal. Default : (500,1500).
#'
#' @return An object moveData
#' @examples
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' anglePar <- c(0,pi,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar,nbCovs=2,zeroInflation=TRUE)
#'
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(0,pi,0.5,0.7) # mean1, mean2, k1, k2
#' stepDist <- "weibull"
#' angleDist <- "wrpcauchy"
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar)
#'
#' # step length only and zero-inflation
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' stepDist <- "gamma"
#' data <- simData(5,2,stepDist,"NULL",stepPar,nbCovs=2,zeroInflation=TRUE)

simData <- function(nbAnimals,nbStates,stepDist=c("gamma","weibull","exp"),
                    angleDist=c("NULL","vm","wrpcauchy"),stepPar,anglePar=NULL,
                    nbCovs=0,zeroInflation=FALSE,obsPerAnimal=c(500,1500))
{
  # check arguments
  stepDist <- match.arg(stepDist)
  stepFun <- paste("r",stepDist,sep="")
  angleDist <- match.arg(angleDist)
  angleFun <- paste("r",angleDist,sep="")

  if(nbAnimals<1) stop("nbAnimals should be at least 1.")
  if(nbStates<1) stop("nbStates should be at least 1.")
  p <- parDef(stepDist,angleDist,nbStates,TRUE,zeroInflation)

  if(length(stepPar)!=p$parSize[1]*nbStates | length(anglePar)!=p$parSize[2]*nbStates)
    stop("Wrong number of parameters")
  stepBounds <- p$bounds[1:(p$parSize[1]*nbStates),]
  if(length(which(stepPar<stepBounds[,1] | stepPar>stepBounds[,2]))>0)
    stop("Check the step length parameters bounds.")
  if(angleDist!="NULL")
  {
    angleBounds <- p$bounds[(p$parSize[1]*nbStates+1):nrow(p$bounds),]
    if(length(which(anglePar<angleBounds[,1] | anglePar>angleBounds[,2]))>0)
      stop("Check the turning angle parameters bounds.")
  }

  if(length(which(obsPerAnimal<0))>0)
    stop("obsPerAnimal should have positive values.")

  # generate regression parameters for transition probabilities
  beta <- matrix(rnorm(nbStates*(nbStates-1)*(nbCovs+1)),nrow=nbCovs+1)
  # initial state distribution
  delta <- rep(1,nbStates)/nbStates

  # format parameters
  wpar <- n2w(c(stepPar,anglePar),p$bounds,beta,delta,nbStates)
  par <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs)
  if(zeroInflation) {
    zeroMass <- par$stepPar[nrow(par$stepPar),]
    stepPar <- par$stepPar[-(nrow(par$stepPar)),]
  }
  else {
    zeroMass <- rep(0,nbStates)
    stepPar <- par$stepPar
  }
  anglePar <- par$anglePar # i.e. NULL if angleDist=="NULL"

  trackData <- NULL
  allCovs <- NULL

  for (zoo in 1:nbAnimals) {
    if(obsPerAnimal[1]!=obsPerAnimal[2])
      nbObs <- sample(obsPerAnimal[1]:obsPerAnimal[2],1)
    else
      nbObs <- obsPerAnimal[1]

    # generate covariate values
    covs <- NULL
    if(nbCovs==1) covs <- data.frame(cov1=rnorm(nbObs))
    if(nbCovs>1) {
      covs <- data.frame(cov1=rnorm(nbObs))
      for(j in 2:nbCovs) {
        c <- data.frame(rnorm(nbObs))
        colnames(c) <- paste("cov",j,sep="")
        covs <- cbind(covs,c)
      }
    }
    allCovs <- rbind(allCovs,covs)

    # generate state sequence Z
    Z <- rep(NA,nbObs)
    Z[1] <- sample(1:nbStates,size=1,prob=delta)
    for (k in 2:nbObs) {
      gamma <- diag(nbStates)

      g <- beta[1,]
      if(nbCovs==1) g <- g + beta[2,]*covs[k]
      if(nbCovs>1) {
        for(j in 1:nbCovs)
          g <- g + beta[j+1,]*covs[k,j]
      }

      gamma[!gamma] <- exp(g)
      gamma <- t(gamma)
      gamma <- gamma/apply(gamma,1,sum)
      Z[k] <- sample(1:nbStates,size=1,prob=gamma[Z[k-1],])
    }

    X <- matrix(nbObs,nrow=nbObs,ncol=2)
    X[1,] <- c(0,0) # initial position of animal

    phi <- 0
    # simulate movement path
    for (k in 1:(nbObs-1)){
      # prepare lists of arguments for step and angle distributions
      stepArgs <- list(1) ; angleArgs <- list(1) # first argument = 1 (one random draw)
      if(nrow(stepPar)==1) stepArgs[[2]] <- stepPar[Z[k]]
      else {
        for(j in 1:nrow(stepPar))
          stepArgs[[j+1]] <- stepPar[j,Z[k]]
      }
      if(angleDist!="NULL") {
        if(nrow(anglePar)==1) angleArgs[[2]] <- anglePar[Z[k]]
        else {
          for(j in 1:nrow(anglePar))
            angleArgs[[j+1]] <- anglePar[j,Z[k]]
        }
      }

      if(stepDist=="gamma") {
        shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
        scale <- stepArgs[[3]]^2/stepArgs[[2]]
        stepArgs[[2]] <- shape
        stepArgs[[3]] <- 1/scale # rgamma expects rate=1/scale
      }

      if(runif(1)>zeroMass[Z[k]]) len <- do.call(stepFun,stepArgs)
      else len <- 0

      if(angleDist!="NULL")
        phi <- phi + do.call(angleFun,angleArgs)

      m <- len*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
      X[k+1,] <- X[k,] + m
    }

    d <- data.frame(ID=rep(zoo,nbObs),x=X[,1],y=X[,2])
    trackData <- rbind(trackData,d)
  }

  # build the data frame to be returned
  data <- data.frame(ID=character(),
                     step=numeric(),
                     angle=numeric(),
                     x=numeric(),
                     y=numeric())

  for (zoo in 1:nbAnimals) {
    ind <- which(trackData$ID==unique(trackData$ID)[zoo])
    nbObs <- length(ind)

    s <- rep(NA,nbObs)
    a <- rep(NA,nbObs)
    compDir <- rep(NA,nbObs)

    for(k in 1:(nbObs-1)) {
      x <- trackData$x[ind]
      y <- trackData$y[ind]
      s[k] <- sqrt((x[k+1]-x[k])^2+(y[k+1]-y[k])^2) # euclidean distance
      coo <- c(x[k+1],y[k+1])-c(x[k],y[k])
      compDir[k] <- Arg(coo[1]+1i*coo[2]) # compass direction
      if (k!=1){
        a[k] <- compDir[k]-compDir[k-1]
        # angles between -pi and pi
        if(a[k] >  pi) a[k] <- a[k]-2*pi
        if(a[k] < -pi) a[k] <- a[k]+2*pi
      }
    }

    if(angleDist=="NULL") a <- rep(NA,nbObs)
    d <- data.frame(ID=trackData$ID[ind],step=s,angle=a,x=x,y=y)
    data <- rbind(data,d)
  }

  if(nbCovs>0) data <- cbind(data,allCovs)
  return(moveData(data))
}
