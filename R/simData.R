
#' Simulation tool
#'
#' Simulates movement data from an HMM.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepDist Name of the distribution of the step lengths (as a character string).
#' Supported distributions are : gamma, weibull, lnorm, exp. Default : gamma.
#' @param angleDist Name of the distribution of the turning angles (as a character string).
#' Supported distributions are : vm, wrpcauchy. Set to \code{"none"} if the angle distribution should
#' not be estimated. Default : vm.
#' @param stepPar Parameters of the step length distribution.
#' @param anglePar Parameters of the turning angle distribution.
#' @param beta Matrix of regression parameters for the transition probabilities.
#' @param covs Covariate values to include in the model, as a dataframe. The number of rows of \code{covs}
#' needs to be a multiple of the number of animals, and the same number of observations will be
#' simulated for each animal. Default : \code{NULL}. Covariates can also be simulated according to a standard
#' normal distribution, by setting \code{covs} to \code{NULL}, and specifying \code{nbCovs}.
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified of
#' \code{covs} is specified.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#' Default : \code{FALSE}. If \code{TRUE}, values for the zero-mass parameters should be
#' included in \code{stepPar}.
#' @param obsPerAnimal Either the number of the number of observations per animal (if single value),
#' or the bounds of the number of observations per animal (if vector of two values). In the latter case,
#' the numbers of obervations generated for each animal are uniformously picked from this interval.
#' Default : \code{c(500,1500)}. \code{obsPerAnimal} does not need to be specified if \code{covs} is
#' specified.
#'
#' @return An object moveData, i.e. a dataframe of :
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{step}{The step lengths}
#' \item{angle}{The turning angles (if any)}
#' \item{x}{Either easting or longitude}
#' \item{y}{Either norting or latitude}
#' \item{...}{Covariates (if any)}
#'
#' @examples
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' anglePar <- c(0,pi,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' obsPerAnimal=c(100,150)
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar,nbCovs=2,zeroInflation=TRUE,
#'                obsPerAnimal=obsPerAnimal)
#'
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(0,pi,0.5,0.7) # mean1, mean2, k1, k2
#' stepDist <- "weibull"
#' angleDist <- "wrpcauchy"
#' data <- simData(5,2,stepDist,angleDist,stepPar,anglePar,obsPerAnimal=obsPerAnimal)
#'
#' # step length only and zero-inflation
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' stepDist <- "gamma"
#' data <- simData(5,2,stepDist,"none",stepPar,nbCovs=2,zeroInflation=TRUE,obsPerAnimal=obsPerAnimal)

simData <- function(nbAnimals,nbStates,stepDist=c("gamma","weibull","lnorm","exp"),
                    angleDist=c("vm","wrpcauchy","none"),stepPar,anglePar=NULL,
                    beta=NULL,covs=NULL,nbCovs=0,zeroInflation=FALSE,obsPerAnimal=c(500,1500))
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
  if(angleDist!="none")
  {
    angleBounds <- p$bounds[(p$parSize[1]*nbStates+1):nrow(p$bounds),]
    if(length(which(anglePar<angleBounds[,1] | anglePar>angleBounds[,2]))>0)
      stop("Check the turning angle parameters bounds.")
  }

  if(length(which(obsPerAnimal<1))>0)
    stop("obsPerAnimal should have positive values.")

  if(!is.null(covs) & nbCovs>0) {
    if(ncol(covs)!=nbCovs)
      warning("covs and nbCovs argument conflicting - nbCovs was set to ncol(covs)")
  }

  if(!is.null(covs)) {
    if(nrow(covs)%%nbAnimals!=0)
      stop("The number of rows in covs should be a multiple of nbAnimals")
  }

  if(!is.null(covs)) {
    # same number of observations for all animals
    obsPerAnimal <- c(nrow(covs)/nbAnimals,nrow(covs)/nbAnimals)

    nbCovs <- ncol(covs)
  }

  if(length(obsPerAnimal)==1)
    obsPerAnimal <- rep(obsPerAnimal,2)
  else if(length(obsPerAnimal)!=2)
    stop("obsPerAnimal should be of length 1 or 2.")

  # generate regression parameters for transition probabilities
  if(is.null(beta))
    beta <- matrix(rnorm(nbStates*(nbStates-1)*(nbCovs+1)),nrow=nbCovs+1)
  # initial state distribution
  delta <- rep(1,nbStates)/nbStates

  # format parameters
  wpar <- n2w(c(stepPar,anglePar),p$bounds,beta,delta,nbStates,estAngleMean=(angleDist!="none"))
  par <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,estAngleMean=(angleDist!="none"),stationary=FALSE)

  if(zeroInflation) {
    zeroMass <- par$stepPar[nrow(par$stepPar),]
    stepPar <- par$stepPar[-(nrow(par$stepPar)),]
  }
  else {
    zeroMass <- rep(0,nbStates)
    stepPar <- par$stepPar
  }
  anglePar <- par$anglePar # i.e. NULL if angleDist=="none"

  trackData <- NULL
  allCovs <- NULL

  # build the data frame to be returned
  data <- data.frame(ID=character(),
                     step=numeric(),
                     angle=numeric(),
                     x=numeric(),
                     y=numeric())

  for (zoo in 1:nbAnimals) {
    if(obsPerAnimal[1]!=obsPerAnimal[2])
      nbObs <- sample(obsPerAnimal[1]:obsPerAnimal[2],1)
    else
      nbObs <- obsPerAnimal[1]

    # generate covariate values
    if(nbCovs>0) {
      if(is.null(covs)) {
        if(nbCovs==1) subCovs <- data.frame(cov1=rnorm(nbObs))
        if(nbCovs>1) {
          subCovs <- data.frame(cov1=rnorm(nbObs))
          for(j in 2:nbCovs) {
            c <- data.frame(rnorm(nbObs))
            colnames(c) <- paste("cov",j,sep="")
            subCovs <- cbind(subCovs,c)
          }
        }
      } else {
        # select covariate values which concern the current animal
        subCovs <- covs[((zoo-1)*obsPerAnimal[1]+1):(zoo*obsPerAnimal[1]),]
      }
      allCovs <- rbind(allCovs,subCovs)
    }

    # generate state sequence Z
    Z <- rep(NA,nbObs)
    Z[1] <- sample(1:nbStates,size=1,prob=delta)
    for (k in 2:nbObs) {
      gamma <- diag(nbStates)

      g <- beta[1,]
      if(nbCovs==1) g <- g + beta[2,]*subCovs[k,1]
      if(nbCovs>1) {
        for(j in 1:nbCovs)
          g <- g + beta[j+1,]*subCovs[k,j]
      }

      gamma[!gamma] <- exp(g)
      gamma <- t(gamma)
      gamma <- gamma/apply(gamma,1,sum)
      Z[k] <- sample(1:nbStates,size=1,prob=gamma[Z[k-1],])
    }

    X <- matrix(nbObs,nrow=nbObs,ncol=2)
    X[1,] <- c(0,0) # initial position of animal

    phi <- 0
    s <- rep(NA,nbObs)
    a <- rep(NA,nbObs)

    # simulate movement path
    for (k in 1:(nbObs-1)){
      # prepare lists of arguments for step and angle distributions
      stepArgs <- list(1) ; angleArgs <- list(1) # first argument = 1 (one random draw)
      for(j in 1:nrow(stepPar))
        stepArgs[[j+1]] <- stepPar[j,Z[k]]

      if(angleDist!="none") {
        for(j in 1:nrow(anglePar))
          angleArgs[[j+1]] <- anglePar[j,Z[k]]
      }

      if(stepDist=="gamma") {
        shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
        scale <- stepArgs[[3]]^2/stepArgs[[2]]
        stepArgs[[2]] <- shape
        stepArgs[[3]] <- 1/scale # rgamma expects rate=1/scale
      }

      if(runif(1)>zeroMass[Z[k]])
        s[k] <- do.call(stepFun,stepArgs)
      else
        s[k] <- 0

      if(angleDist!="none" & s[k]>0) {
        a[k] <- do.call(angleFun,angleArgs)
        if(a[k] >  pi) a[k] <- a[k]-2*pi
        if(a[k] < -pi) a[k] <- a[k]+2*pi
        phi <- phi + a[k]
      }
      else if(s[k]==0) {
        a[k] <- NA # angle = NA if step = 0
      }

      m <- s[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
      X[k+1,] <- X[k,] + m
    }

    a[1] <- NA # the first angle value is arbitrary
    d <- data.frame(ID=rep(zoo,nbObs),step=s,angle=a,x=X[,1],y=X[,2])
    data <- rbind(data,d)
  }

  # if covs provided as argument
  if(!is.null(covs) & is.null(allCovs))
    allCovs <- covs

  if(nbCovs>0)
    data <- cbind(data,allCovs)
  return(moveData(data))
}
