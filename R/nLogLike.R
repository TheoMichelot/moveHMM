
#' Negative log-likelihood function
#'
#' @param nbStates Number of states of the HMM.
#' @param wpar Vector of working parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in wpar. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param parSize Vector of two values : c(number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution).
#' @param data An object moveData.
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Defaults to NULL.
#'
#' @return The negative log-likelihood of wpar given data.
#'
#' @examples
#' nbStates <- 2
#' stepPar <- matrix(c(15,50,80,
#'                     10,20,30),
#'                   byrow=T,ncol=3)
#' anglePar <- matrix(c(pi,0,pi/3,
#'                      0.7,2,1),
#'                    byrow=T,ncol=3)
#' data <- simData(2,3,"gamma","vm",stepPar,anglePar,nbCov=2)
#'
#' mu0 <- c(20,80)
#' sigma0 <- c(20,40)
#' kappa0 <- c(1,1)
#' bounds <- matrix(c(0,Inf,
#'                    0,Inf,
#'                    0,Inf,
#'                    0,Inf,
#'                    0,Inf,
#'                    0,Inf),
#'                  ncol=2,byrow=T)
#' nbCovs <- ncol(data[[1]]$covs)
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=T)
#' delta0 <- c(1,1)/2
#' par0 <- c(mu0,sigma0,kappa0)
#' wpar <- n2w(par0,bounds,beta0,delta0,nbStates)
#' parSize <- c(2,1)
#' angleMean <- c(pi,0)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#'
#' # seems to work with simulated data (non-NaN result)
#' l <- nLogLike(nbStates,wpar,bounds,parSize,data,stepDist,angleDist,angleMean)
#'
#' trackData <- read.csv("~/Dropbox/Theo/real_data/two_lions.txt",sep="\t")[,c(1,10,11,12,13)]
#' data <- prepData(trackData,'euclidean')
#'
#' # does not seem to work on real data (NaN result)
#' l <- nLogLike(nbStates,wpar,bounds,parSize,data,stepDist,angleDist,angleMean)
nLogLike <- function(nbStates,wpar,bounds,parSize,data,stepDist,angleDist,angleMean=NULL)
{
  llk <- 0
  nbAnimals <- length(data)
  if(!is.null(data[[1]]$covs)) nbCovs <- ncol(data[[1]]$covs)
  else nbCovs <- 0

  par <- w2n(wpar,bounds,parSize,nbStates,nbCovs)
  if(!is.null(angleMean)) # if the turning angles' mean is not estimated
    par$anglePar <- rbind(angleMean,par$anglePar)

  for(zoo in 1:nbAnimals) {
    nbObs <- length(data[[zoo]]$step)
    covs <- data[[zoo]]$covs

    trMat <- trMatrix(nbStates,nbObs,par$beta,covs)
    allProbs <- allProbs(data[[zoo]],nbStates,stepDist,angleDist,par$stepPar,par$anglePar)

    lscale <- nLogLike_rcpp(trMat,par$delta,allProbs) # call to C++ function
    llk <- llk + lscale
  }

  return(-llk)
}
