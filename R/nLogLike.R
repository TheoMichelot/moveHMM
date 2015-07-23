
#' Negative log-likelihood function
#'
#' @param wpar Vector of working parameters.
#' @param nbStates Number of states of the HMM.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in wpar. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param parSize Vector of two values : c(number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution).
#' @param data An object moveData.
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values. Defaults to "NULL"
#' if the turning angles distributions is not estimated.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Defaults to NULL.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return The negative log-likelihood of wpar given data.
#'
#' @examples
#' nbStates <- 2
#' stepPar <- matrix(c(15,50,80,
#'                     10,20,30),
#'                   byrow=TRUE,ncol=3)
#' anglePar <- matrix(c(pi,0,pi/3,
#'                      0.7,2,1),
#'                    byrow=TRUE,ncol=3)
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
#'                  ncol=2,byrow=TRUE)
#' nbCovs <- ncol(data[[1]]$covs)
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=TRUE)
#' delta0 <- c(1,1)/2
#' par0 <- c(mu0,sigma0,kappa0)
#' wpar <- n2w(par0,bounds,beta0,delta0,nbStates)
#' parSize <- c(2,1)
#' angleMean <- c(pi,0)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#'
#' l <- nLogLike(wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,FALSE)

nLogLike <- function(wpar,nbStates,bounds,parSize,data,stepDist=c("gamma","weibull","exp"),
                     angleDist=c("NULL","vm","wrpcauchy"),angleMean=NULL,zeroInflation=FALSE)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  if(nbStates<1) stop("nbStates must be at least 1.")
  if(length(wpar)!=sum(parSize)*nbStates+nbStates*(nbStates-1)*(nbCovs+1)+nbStates-1)
    stop("Wrong number of parameters in wpar.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data[[1]]$step)) stop("Missing field(s) in data.")

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
    allProbs <- allProbs(data[[zoo]],nbStates,stepDist,angleDist,par$stepPar,par$anglePar,
                         zeroInflation)

    lscale <- nLogLike_rcpp(trMat,par$delta,allProbs) # call to C++ function
    llk <- llk + lscale
  }

  return(-llk)
}
